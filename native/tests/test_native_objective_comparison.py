from pathlib import Path
import numpy as np
import pandas as pd
import pytest
from scipy.special import logsumexp

from spea.native_objective_comparison import (
    encode_objective, save_common_borders, fit_native_objective,
    block_observation_mass, held_diagnostics, area_capture,
    select_complete_candidates)


def example():
    q = np.arange(8,dtype=np.int64)
    p = q[[1,2,6]]
    area = np.array([1,2,4,3,2,5,1,2.],float)
    frame = pd.DataFrame({'x':q.astype('float32'), 'z':np.sin(q).astype('float32')})
    return p,q,area,frame


@pytest.mark.parametrize('objective',['bce','poisson'])
def test_actual_first_newton_leaf_damping(objective,tmp_path):
    p,q,area,frame = example()
    encoding = encode_objective(objective,p,q,area,.03)
    shared = pd.concat([frame.iloc[p],frame],ignore_index=True)
    borders = tmp_path/'borders.tsv'
    save_common_borders(shared,borders)
    fit_frame = frame.iloc[encoding.native].reset_index(drop=True)
    candidate = dict(iterations=1,depth=1,learning_rate=.1)
    model,audit = fit_native_objective(fit_frame,encoding,objective,candidate,borders,4,
                                     leaf_iterations=1,threads=1)
    leaves = model.calc_leaf_indexes(fit_frame)[:,0]
    expected = np.empty(len(fit_frame))
    for leaf in np.unique(leaves):
        mask = leaves == leaf
        w,y = encoding.weight[mask],encoding.target[mask]
        if objective == 'bce':
            numerator = np.dot(w,y-.5); curvature = .25*w.sum()
        else:
            numerator = np.dot(w,y-1.); curvature = w.sum()
        expected[mask] = .1*numerator/(curvature+encoding.effective_leaf_damping)
    actual = model.predict(fit_frame,prediction_type='RawFormulaVal')
    np.testing.assert_allclose(actual,expected,atol=2e-7,rtol=2e-7)
    assert audit['training_metric_error'] < 1e-6
    assert encoding.native_l2_leaf_reg*encoding.weight.mean() == pytest.approx(.09)


def test_poisson_is_additive_and_background_area_units_cancel():
    p,q,area,_ = example();eta=np.linspace(-.5,.5,len(q))
    first=encode_objective('poisson',p,q,area,.02)
    second=encode_objective('poisson',p,q,area*1e6,.02)
    np.testing.assert_allclose(first.target,second.target,rtol=1e-15)
    np.testing.assert_allclose(first.weight,second.weight,rtol=1e-15)
    native=np.dot(first.weight,np.exp(eta)-first.target*eta)
    additive=np.dot(first.reference_on_q,np.exp(eta))-first.empirical_on_q@eta
    assert native == pytest.approx(additive)
    assert first.empirical_on_q.sum() == len(p)
    assert first.reference_on_q.sum() == pytest.approx(len(p))


@pytest.mark.parametrize('objective',['bce','poisson'])
def test_native_links_serialization_and_nan_retention(objective,tmp_path):
    from catboost import CatBoostClassifier,CatBoostRegressor
    p,q,area,frame = example();frame.loc[3,'z']=np.nan
    encoding=encode_objective(objective,p,q,area,.02)
    borders=tmp_path/'borders.tsv';save_common_borders(pd.concat([frame.iloc[p],frame]),borders)
    model,audit=fit_native_objective(frame.iloc[encoding.native].reset_index(drop=True),
        encoding,objective,dict(iterations=12,depth=2,learning_rate=.05),borders,8,threads=1)
    raw=model.predict(frame,prediction_type='RawFormulaVal')
    assert np.isfinite(raw).all()
    if objective=='poisson':
        np.testing.assert_allclose(model.predict(frame),np.exp(raw),rtol=1e-12)
    else:
        from scipy.special import expit
        np.testing.assert_allclose(model.predict_proba(frame)[:,1],expit(raw),rtol=1e-12)
    path=tmp_path/'model.cbm';model.save_model(str(path))
    reload=(CatBoostClassifier() if objective=='bce' else CatBoostRegressor()).load_model(str(path))
    np.testing.assert_array_equal(raw,reload.predict(frame,prediction_type='RawFormulaVal'))
    assert audit['tree_count']==12 and not audit['unpenalized_intercept_profiled']


def test_missing_occupied_support_and_duplicate_area_fail():
    with pytest.raises(ValueError,match='occupied cell'):
        encode_objective('poisson',[0],[1,2],[1,1],.01)
    with pytest.raises(ValueError,match='sorted positive quadrature'):
        encode_objective('poisson',[1],[1,1],[1,1],.01)


def test_uniform_diagnostics_and_fractional_ties():
    eta=np.zeros(4);pe=np.zeros(3);area=np.array([1.,5.,2.,3.]);blocks=[0,0,1]
    result,scores,captures=held_diagnostics(pe,eta,area,blocks)
    assert result['mean_block_log_score']==pytest.approx(0.,abs=1e-14)
    assert result['area_weighted_source_auc']==pytest.approx(.5)
    assert result['source_average_precision_prior_half']==pytest.approx(.5)
    assert result['source_logloss_prior_half']==pytest.approx(np.log(2.))
    assert result['source_brier_prior_half']==pytest.approx(.25)
    assert result['entropy_effective_area_fraction']==pytest.approx(1.)
    for c in captures:
        assert c['held_presence_capture']==pytest.approx(c['area_fraction'])
        assert c['selected_area_fraction']==pytest.approx(c['area_fraction'])


def test_ranking_does_not_capture_density_overconfidence():
    eta=np.array([-2.,0.,1.,3.]);area=np.array([1.,2.,1.,2.]);pe=eta[[0,3]]
    a,_,_=held_diagnostics(pe,eta,area,[0,1])
    b,_,_=held_diagnostics(pe*10,eta*10,area,[0,1])
    assert a['area_weighted_source_auc']==pytest.approx(b['area_weighted_source_auc'])
    assert a['capture_top_10pct_area']==pytest.approx(b['capture_top_10pct_area'])
    assert b['mean_block_log_score'] < a['mean_block_log_score']-10


def test_auc_matches_independent_area_pair_ranking_and_block_mass():
    eta=np.array([0.,1.,1.,3.]);area=np.array([1.,2.,3.,4.]);pe=eta[[1,3,3]]
    blocks=np.array([0,1,1]);a=block_observation_mass(blocks)
    result,_,_=held_diagnostics(pe,eta,area,blocks)
    oracle=sum(ai*np.sum(area*((e>eta)+.5*(e==eta)))/area.sum() for ai,e in zip(a,pe))
    assert result['area_weighted_source_auc']==pytest.approx(oracle)
    np.testing.assert_allclose(a,[.5,.25,.25])


def test_partial_inner_failure_cannot_win_selection():
    matrix=np.array([[100.,np.nan],[.2,.3],[-1.,-.5]])
    chosen,reason=select_complete_candidates(matrix,['partial','complete','weak'])
    assert chosen=='complete'
    assert select_complete_candidates(np.array([[10.],[12.]]),['baseline','other'])[0]=='baseline'
    with pytest.raises(ValueError,match='No candidate completed'):
        select_complete_candidates(np.full((2,2),np.nan),['a','b'])
