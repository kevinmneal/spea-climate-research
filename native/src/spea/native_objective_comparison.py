"""Matched native CatBoost objectives and common held-landscape diagnostics.

These are pooled relative-recording-density fits, not occupancy models or
separately conditioned regional likelihoods. All cutoffs are training-only.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.special import expit, logsumexp

from .genus_tuning import select_candidate
from .native_benchmark import relative_scores, block_score
from .objective_measures import binary_class_weights, poisson_quadrature_encoding


@dataclass
class ObjectiveEncoding:
    native: np.ndarray
    target: np.ndarray
    weight: np.ndarray
    empirical_on_q: np.ndarray
    reference_on_q: np.ndarray
    presence_count: int
    leaf_damping_per_presence: float
    native_l2_leaf_reg: float
    effective_leaf_damping: float


def encode_objective(objective, presence_native, q_native, q_area,
                     leaf_damping_per_presence):
    """Fix event mass to n_presence, and integration mass to the same total.

    CatBoost scales its native leaf penalty by mean object weight. Supplying
    lambda/mean(weight) therefore fixes effective damping to r*n_presence.
    This convention is a leaf-update control, not equivalent coefficient priors
    or proof of identical optimization paths across the two losses.
    """
    p = np.asarray(presence_native, dtype=np.int64)
    q = np.asarray(q_native, dtype=np.int64)
    area = np.asarray(q_area, dtype=float)
    r = float(leaf_damping_per_presence)
    if (p.ndim != 1 or q.ndim != 1 or not len(p) or not len(q)
            or len(np.unique(p)) != len(p) or np.any(np.diff(q) <= 0)
            or area.shape != q.shape or not np.isfinite(area).all()
            or np.any(area <= 0) or not np.isfinite(r) or r < 0):
        raise ValueError('Unique occupied cells, sorted positive quadrature and valid damping required')
    positions = np.searchsorted(q, p)
    if np.any(positions == len(q)) or not np.array_equal(q[positions], p):
        raise ValueError('Every occupied cell must occur on positive quadrature support')
    reference = area / area.max()
    reference *= len(p) / reference.sum()
    empirical = np.zeros(len(q), dtype=float)
    empirical[positions] = 1.
    if objective == 'bce':
        native = np.r_[p, q]
        target = np.r_[np.ones(len(p)), np.zeros(len(q))]
        weight = binary_class_weights(np.ones(len(p)), area, loss_mass=len(p))
    elif objective == 'poisson':
        native = q.copy()
        target, weight = poisson_quadrature_encoding(empirical, reference)
    else:
        raise ValueError('Expected bce or poisson')
    damping = r * len(p)
    native_l2 = damping / weight.mean()
    return ObjectiveEncoding(native, target, weight, empirical, reference,
                             len(p), r, float(native_l2), float(damping))


def save_common_borders(training_frame, path, *, border_count=128, seed=0):
    """Freeze numeric cuts from the shared p+q training table, without labels.

    Both losses receive the same file. Its source is restricted to training
    geography, including training-only NaNs handled with the existing Min rule.
    """
    from catboost import Pool
    path = Path(path)
    if path.exists():
        raise FileExistsError(path)
    pool = Pool(training_frame)
    pool.quantize(border_count=border_count, feature_border_type='GreedyLogSum',
                  nan_mode='Min', random_seed=int(seed), task_type='CPU')
    pool.save_quantization_borders(str(path))


def fit_native_objective(frame, encoding, objective, candidate, borders, seed,
                         *, callbacks=None, leaf_iterations=10, threads=2):
    """Fit fixed-budget native loss; no held rows or statistical early stop.

    The intercept begins at zero for both normalized measures. It can change
    through tree leaves; we do not claim exact profiling of an unpenalized bias.
    """
    from catboost import CatBoostClassifier, CatBoostRegressor
    if len(frame) != len(encoding.native):
        raise ValueError('Feature and objective rows must align')
    if not np.array_equal(np.isfinite(encoding.target), np.ones(len(frame), bool)):
        raise ValueError('Finite native targets required')
    parameters = dict(
        loss_function='Logloss' if objective == 'bce' else 'Poisson',
        iterations=int(candidate['iterations']), depth=int(candidate['depth']),
        learning_rate=float(candidate['learning_rate']),
        l2_leaf_reg=encoding.native_l2_leaf_reg, border_count=128,
        feature_border_type='GreedyLogSum', input_borders=str(borders),
        bootstrap_type='No', random_strength=0., boosting_type='Plain',
        leaf_estimation_method='Newton', leaf_estimation_iterations=int(leaf_iterations),
        leaf_estimation_backtracking='AnyImprovement' if leaf_iterations > 1 else 'No',
        boost_from_average=False, one_hot_max_size=255, nan_mode='Min',
        random_seed=int(seed), thread_count=int(threads), task_type='CPU',
        allow_writing_files=False, verbose=False, use_best_model=False)
    constructor = CatBoostClassifier if objective == 'bce' else CatBoostRegressor
    model = constructor(**parameters)
    model.fit(frame, encoding.target, sample_weight=encoding.weight,
              callbacks=callbacks)
    if model.tree_count_ != candidate['iterations']:
        raise ValueError('Unexpected truncated fit')
    eta = np.asarray(model.predict(frame, prediction_type='RawFormulaVal',
                                   thread_count=threads), dtype=float)
    if not np.isfinite(eta).all():
        raise ValueError('Nonfinite training margins')
    if objective == 'bce':
        terms = np.logaddexp(0., eta) - encoding.target * eta
    else:
        with np.errstate(over='raise', invalid='raise'):
            terms = np.exp(eta) - encoding.target * eta
    loss = float(np.dot(encoding.weight, terms) / encoding.weight.sum())
    from catboost.utils import eval_metric
    library_loss = float(eval_metric(encoding.target, eta, parameters['loss_function'],
                                    weight=encoding.weight, thread_count=threads)[0])
    history_loss = float(next(iter(model.get_evals_result()['learn'].values()))[-1])
    # The engine's incremental learn history can differ from its metric freshly
    # recomputed on serialized-model margins. Verify the latter and retain both;
    # do not relax the identity to accommodate the internal-history discrepancy.
    if not np.isclose(loss, library_loss, rtol=2e-6, atol=2e-6):
        raise ValueError('Native training metric differs from explicit weighted loss')
    audit = dict(objective=objective, presence_mass=encoding.presence_count,
                 object_rows=len(frame), total_object_weight=float(encoding.weight.sum()),
                 native_l2_leaf_reg=encoding.native_l2_leaf_reg,
                 effective_leaf_damping=encoding.effective_leaf_damping,
                 leaf_damping_per_presence=encoding.leaf_damping_per_presence,
                 maximum_target=float(encoding.target.max()),
                 training_loss=loss, library_training_loss=library_loss,
                 incremental_training_history_loss=history_loss,
                 history_recomputed_difference=history_loss-library_loss,
                 training_metric_error=abs(loss-library_loss),
                 raw_margin_minimum=float(eta.min()), raw_margin_maximum=float(eta.max()),
                 tree_count=int(model.tree_count_), parameters=model.get_all_params(),
                 no_held_data_in_fit=True, statistical_early_stopping=False,
                 unpenalized_intercept_profiled=False)
    return model, audit


def block_observation_mass(blocks):
    blocks = np.asarray(blocks)
    if blocks.ndim != 1 or not len(blocks):
        raise ValueError('Nonempty aligned block IDs required')
    unique, inverse, counts = np.unique(blocks, return_inverse=True, return_counts=True)
    return 1. / (len(unique) * counts[inverse])


def area_capture(presence_eta, target_eta, area, presence_mass, fractions=(.05, .1, .2)):
    """Cut by physical area, fractionally admitting all exact boundary ties.

    Every cell with the boundary score receives the same inclusion fraction.
    This avoids arbitrary native-ID ties and preserves exact requested area.
    """
    pe, eta, area, mass = (np.asarray(x, float) for x in
                           (presence_eta, target_eta, area, presence_mass))
    if (pe.shape != mass.shape or eta.shape != area.shape or not len(pe)
            or not len(eta) or not all(np.isfinite(x).all() for x in (pe,eta,area,mass))
            or np.any(area <= 0) or np.any(mass < 0) or mass.sum() <= 0):
        raise ValueError('Finite aligned scores, positive area and observation mass required')
    order = np.argsort(-eta, kind='stable')
    cumulative = np.cumsum(area[order])
    total = float(area.sum()); result = []
    for fraction in fractions:
        if not 0 < fraction < 1:
            raise ValueError('Landscape fractions must lie between zero and one')
        index = min(int(np.searchsorted(cumulative, fraction*total)), len(eta)-1)
        boundary = float(eta[order[index]])
        higher, tied = eta > boundary, eta == boundary
        admit = float(np.clip((fraction*total-area[higher].sum())/area[tied].sum(),0.,1.))
        chosen = higher.astype(float) + admit*tied
        capture = float(np.dot(mass, (pe>boundary)+admit*(pe==boundary))/mass.sum())
        result.append(dict(area_fraction=float(fraction), boundary_raw_eta=boundary,
                           boundary_tie_fraction=admit, held_presence_capture=capture,
                           selected_area_fraction=float(np.dot(area,chosen)/total)))
    return result


def held_diagnostics(presence_eta, target_eta, area, blocks, novel_presence=None):
    """One common density score plus explicitly artificial source diagnostics."""
    from sklearn.metrics import roc_auc_score, average_precision_score
    pe, eta, area = (np.asarray(x, float) for x in (presence_eta,target_eta,area))
    mass = block_observation_mass(blocks)
    if pe.shape != mass.shape or eta.shape != area.shape or np.any(area <= 0):
        raise ValueError('Held observation and landscape arrays must align')
    scores, norm = relative_scores(pe, eta, area)
    relative = eta - norm
    if not np.isfinite(scores).all() or not np.isfinite(area).all():
        raise ValueError('Finite held observations and area required')
    availability = area / area.sum()
    y = np.r_[np.ones(len(pe)), np.zeros(len(eta))]
    joined_score = np.r_[scores,relative]
    source_weight = np.r_[mass,availability]
    probability = expit(joined_score)
    source_logloss = float(np.dot(source_weight,
                                 np.logaddexp(0.,joined_score)-y*joined_score)/2.)
    brier = float(np.dot(source_weight, (probability-y)**2)/2.)
    density_mass = np.exp(relative+np.log(availability))
    kl = float(np.dot(density_mass, relative))
    captures = area_capture(pe,eta,area,mass)
    result = dict(mean_block_log_score=block_score(scores,blocks),
                  mean_cell_log_score=float(scores.mean()),
                  area_weighted_source_auc=float(roc_auc_score(y,joined_score,sample_weight=source_weight)),
                  source_average_precision_prior_half=float(average_precision_score(y,joined_score,sample_weight=source_weight)),
                  source_logloss_prior_half=source_logloss,
                  source_brier_prior_half=brier,
                  predicted_mass_sum=float(density_mass.sum()),
                  density_kl_from_uniform_area=kl,
                  entropy_effective_area_fraction=float(np.exp(-max(kl,0.))),
                  target_log_relative_minimum=float(relative.min()),
                  target_log_relative_maximum=float(relative.max()),
                  target_log_normalizer=norm)
    for c in captures:
        result[f"capture_top_{int(round(100*c['area_fraction']))}pct_area"] = c['held_presence_capture']
    if novel_presence is not None:
        novel = np.asarray(novel_presence,bool)
        if novel.shape != pe.shape:
            raise ValueError('Novelty mask must align with held observations')
        result['novel_presence_fraction'] = float(novel.mean())
        for name, mask in [('within_training_marginals',~novel),('outside_training_marginals',novel)]:
            result[name+'_cells'] = int(mask.sum())
            result[name+'_block_score'] = (block_score(scores[mask],np.asarray(blocks)[mask])
                                           if mask.any() else None)
    return result, scores, captures


def select_complete_candidates(scores, candidate_ids, minimum_folds=2):
    """A candidate failing any eligible fold is excluded from selection."""
    values = np.asarray(scores,float)
    if values.ndim != 2 or values.shape[0] != len(candidate_ids):
        raise ValueError('Candidate by eligible-inner-fold scores required')
    if values.shape[1] < minimum_folds:
        return candidate_ids[0], 'insufficient_inner_folds_keep_baseline'
    complete = np.isfinite(values).all(axis=1)
    if not complete.any():
        raise ValueError('No candidate completed every eligible inner fold')
    ids = [c for c, use in zip(candidate_ids,complete) if use]
    selected, reason = select_candidate(values[complete],ids,minimum_folds)
    return selected, reason
