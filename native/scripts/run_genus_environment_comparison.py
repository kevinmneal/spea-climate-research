#!/usr/bin/env python3
"""Fixed environmental-block comparison with ordinary weighted CatBoost Logloss.

No empirical work occurs on import or --plan-only. All training rows, class
weights, seeds, feature groups and full-map choices are fixed before fitting.
"""
from __future__ import annotations
import os
for variable in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:
    os.environ.setdefault(variable,'1')
import argparse,json,resource,sys,time,traceback,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import logsumexp
import catboost
from catboost import CatBoostClassifier
from spea.provenance import ROOT,sha256,write_json,utcnow
from spea.environment_benchmark import NativeEnvironment,case_control_weights
from spea.genus_model import geographic_blocks,stratified_quadrature
from spea.native_benchmark import relative_scores,block_score


class BudgetExceeded(RuntimeError):pass


class Budget:
    def __init__(self,cfg,prior_cpu=0.):
        self.started=time.process_time();self.prior_cpu=prior_cpu
        self.limit_seconds=60*cfg['cpu_minutes'];self.limit_memory=cfg['peak_gib']*1024**3
        self.callback_reason=None
    def cpu(self):return time.process_time()-self.started+self.prior_cpu
    def check(self):
        if self.cpu()>=self.limit_seconds:raise BudgetExceeded('Predeclared CPU budget exhausted')
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>=self.limit_memory:raise BudgetExceeded('Predeclared resident-memory budget exhausted')
    def after_iteration(self,info):
        # This is only a hard resource guard, never statistical early stopping.
        try:self.check();return True
        except BudgetExceeded as error:self.callback_reason=str(error);return False


def digest_array(values):
    import hashlib
    a=np.ascontiguousarray(values)
    return hashlib.sha256(a.dtype.str.encode()+str(a.shape).encode()+a.tobytes()).hexdigest()


def training_reference(frame,n_presence,keys):
    """Marginal bounds use training quadrature only; no held covariates."""
    q=frame.iloc[n_presence:];numeric=[k for k in keys if k!='cec_level2'];lower=[];upper=[]
    for key in numeric:
        x=q[key].to_numpy(float);finite=x[np.isfinite(x)]
        lower.append(float(finite.min()) if len(finite) else np.nan)
        upper.append(float(finite.max()) if len(finite) else np.nan)
    seen=sorted(q.cec_level2.unique().tolist()) if 'cec_level2' in q else []
    return {'numeric_keys':numeric,'minimum':np.asarray(lower),'maximum':np.asarray(upper),'seen_cec':seen}


def novelty(frame,reference):
    x=frame[reference['numeric_keys']].to_numpy(float);finite=np.isfinite(x)
    supported=np.isfinite(reference['minimum'])&np.isfinite(reference['maximum'])
    outside=finite&((x<reference['minimum'])|(x>reference['maximum'])|~supported)
    count=outside.sum(axis=1).astype('int16');missing=(~finite).sum(axis=1).astype('int16')
    unseen=~frame.cec_level2.isin(reference['seen_cec']).to_numpy() if 'cec_level2' in frame else np.zeros(len(frame),bool)
    return count>0,count,missing,unseen


def predict(model,ids,keys,reference,reader,chunk,budget):
    n=len(ids);eta=np.empty(n,'float64');novel=np.empty(n,bool);counts=np.empty(n,'int16');missing=np.empty(n,'int16');unseen=np.empty(n,bool)
    needed=keys+(['cec_level2'] if 'cec_level2' in reader.fields and 'cec_level2' not in keys else [])
    for start in range(0,n,chunk):
        budget.check();stop=min(start+chunk,n);frame=reader.read(ids[start:stop],needed)
        eta[start:stop]=model.predict(frame[keys],prediction_type='RawFormulaVal',thread_count=2)
        novel[start:stop],counts[start:stop],missing[start:stop],unseen[start:stop]=novelty(frame,reference)
    if not np.isfinite(eta).all():raise ValueError('Nonfinite raw classifier log odds')
    return {'raw_eta':eta,'novel_numeric':novel,'numeric_novel_count':counts,'numeric_missing_count':missing,'unseen_cec':unseen}


def fit_model(train_frame,n_presence,weights,keys,cfg,seed,budget):
    budget.check();budget.callback_reason=None
    labels=np.r_[np.ones(n_presence,dtype='int8'),np.zeros(len(train_frame)-n_presence,dtype='int8')]
    model=CatBoostClassifier(**cfg['parameters'],random_seed=int(seed))
    categories=[k for k in keys if k=='cec_level2']
    # Do not pass held data, eval_set, custom losses or stopping-round parameters.
    model.fit(train_frame[keys],labels,sample_weight=weights,cat_features=categories,callbacks=[budget])
    if budget.callback_reason:raise BudgetExceeded(budget.callback_reason)
    budget.check()
    if model.tree_count_!=cfg['parameters']['iterations']:raise ValueError('Unexpected truncated model')
    return model


def save_model(model,stem,out,keys,reference,weights,seed):
    model.save_model(str(out/'models'/(stem+'.cbm')))
    model.save_borders(str(out/'models'/(stem+'_borders.tsv')))
    write_json(out/'models'/(stem+'_metadata.json'),{
        'feature_names':keys,'categorical_features':[k for k in keys if k=='cec_level2'],
        'parameters':model.get_all_params(),'tree_count':int(model.tree_count_),'seed':int(seed),
        'class_object_weights_sha256':digest_array(weights),'background_numeric_keys':reference['numeric_keys'],
        'background_minimum':[float(x) if np.isfinite(x) else None for x in reference['minimum']],
        'background_maximum':[float(x) if np.isfinite(x) else None for x in reference['maximum']],
        'training_background_cec_categories':reference['seen_cec'],
        'novelty_definition':'Outside sampled training-quadrature marginal numeric ranges; finite values with no finite training support also flagged. Missing values counted separately. CEC unseen diagnostics do not enter non-CEC models.',
        'target':'ordinary weighted presence/background classifier raw log odds; relative recorded-cell density, not occupancy',
        'no_held_data_in_fit':True,'statistical_early_stopping':False})


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--run-id',required=True);ap.add_argument('--config',default='config/genus_environment_comparison.json');ap.add_argument('--plan-only',action='store_true')
    mode=ap.add_mutually_exclusive_group();mode.add_argument('--cv-only',action='store_true');mode.add_argument('--full-from-run',help='New4-fit projection run using a completed100-fit CV run; shared104-fit/CPU budget')
    args=ap.parse_args()
    if Path(args.run_id).name!=args.run_id:raise ValueError('Run ID must be one path component')
    cfgpath=ROOT/args.config;cfg=json.loads(cfgpath.read_text());expected_cv=len(cfg['species'])*len(cfg['folds'])*len(cfg['families'])
    if expected_cv!=100 or cfg['budget']['maximum_fits']!=104:raise ValueError('Unexpected frozen design size')
    effective={**cfg,'execution_mode':'full_only' if args.full_from_run else 'cv_only' if args.cv_only else 'cv_and_full',
      'cv_seed_formula':'seed+1000*species_index+fold','full_data':{'family':'broad_soil','presence_rule':'All canonical preparation cells whose point-native cell lies in unchanged primaryM; outsideM retained in membership audit',
      'quadrature_nonpresence_per_25km_block':16,'quadrature_stratum_width_m':25000,'seed_formula':'seed+1000*species_index+100',
      'projection_grid':'Every primary native cell from preparation; exact saved per-row geodesic areas; no coarse downsampling',
      'display_normalization':'Subtract log(area-weighted mean exp(raw classifier log odds)) across whole primaryM',
      'fixed_before_results':True}}
    if args.plan_only:print(json.dumps(effective,indent=2));return
    if catboost.__version__!='1.2.10':raise ValueError('Pinned CatBoost1.2.10 required')
    prep=ROOT/'runs'/cfg['preparation_run'];environment=ROOT/'runs'/cfg['environment_run'];env_path=environment/'registry.json'
    if not env_path.exists():raise FileNotFoundError('Environmental feature registry not ready; do not fit partial families')
    cfg_digest=sha256(cfgpath);prior_cpu=0.;prior_fits=0;prior=None
    if args.full_from_run:
        prior=ROOT/'runs'/args.full_from_run;prior_record=json.loads((prior/'run.json').read_text())
        if prior_record['status']!='completed' or prior_record['cv_fit_attempts']!=100 or prior_record['full_data_fit_attempts']!=0:raise ValueError('Full-data continuation requires completed100-fit CV-only run')
        if prior_record['input_hashes'].get(args.config)!=cfg_digest or prior_record['input_hashes'].get(str(env_path.relative_to(ROOT)))!=sha256(env_path):raise ValueError('CV source/config differs')
        prior_cpu=prior_record['cpu_seconds'];prior_fits=prior_record['fit_attempts']
    out=ROOT/'runs'/args.run_id;out.mkdir(exist_ok=False)
    for folder in ['models','predictions','projection','training_design']:(out/folder).mkdir()
    start=time.monotonic();budget=Budget(cfg['budget'],prior_cpu);inputs={};checks=[];fits=[];cell_scores=[];projection_membership=[]
    metrics=[dict(species=s,fold=f,family=k,status='not_run') for s in cfg['species'] for f in cfg['folds'] for k in cfg['families']] if not args.full_from_run else []
    projections=[dict(species=s,family='broad_soil',status='not_requested' if args.cv_only else 'not_run') for s in cfg['species']]
    record={'status':'running','created_utc':utcnow(),'command_argv':sys.argv,'config':cfg,'effective_config':effective,'input_hashes':inputs,'fit_attempts':0,'cv_fit_attempts':0,'full_data_fit_attempts':0,'prior_fit_attempts':prior_fits,'prior_cpu_seconds':prior_cpu,'formal_significance_tests':0,'parameter_searches':0,'occurrence_changes':0}
    write_json(out/'effective_config.json',effective);write_json(out/'run.json',record)
    def track(path,expected=None):
        digest=sha256(path)
        if expected is not None and digest!=expected:raise ValueError('Changed input:'+str(path))
        inputs[str(path.relative_to(ROOT))]=digest;return path
    def check(name,value):
        checks.append(dict(name=name,passed=bool(value)))
        if not value:raise ValueError(name)
    def outputs_hash():return {str(p.relative_to(out)):sha256(p) for p in out.rglob('*') if p.is_file() and p.name!='run.json'}
    def flush():
        for name,rows in [('metrics.csv',metrics),('cell_scores.csv.gz',cell_scores),('fit_audit.csv',fits),('projection_summary.csv',projections),('projection_training_membership.csv.gz',projection_membership)]:
            if rows:pd.DataFrame(rows).to_csv(out/name,index=False)
        record.update(elapsed_seconds=time.monotonic()-start,cpu_seconds=budget.cpu(),new_cpu_seconds=budget.cpu()-prior_cpu,peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write_json(out/'run.json',record)
    def begin_attempt(audit):
        budget.check()
        if prior_fits+record['fit_attempts']>=cfg['budget']['maximum_fits']:raise BudgetExceeded('Predeclared fit-count budget exhausted')
        record['fit_attempts']+=1;record['cv_fit_attempts' if audit['kind']=='cv' else 'full_data_fit_attempts']+=1
        audit['status']='running';fits.append(audit);flush()
    def fit_audit_failure(audit,error):audit.update(status='failed',error=str(error),traceback=traceback.format_exc())
    try:
        pm=json.loads(track(prep/'run.json').read_text());em=json.loads(track(environment/'run.json').read_text())
        if pm['status']!='completed' or em['status']!='completed':raise ValueError('Incomplete immutable preparation')
        for name,h in pm['output_hashes'].items():track(prep/name,h)
        env=json.loads(track(env_path,em['output_hashes']['registry.json']).read_text())
        check('environment_registry_completed',env['status']=='completed');check('feature_families_exact',list(env['families'])==cfg['families'])
        for field in env['fields']:
            if 'physical_path' in field:
                path=ROOT/field['physical_path'];track(path,em['input_hashes'][field['physical_path']])
        for path in [cfgpath,Path(__file__),ROOT/'src/spea/environment_benchmark.py',ROOT/'src/spea/genus_model.py',ROOT/'src/spea/genus_spatial.py',ROOT/'src/spea/native_benchmark.py',ROOT/'uv.lock',ROOT/'pyproject.toml']:track(path)
        if prior:track(prior/'run.json')
        with zipfile.ZipFile(out/'source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as z:
            for rel in inputs:
                if rel.endswith(('.py','.json','.toml','.lock')):z.write(ROOT/rel,rel)
        reader=NativeEnvironment(ROOT,env['window_grid']['shape'],env['fields']);shape=reader.shape;rowarea=np.load(prep/'row_area_m2.npy');cells=pd.read_csv(prep/'cells.csv.gz')
        base=json.loads((prep/'feature_registry.json').read_text());check('native_grid_exact',env['window_grid']==base['window_grid'])
        check('positive_native_areas',rowarea.shape==(shape[0],) and np.isfinite(rowarea).all() and (rowarea>0).all())
        allkeys=env['families']['broad_soil_ecoregion'];chunk=cfg['budget']['prediction_chunk_size']
        if not args.full_from_run:
            for si,short in enumerate(cfg['species']):
                for fold in cfg['folds']:
                    prefix=f'{short}_fold{fold}';splitpath=prep/'splits'/(prefix+'.npz')
                    with np.load(splitpath) as split:
                        p=split['train_native'];q=split['q_native'];tp=split['test_native'];qw=split['q_weights_m2'];tq=split['target_native'];tw=split['target_area_m2'];tr=split['train_cell_rows'];trows=split['test_cell_rows']
                    check(prefix+'_row_identity',np.array_equal(p,cells.loc[tr,'native_flat']) and np.array_equal(tp,cells.loc[trows,'native_flat']))
                    check(prefix+'_native_disjoint',not np.isin(q,tq).any() and np.isin(p,q).all() and np.isin(tp,tq).all())
                    check(prefix+'_sorted_target',np.all(np.diff(tq)>0));check(prefix+'_target_areas_exact',np.array_equal(tw,rowarea[tq//shape[1]]))
                    weights=case_control_weights(len(p),qw);seed=cfg['seed']+1000*si+fold
                    check(prefix+'_balanced_classes',np.isclose(weights[:len(p)].sum(),weights[len(p):].sum()) and np.isclose(weights.sum(),len(p)+len(q)))
                    np.savez_compressed(out/'training_design'/(prefix+'.npz'),train_native=p,q_native=q,train_cell_rows=tr,q_weights_m2=qw,class_object_weights=weights)
                    train=reader.read(np.r_[p,q],allkeys)
                    for family in cfg['families']:
                        budget.check();keys=env['families'][family];stem=prefix+'_'+family;t=time.monotonic();metric=next(r for r in metrics if r['species']==short and r['fold']==fold and r['family']==family)
                        audit=dict(kind='cv',species=short,fold=fold,family=family,seed=seed,n_train=len(p),n_quadrature=len(q),row_identity_sha256=digest_array(np.r_[p,q]),class_object_weights_sha256=digest_array(weights));begin_attempt(audit)
                        try:
                            reference=training_reference(train,len(p),keys);model=fit_model(train,len(p),weights,keys,cfg,seed,budget);fit_time=time.monotonic()-t
                            save_model(model,stem,out,keys,reference,weights,seed)
                            predp=predict(model,tp,keys,reference,reader,chunk,budget);predq=predict(model,tq,keys,reference,reader,chunk,budget)
                            check(stem+'_held_point_prediction_matches_target',np.array_equal(predp['raw_eta'],predq['raw_eta'][np.searchsorted(tq,tp)]))
                            scores,norm=relative_scores(predp['raw_eta'],predq['raw_eta'],tw)
                            check(stem+'_unit_area_mean',abs(float(logsumexp(predq['raw_eta']-norm+np.log(tw))-np.log(tw.sum())))<1e-10)
                            np.savez_compressed(out/'predictions'/(stem+'.npz'),test_native=tp,test_cell_rows=trows,target_native=tq,target_area_m2=tw,presence_raw_eta=predp['raw_eta'],target_raw_eta=predq['raw_eta'],presence_log_relative_density=scores,target_log_normalizer=np.asarray(norm),novel_presence=predp['novel_numeric'],novel_target=predq['novel_numeric'],presence_numeric_novel_count=predp['numeric_novel_count'],target_numeric_novel_count=predq['numeric_novel_count'],presence_numeric_missing_count=predp['numeric_missing_count'],target_numeric_missing_count=predq['numeric_missing_count'],unseen_cec_presence=predp['unseen_cec'],unseen_cec_target=predq['unseen_cec'],training_feature_min=reference['minimum'],training_feature_max=reference['maximum'],training_numeric_feature_names=np.asarray(reference['numeric_keys']))
                            tc=cells.loc[trows];blocks=tc.block50.to_numpy();metric.update(status='passed',n_train=len(p),n_test=len(tp),n_test_blocks=len(np.unique(blocks)),n_quadrature=len(q),target_native_cells=len(tq),target_area_km2=float(tw.sum()/1e6),mean_cell_log_score=float(scores.mean()),mean_block_log_score=block_score(scores,blocks),novel_test_fraction=float(predp['novel_numeric'].mean()),novel_target_fraction=float(predq['novel_numeric'].mean()),novel_test_numeric_feature_cells=int(predp['numeric_novel_count'].sum()),novel_target_numeric_feature_cells=int(predq['numeric_novel_count'].sum()),unseen_cec_test_fraction=float(predp['unseen_cec'].mean()),unseen_cec_target_fraction=float(predq['unseen_cec'].mean()),missing_numeric_test_fraction=float((predp['numeric_missing_count']>0).mean()),target_log_normalizer=norm,fit_seconds=fit_time,evaluation_seconds=time.monotonic()-t-fit_time)
                            for i,(index,score) in enumerate(zip(trows,scores)):
                                c=cells.loc[index];cell_scores.append(dict(species=short,fold=fold,family=family,source_grid_cell_id=int(c.source_grid_cell_id),native_flat=int(c.native_flat),longitude=float(c.longitude),latitude=float(c.latitude),block50=int(c.block50),score=float(score),novel_numeric=bool(predp['novel_numeric'][i]),unseen_cec=bool(predp['unseen_cec'][i]),raw_eta=float(predp['raw_eta'][i])))
                            audit.update(status='passed',tree_count=int(model.tree_count_),fit_seconds=fit_time);del model,predp,predq
                        except BudgetExceeded as error:fit_audit_failure(audit,error);metric.update(status='budget_aborted',error=str(error));raise
                        except Exception as error:fit_audit_failure(audit,error);metric.update(status='failed',error=str(error))
                        finally:
                            audit['total_seconds']=time.monotonic()-t;flush();print(json.dumps({k:audit.get(k) for k in ['kind','species','fold','family','status','fit_seconds','total_seconds']}),flush=True)
                    del train
        if not args.cv_only:
            for si,short in enumerate(cfg['species']):
                budget.check();ids=np.load(prep/(short+'_primary_native.npy'));part=cells[cells.working_taxon.eq('Spea '+short)];inside=part.native_flat.isin(ids).to_numpy();p=part.loc[inside,'native_flat'].to_numpy('int64');seed=cfg['seed']+1000*si+100;family='broad_soil';keys=env['families'][family];stem=short+'_broad_soil_full';t=time.monotonic();summary=next(r for r in projections if r['species']==short)
                for (index,c),admitted in zip(part.iterrows(),inside):projection_membership.append(dict(species=short,cell_row=int(index),native_flat=int(c.native_flat),source_grid_cell_id=int(c.source_grid_cell_id),inside_primary_M=bool(admitted)))
                check(short+'_unique_full_presence',len(p)>0 and len(np.unique(p))==len(p));area=rowarea[ids//shape[1]];blocks=np.empty(len(ids),'int64')
                for startrow in range(0,len(ids),chunk):
                    budget.check();end=startrow+chunk;blocks[startrow:end]=geographic_blocks(ids[startrow:end],shape,env['window_grid']['transform'],25000)
                q,qw,pi,qblocks,counts=stratified_quadrature(ids,p,blocks,area,16,seed)
                weights=case_control_weights(len(p),qw)
                np.savez_compressed(out/'training_design'/(stem+'.npz'),train_native=p,q_native=q,q_weights_m2=qw,q_inclusion_probability=pi,q_block25=qblocks,stratum_counts=counts,class_object_weights=weights)
                train=reader.read(np.r_[p,q],allkeys);reference=training_reference(train,len(p),keys)
                audit=dict(kind='full',species=short,fold=-1,family=family,seed=seed,n_train=len(p),n_quadrature=len(q),row_identity_sha256=digest_array(np.r_[p,q]),class_object_weights_sha256=digest_array(weights));begin_attempt(audit)
                try:
                    model=fit_model(train,len(p),weights,keys,cfg,seed,budget);fit_time=time.monotonic()-t;save_model(model,stem,out,keys,reference,weights,seed);pred=predict(model,ids,keys,reference,reader,chunk,budget)
                    _,norm=relative_scores(np.array([0.]),pred['raw_eta'],area);rel=pred['raw_eta']-norm
                    check(stem+'_unit_area_mean',abs(float(logsumexp(rel+np.log(area))-np.log(area.sum())))<1e-10)
                    np.savez_compressed(out/'projection'/(stem+'.npz'),native_flat=ids,area_m2=area,raw_eta=pred['raw_eta'],log_relative_density=rel,novel_numeric=pred['novel_numeric'],numeric_novel_count=pred['numeric_novel_count'],numeric_missing_count=pred['numeric_missing_count'],unseen_cec=pred['unseen_cec'],display_log_normalizer=np.asarray(norm))
                    summary.update(status='passed',training_cells=len(p),canonical_cells=len(part),canonical_cells_outside_M=int((~inside).sum()),quadrature_cells=len(q),quadrature_mass_relative_error=float(qw.sum()/area.sum()-1),projection_native_cells=len(ids),projection_area_km2=float(area.sum()/1e6),novel_fraction=float(pred['novel_numeric'].mean()),novel_numeric_feature_cells=int(pred['numeric_novel_count'].sum()),unseen_cec_fraction=float(pred['unseen_cec'].mean()),missing_numeric_fraction=float((pred['numeric_missing_count']>0).mean()),display_log_normalizer=norm,fit_seconds=fit_time)
                    audit.update(status='passed',tree_count=int(model.tree_count_),fit_seconds=fit_time);del model,pred,train
                except BudgetExceeded as error:fit_audit_failure(audit,error);summary.update(status='budget_aborted',error=str(error));raise
                except Exception as error:fit_audit_failure(audit,error);summary.update(status='failed',error=str(error))
                finally:audit['total_seconds']=time.monotonic()-t;flush();print(json.dumps(audit,default=str),flush=True)
        for rel,h in inputs.items():check('unchanged:'+rel,sha256(ROOT/rel)==h)
        expected_here=4 if args.full_from_run else 100 if args.cv_only else 104
        check('declared_fit_count',record['fit_attempts']==expected_here)
        record.update(status='completed' if all(x['status']=='passed' for x in fits) else 'completed_with_failures',finished_utc=utcnow())
    except BudgetExceeded as error:record.update(status='budget_exhausted',error=str(error),finished_utc=utcnow())
    except BaseException as error:record.update(status='interrupted_or_failed',error=str(error),traceback=traceback.format_exc(),finished_utc=utcnow());raise
    finally:
        flush();write_json(out/'checks.json',dict(passed=sum(c['passed'] for c in checks),failed=sum(not c['passed'] for c in checks),checks=checks));record['output_hashes']=outputs_hash();write_json(out/'run.json',record)
    print(json.dumps({k:record[k] for k in ['status','fit_attempts','cv_fit_attempts','full_data_fit_attempts','cpu_seconds','peak_resident_bytes']}),flush=True)

if __name__=='__main__':main()
