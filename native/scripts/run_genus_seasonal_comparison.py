#!/usr/bin/env python3
"""Frozen nested seasonal representations with bounded spatial BCE selection.

All 60 selections precede outer evaluation. Import and --plan-only fit nothing.
No occurrences, spatial domains, genetic groups or scientific gates are changed.
"""
from __future__ import annotations
import os
for variable in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','MKL_NUM_THREADS']:
    os.environ.setdefault(variable,'1')
import argparse
import hashlib
import json
import resource
import sys
import time
import traceback
import zipfile
from pathlib import Path
import catboost
import numpy as np
import pandas as pd
from spea.provenance import ROOT,sha256,write_json,utcnow
from spea.genus_model import geographic_blocks,stratified_quadrature
from spea.genus_spatial import exact_pair_min_clearance
from spea.genus_tuning import nested_roles
from spea.native_benchmark import relative_scores,block_score
from spea.native_objective_comparison import (
    encode_objective,save_common_borders,fit_native_objective,
    held_diagnostics,select_complete_candidates)
from spea.seasonal_model import MonthlyReader,FAMILIES,prefix_borders


class BudgetExceeded(RuntimeError):
    pass


class Budget:
    def __init__(self,cfg):
        self.start=time.process_time();self.cfg=cfg;self.reason=None
    def cpu(self):
        return time.process_time()-self.start
    def check(self):
        if self.cpu()>=60*self.cfg['cpu_minutes']:
            raise BudgetExceeded('Declared CPU budget exhausted')
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>=self.cfg['peak_gib']*1024**3:
            raise BudgetExceeded('Declared resident-memory budget exhausted')
    def after_iteration(self,info):
        try:
            self.check();return True
        except BudgetExceeded as error:
            self.reason=str(error);return False


def array_digest(value):
    a=np.ascontiguousarray(value)
    return hashlib.sha256(a.dtype.str.encode()+str(a.shape).encode()+a.tobytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',default='config/genus_seasonal_comparison.json')
    ap.add_argument('--run-id',required=True);ap.add_argument('--plan-only',action='store_true')
    args=ap.parse_args()
    if Path(args.run_id).name!=args.run_id:
        raise ValueError('Run ID must be one path component')
    cfgpath=ROOT/args.config;cfg=json.loads(cfgpath.read_text())
    if args.plan_only:
        print(json.dumps(cfg,indent=2));return
    if catboost.__version__!='1.2.10':
        raise ValueError('Require audited CatBoost 1.2.10')
    if (cfg['families']!=list(FAMILIES) or cfg['family_dimensions']!=[2,4,28]
            or len(cfg['species'])!=4 or len(cfg['candidates'])!=4):
        raise ValueError('Unexpected bounded seasonal design')
    out=ROOT/'runs'/args.run_id;out.mkdir(exist_ok=False)
    for folder in ['splits','borders','models','predictions']:
        (out/folder).mkdir()
    start=time.monotonic();budget=Budget(cfg['budget'])
    inputs={};checks=[];fit_rows=[];inner_rows=[];selection_rows=[];outer_rows=[];capture_rows=[];cell_rows=[]
    record=dict(status='preparing',created_utc=utcnow(),command_argv=sys.argv,config=cfg,
        input_hashes=inputs,fit_attempts=0,inner_fit_attempts=0,outer_fit_attempts=0,
        empirical_fits=0,full_data_fits=0,formal_significance_tests=0,model_input_changes=0,
        scientific_gate_changes=0,new_genetic_data_or_analysis=False)
    write_json(out/'effective_config.json',cfg)

    def check(name,value):
        checks.append(dict(name=name,passed=bool(value)))
        if not value:
            raise ValueError(name)
    def track(path,expected=None):
        path=Path(path);rel=str(path.relative_to(ROOT));digest=sha256(path)
        check('hash:'+rel,expected is None or digest==expected)
        if rel in inputs:
            check('repeat_input:'+rel,inputs[rel]==digest)
        inputs[rel]=digest;return path
    manifests={}
    def asset(directory,name):
        directory=Path(directory)
        if directory not in manifests:
            manifests[directory]=json.loads(track(directory/'run.json').read_text())
            expected=('completed_diagnostic_no_fits'
                      if directory==ROOT/'runs'/cfg['spatial_run'] else 'completed')
            check('source_completed:'+str(directory.relative_to(ROOT)),manifests[directory]['status']==expected)
        return track(directory/name,manifests[directory]['output_hashes'][name])
    def flush():
        for name,rows in [('fit_audit.csv',fit_rows),('inner_metrics.csv',inner_rows),
                          ('selection.csv',selection_rows),('metrics.csv',outer_rows),
                          ('area_capture.csv',capture_rows),('cell_scores.csv.gz',cell_rows)]:
            if rows:
                pd.DataFrame(rows).to_csv(out/name,index=False)
        record.update(cpu_seconds=budget.cpu(),elapsed_seconds=time.monotonic()-start,
                      peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write_json(out/'run.json',record)
    def source_snapshot():
        with zipfile.ZipFile(out/'source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as z:
            for rel in inputs:
                if rel.endswith(('.py','.json','.toml','.lock')):
                    z.write(ROOT/rel,rel)

    try:
        prep=ROOT/'runs'/cfg['preparation_run'];environment=ROOT/'runs'/cfg['environment_run']
        spatial=ROOT/'runs'/cfg['spatial_run'];old=ROOT/'runs'/cfg['inner_reference_run']
        cells=pd.read_csv(asset(prep,'cells.csv.gz'))
        base=json.loads(asset(prep,'feature_registry.json').read_text())
        prep_cfg=json.loads(asset(prep,'effective_config.json').read_text())
        env=json.loads(asset(environment,'registry.json').read_text())
        check('same_native_grid',env['window_grid']==base['window_grid'])
        reader=MonthlyReader(ROOT,env);shape=reader.shape
        for paths in reader.paths.values():
            for path in paths:
                rel=str(path.relative_to(ROOT));track(path,manifests[environment]['input_hashes'][rel])
        old_record=json.loads(track(old/'run.json').read_text())
        check('old_nested_completed',old_record['status']=='completed')
        for rel in [str((prep/'effective_config.json').relative_to(ROOT)),
                    str((prep/'cells.csv.gz').relative_to(ROOT)),
                    str((spatial/'cell_fold_diagnostics.csv.gz').relative_to(ROOT))]:
            if rel not in inputs:
                track(ROOT/rel)
            check('current_inputs_match_old_nested:'+rel,inputs[rel]==old_record['input_hashes'][rel])
        track(cfgpath)
        for rel in ['scripts/run_genus_seasonal_comparison.py','src/spea/seasonal_model.py',
                    'src/spea/native_objective_comparison.py','src/spea/objective_measures.py',
                    'src/spea/genus_model.py','src/spea/genus_spatial.py','src/spea/genus_tuning.py',
                    'src/spea/native_benchmark.py','tests/test_seasonal_model.py',
                    'tests/test_native_objective_comparison.py','pyproject.toml','uv.lock']:
            track(ROOT/rel)
        rowarea=np.load(asset(prep,'row_area_m2.npy'))
        role_table=pd.read_csv(asset(spatial,'cell_fold_diagnostics.csv.gz'))
        safepath=asset(spatial,'safe_background_global.npz')
        heldpath=asset(spatial,'test_background_global.npz')
        designs=[];supported=[];actual={}
        for species in cfg['species']:
            primary=np.load(asset(prep,species+'_primary_native.npy'))
            part=cells[cells.working_taxon.eq('Spea '+species)]
            roles=role_table[role_table.working_taxon.eq('Spea '+species)&role_table.variant.eq('primary')]
            lookup={f:roles[roles.fold.eq(f)].set_index('source_grid_cell_id').loc[part.source_grid_cell_id]
                    for f in cfg['folds']}
            with np.load(safepath) as source:
                safe={f:source[str(f)].ravel()[primary] for f in cfg['folds']}
            with np.load(heldpath) as source:
                held={f:source[str(f)].ravel()[primary] for f in cfg['folds']}
            blocks25=geographic_blocks(primary,shape,env['window_grid']['transform'],
                                       prep_cfg['quadrature']['stratum_width_m'])
            for outer in cfg['folds']:
                outerpath=asset(prep,f'splits/{species}_fold{outer}.npz')
                with np.load(outerpath) as s:
                    payload={k:s[k] for k in ['train_native','test_native','train_cell_rows',
                              'test_cell_rows','q_native','q_weights_m2','q_inclusion_probability',
                              'target_native','target_area_m2']}
                prefix=f'{species}_fold{outer}'
                check(prefix+':occupied_role_identity',np.array_equal(payload['train_native'],
                    part.loc[lookup[outer].train_1km.to_numpy(),'native_flat'].to_numpy()))
                check(prefix+':held_domain_identity',np.array_equal(payload['target_native'],primary[held[outer]]))
                check(prefix+':event_inclusive_unique',np.all(np.diff(payload['q_native'])>0)
                      and np.isin(payload['train_native'],payload['q_native']).all())
                gap=exact_pair_min_clearance(cells.loc[payload['train_cell_rows'],'longitude'].to_numpy(),
                    cells.loc[payload['train_cell_rows'],'latitude'].to_numpy(),
                    cells.loc[payload['train_cell_rows'],'radius_unknown1km_m'].to_numpy(),
                    cells.loc[payload['test_cell_rows'],'longitude'].to_numpy(),
                    cells.loc[payload['test_cell_rows'],'latitude'].to_numpy(),
                    cells.loc[payload['test_cell_rows'],'radius_unknown5km_m'].to_numpy())
                check(prefix+':support_gap50km',np.all(gap>=50000-1e-6))
                np.savez_compressed(out/'splits'/(prefix+'.npz'),**payload)
                for inner in cfg['folds']:
                    if inner==outer:
                        continue
                    tr,va=nested_roles(lookup[outer].train_1km,lookup[inner].train_1km,
                                       lookup[inner].test_common_support)
                    trainrows=part.index[tr].to_numpy();valrows=part.index[va].to_numpy()
                    p=cells.loc[trainrows,'native_flat'].to_numpy('int64')
                    vp=cells.loc[valrows,'native_flat'].to_numpy('int64')
                    allowed=safe[outer]&safe[inner];target_mask=safe[outer]&held[inner]
                    eligible=primary[allowed];target=primary[target_mask]
                    prefix=f'{species}_outer{outer}_inner{inner}'
                    check(prefix+':events_in_domains',np.isin(p,eligible).all() and np.isin(vp,target).all())
                    check(prefix+':inner_domains_disjoint',not np.isin(eligible,target).any())
                    check(prefix+':domains_avoid_outer_holdout',not np.any((allowed|target_mask)&held[outer]))
                    gap=exact_pair_min_clearance(cells.loc[trainrows,'longitude'].to_numpy(),
                        cells.loc[trainrows,'latitude'].to_numpy(),cells.loc[trainrows,'radius_unknown1km_m'].to_numpy(),
                        cells.loc[valrows,'longitude'].to_numpy(),cells.loc[valrows,'latitude'].to_numpy(),
                        cells.loc[valrows,'radius_unknown5km_m'].to_numpy())
                    check(prefix+':support_gap50km',np.all(gap>=50000-1e-6))
                    nblocks=cells.loc[valrows,'block50'].nunique()
                    use=len(p)>=cfg['minimum_inner_training_cells'] and nblocks>=cfg['minimum_inner_validation_blocks']
                    design=dict(species=species,prefix=prefix,outer_fold=outer,inner_fold=inner,train_cells=len(p),
                        validation_cells=len(vp),validation_blocks=int(nblocks),eligible=bool(use),
                        min_support_gap_m=float(gap.min()) if gap.size else None)
                    designs.append(design)
                    if not use:
                        continue
                    q,w,pi,_,_=stratified_quadrature(eligible,p,blocks25[allowed],
                        rowarea[eligible//shape[1]],prep_cfg['quadrature']['nonpresence_per_stratum'],
                        prep_cfg['quadrature']['seed']+1000*prep_cfg['species'].index(species)+10*outer+inner)
                    payload=dict(train_native=p,validation_native=vp,train_cell_rows=trainrows,
                        validation_cell_rows=valrows,q_native=q,q_weights_m2=w,q_inclusion_probability=pi,
                        target_native=target,target_area_m2=rowarea[target//shape[1]])
                    with np.load(asset(old,'inner_splits/'+prefix+'.npz')) as cached:
                        for key,value in payload.items():
                            check(prefix+':exact_cached_'+key,np.array_equal(value,cached[key]))
                    np.savez_compressed(out/'splits'/(prefix+'.npz'),**payload);supported.append(design)
            actual[species]=[sum(d['outer_fold']==f and d['species']==species for d in supported) for f in cfg['folds']]
            check(species+':same_eligible_inner_folds',actual[species]==cfg['expected_inner_folds'][species])
            del primary,safe,held,blocks25
            budget.check()
        pd.DataFrame(designs).to_csv(out/'inner_design.csv',index=False)
        check('78_inner_splits',len(supported)==78)
        source_snapshot()
        record.update(status='inner_search',eligible_inner_folds=actual,design_frozen_utc=utcnow())
        flush();print(json.dumps(dict(stage='design_frozen',eligible=actual)),flush=True)

        def fit_and_predict(split,frame,target_frame,candidate,family,stem,borders,seed,kind,species,prefix):
            budget.check();budget.reason=None
            if (record['fit_attempts']>=cfg['budget']['maximum_fits'] or
                    record[kind+'_fit_attempts']>=cfg['budget']['maximum_'+kind+'_fits']):
                raise BudgetExceeded('Declared fit-count budget exhausted')
            encoding=encode_objective('bce',split['train_native'],split['q_native'],
                                      split['q_weights_m2'],candidate['leaf_damping_per_presence'])
            audit=dict(stem=stem,split=prefix,species=species,family=family,kind=kind,objective='bce',
                candidate=candidate['id'],seed=seed,status='running',started_utc=utcnow(),
                native_identity_sha256=array_digest(encoding.native),weights_sha256=array_digest(encoding.weight),
                target_sha256=array_digest(encoding.target),borders_sha256=sha256(borders))
            record['fit_attempts']+=1;record[kind+'_fit_attempts']+=1
            record['empirical_fits']+=1;fit_rows.append(audit)
            cpu=time.process_time();wall=time.monotonic()
            try:
                model,details=fit_native_objective(frame,encoding,'bce',candidate,borders,seed,
                    callbacks=[budget],threads=cfg['budget']['threads'])
                if budget.reason:
                    raise BudgetExceeded(budget.reason)
                budget.check();fit_seconds=time.monotonic()-wall
                model.save_model(str(out/'models'/(stem+'.cbm')))
                eta=np.empty(len(target_frame))
                for startrow in range(0,len(eta),cfg['budget']['prediction_chunk_size']):
                    budget.check();end=startrow+cfg['budget']['prediction_chunk_size']
                    eta[startrow:end]=model.predict(target_frame.iloc[startrow:end],
                        prediction_type='RawFormulaVal',thread_count=cfg['budget']['threads'])
                if not np.isfinite(eta).all():
                    raise ValueError('Nonfinite held raw margins')
                details.update(candidate=candidate,seed=seed,family=family,split=prefix,species=species,
                    feature_names=FAMILIES[family],borders_path=str(borders.relative_to(out)),
                    borders_sha256=sha256(borders),native_identity_sha256=audit['native_identity_sha256'],
                    weights_sha256=audit['weights_sha256'],target_sha256=audit['target_sha256'])
                write_json(out/'models'/(stem+'_metadata.json'),details)
                audit.update(status='passed',finished_utc=utcnow(),fit_seconds=fit_seconds,
                    total_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu,
                    native_l2_leaf_reg=details['native_l2_leaf_reg'],effective_leaf_damping=details['effective_leaf_damping'],
                    training_metric_error=details['training_metric_error'])
                return eta
            except Exception as error:
                audit.update(status='failed',error=str(error),traceback=traceback.format_exc(),finished_utc=utcnow(),
                    total_seconds=time.monotonic()-wall,cpu_seconds=time.process_time()-cpu)
                if budget.reason:
                    raise BudgetExceeded(budget.reason) from error
                if isinstance(error,BudgetExceeded):
                    raise
                return None

        def borders_for(prefix,shared,seed):
            full=out/'borders'/(prefix+'_monthly.tsv')
            save_common_borders(shared,full,seed=seed)
            result={'monthly':full}
            for family in ['annual','variability']:
                path=out/'borders'/(prefix+'_'+family+'.tsv')
                prefix_borders(full,path,len(FAMILIES[family]));result[family]=path
            return result

        for design in supported:
            species=design['species'];outer=design['outer_fold'];inner=design['inner_fold'];prefix=design['prefix']
            with np.load(out/'splits'/(prefix+'.npz')) as source:
                split={k:source[k] for k in source.files}
            shared=reader.read(np.r_[split['train_native'],split['q_native']])
            target_frame=reader.read(split['target_native'])
            positions=np.searchsorted(split['target_native'],split['validation_native'])
            check(prefix+':held_point_positions',np.array_equal(split['target_native'][positions],split['validation_native']))
            blocks=cells.loc[split['validation_cell_rows'],'block50'].to_numpy()
            seed=cfg['model_seed']+1000*cfg['species'].index(species)+10*outer+inner
            borders=borders_for(prefix,shared,seed)
            auditidx=np.unique(np.linspace(0,len(target_frame)-1,min(256,len(target_frame)),dtype=int))
            for family in cfg['families']:
                keys=FAMILIES[family];frame=shared[keys];targetfamily=target_frame[keys]
                for candidate in cfg['candidates']:
                    stem=prefix+'_'+family+'_'+candidate['id']
                    eta=fit_and_predict(split,frame,targetfamily,candidate,family,stem,borders[family],seed,'inner',species,prefix)
                    row=dict(species=species,family=family,outer_fold=outer,inner_fold=inner,candidate=candidate['id'],
                        validation_cells=len(positions),validation_blocks=design['validation_blocks'],seed=seed,
                        stem=stem,split=prefix,status='failed' if eta is None else 'passed')
                    if eta is not None:
                        scores,norm=relative_scores(eta[positions],eta,split['target_area_m2'])
                        row.update(mean_block_log_score=block_score(scores,blocks),mean_cell_log_score=float(scores.mean()),normalizer=norm)
                        np.savez_compressed(out/'predictions'/(stem+'.npz'),presence_raw_eta=eta[positions],
                            presence_log_relative_density=scores,target_log_normalizer=np.asarray(norm),
                            audit_target_positions=auditidx,audit_target_native=split['target_native'][auditidx],
                            audit_target_raw_eta=eta[auditidx])
                    inner_rows.append(row)
            flush();print(json.dumps(dict(stage='inner_done',species=species,outer=outer,inner=inner,
                fit_attempts=record['fit_attempts'],cpu_seconds=round(budget.cpu(),1))),flush=True)
            del shared,target_frame,frame,targetfamily
        candidate_ids=[c['id'] for c in cfg['candidates']]
        inner_frame=pd.DataFrame(inner_rows);selected={}
        for species in cfg['species']:
            for outer in cfg['folds']:
                folds=[d['inner_fold'] for d in supported if d['outer_fold']==outer and d['species']==species]
                for family in cfg['families']:
                    matrix=np.full((len(candidate_ids),len(folds)),np.nan)
                    for i,candidate in enumerate(candidate_ids):
                        for j,inner in enumerate(folds):
                            row=inner_frame[(inner_frame.species==species)&(inner_frame.outer_fold==outer)&
                                (inner_frame.inner_fold==inner)&(inner_frame.family==family)&(inner_frame.candidate==candidate)].iloc[0]
                            if row.status=='passed':
                                matrix[i,j]=row.mean_block_log_score
                    chosen,reason=select_complete_candidates(matrix,candidate_ids,cfg['minimum_valid_inner_folds'])
                    selected[(species,outer,family)]=next(c for c in cfg['candidates'] if c['id']==chosen)
                    for i,candidate in enumerate(candidate_ids):
                        complete=bool(np.isfinite(matrix[i]).all())
                        selection_rows.append(dict(species=species,outer_fold=outer,family=family,candidate=candidate,
                            selected=candidate==chosen,eligible_inner_folds=len(folds),complete=complete,
                            inner_mean=float(matrix[i].mean()) if complete else None,selection_reason=reason))
        flush();record.update(status='outer_evaluation',selection_frozen_utc=utcnow(),selection_sha256=sha256(out/'selection.csv'))
        flush();print(json.dumps(dict(stage='all_60_selections_frozen',sha256=record['selection_sha256'])),flush=True)
        for species in cfg['species']:
            for fold in cfg['folds']:
                prefix=f'{species}_fold{fold}'
                with np.load(out/'splits'/(prefix+'.npz')) as source:
                    split={k:source[k] for k in source.files}
                p=split['train_native'];target=split['target_native'];tp=split['test_native']
                shared=reader.read(np.r_[p,split['q_native']]);target_frame=reader.read(target)
                positions=np.searchsorted(target,tp);check(prefix+':held_point_positions',np.array_equal(target[positions],tp))
                blocks=cells.loc[split['test_cell_rows'],'block50'].to_numpy()
                seed=cfg['model_seed']+1000*cfg['species'].index(species)+fold
                borders=borders_for(prefix,shared,seed)
                for family in cfg['families']:
                    keys=FAMILIES[family];frame=shared[keys];targetfamily=target_frame[keys]
                    lower=frame.iloc[len(p):].min().to_numpy();upper=frame.iloc[len(p):].max().to_numpy()
                    novelq=((targetfamily.to_numpy()<lower)|(targetfamily.to_numpy()>upper)).any(axis=1)
                    novelp=novelq[positions];fixed_eta=None;fixed_stem=None
                    for variant,candidate in [('fixed',cfg['candidates'][0]),('selected',selected[(species,fold,family)])]:
                        stem=prefix+'_'+family+'_'+variant;alias=None
                        if variant=='selected' and candidate['id']=='baseline' and fixed_eta is not None:
                            eta=fixed_eta;alias=fixed_stem
                        else:
                            eta=fit_and_predict(split,frame,targetfamily,candidate,family,stem,borders[family],seed,'outer',species,prefix)
                        row=dict(species=species,fold=fold,family=family,variant=variant,candidate=candidate['id'],
                            stem=stem,split=prefix,n_train=len(p),n_quadrature=len(split['q_native']),n_test=len(tp),
                            n_test_blocks=len(np.unique(blocks)),target_native_cells=len(target),
                            target_area_km2=float(split['target_area_m2'].sum()/1e6),model_alias=alias,
                            status='failed' if eta is None else 'passed')
                        if eta is not None:
                            diagnostics,scores,captures=held_diagnostics(eta[positions],eta,split['target_area_m2'],blocks,novelp)
                            row.update(diagnostics,novel_target_fraction=float(novelq.mean()))
                            np.savez_compressed(out/'predictions'/(stem+'.npz'),target_raw_eta=eta,
                                presence_raw_eta=eta[positions],presence_log_relative_density=scores,
                                target_log_normalizer=np.asarray(diagnostics['target_log_normalizer']),
                                novel_presence=novelp,novel_target=novelq,training_feature_minimum=lower,
                                training_feature_maximum=upper)
                            for capture in captures:
                                capture_rows.append(dict(species=species,fold=fold,family=family,variant=variant,**capture))
                            for cell,score,novel in zip(split['test_cell_rows'],scores,novelp):
                                cell_rows.append(dict(species=species,fold=fold,family=family,variant=variant,cell_row=int(cell),
                                    native_flat=int(cells.loc[cell,'native_flat']),block50=int(cells.loc[cell,'block50']),
                                    score=float(score),novel=bool(novel)))
                            if variant=='fixed':
                                fixed_eta=eta;fixed_stem=stem
                        outer_rows.append(row)
                flush();budget.check()
                if sum(p.stat().st_size for p in out.rglob('*') if p.is_file())>cfg['budget']['output_gb']*1024**3:
                    raise BudgetExceeded('Declared output-size budget exhausted')
                print(json.dumps(dict(stage='outer_done',species=species,fold=fold,fit_attempts=record['fit_attempts'],
                    cpu_seconds=round(budget.cpu(),1))),flush=True)
                del shared,target_frame,frame,targetfamily
        check('all_declared_inner_slots',len(inner_rows)==936)
        check('all_declared_outer_slots',len(outer_rows)==120)
        check('all_outer_evaluations_passed',all(r['status']=='passed' for r in outer_rows))
        check('selection_unchanged',sha256(out/'selection.csv')==record['selection_sha256'])
        for rel,digest in inputs.items():
            check('unchanged:'+rel,sha256(ROOT/rel)==digest)
        check('within_fit_budget',record['fit_attempts']<=cfg['budget']['maximum_fits']);budget.check()
        output_bytes=sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
        check('within_output_budget',output_bytes<=cfg['budget']['output_gb']*1024**3)
        record.update(status='completed' if all(r['status']=='passed' for r in fit_rows) else 'completed_with_inner_failures',
            finished_utc=utcnow(),output_bytes=output_bytes,outer_aliases=sum(r.get('model_alias') is not None for r in outer_rows))
    except BudgetExceeded as error:
        record.update(status='budget_exhausted',error=str(error),finished_utc=utcnow())
    except BaseException as error:
        record.update(status='interrupted_or_failed',error=str(error),traceback=traceback.format_exc(),finished_utc=utcnow())
        raise
    finally:
        flush();write_json(out/'checks.json',dict(passed=sum(r['passed'] for r in checks),
            failed=sum(not r['passed'] for r in checks),checks=checks))
        record['output_hashes']={str(p.relative_to(out)):sha256(p) for p in out.rglob('*') if p.is_file() and p.name!='run.json'}
        write_json(out/'run.json',record)
    print(json.dumps({k:record[k] for k in ['status','fit_attempts','inner_fit_attempts',
        'outer_fit_attempts','cpu_seconds','peak_resident_bytes']}),flush=True)


if __name__=='__main__':
    main()
