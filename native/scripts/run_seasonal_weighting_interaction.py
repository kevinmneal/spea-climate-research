#!/usr/bin/env python3
"""Matched fixed-shallow seasonal representation × saved spatial weights."""
import os
for key in ['OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'MKL_NUM_THREADS']:
    os.environ.setdefault(key, '1')

import argparse
from dataclasses import replace
import json
from pathlib import Path
import resource
import sys
import time
import traceback
import zipfile

import catboost
import numpy as np
import pandas as pd

from run_climate_statistics_comparison import Budget
from spea.seasonal_model import FAMILIES, MONTHLY, MonthlyReader, prefix_borders
from spea.spatial_thinning import inherited_cell_rows
from spea.native_objective_comparison import encode_objective, fit_native_objective, held_diagnostics
from spea.objective_measures import binary_class_weights
from spea.provenance import ROOT, sha256, utcnow, write_json


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', default='config/seasonal_weighting_interaction_v1.json')
    ap.add_argument('--run-id', required=True)
    args = ap.parse_args()
    if Path(args.run_id).name != args.run_id:
        raise ValueError('Run ID must be one component')
    config_path = ROOT / args.config
    cfg = json.loads(config_path.read_text())
    parent = ROOT / 'runs' / cfg['parent_run']
    soft = ROOT / 'runs' / cfg['soft_run']
    weights_dir = ROOT / 'runs' / cfg['weight_run']
    prec, srec, wrec = [json.loads((p/'run.json').read_text()) for p in [parent, soft, weights_dir]]
    if (catboost.__version__ != '1.2.10' or prec['status'] != 'completed' or srec['status'] != 'completed'
            or wrec['status'] != 'completed_diagnostic_no_fits' or wrec['unsupported_slots'] != 0
            or cfg['species'] != prec['config']['species'] or cfg['species'] != srec['config']['species']
            or cfg['families'] != list(FAMILIES) or cfg['new_fit_families'] != ['annual', 'variability']
            or cfg['sampling_arms'] != ['equal_cells', 'soft_geo_10km']
            or cfg['model'] != prec['config']['model'] or cfg['model'] != srec['config']['model']
            or cfg['model_seed'] != prec['config']['model_seed'] or cfg['model_seed'] != srec['config']['model_seed']):
        raise ValueError('Frozen matched design or parents changed')
    out = ROOT / 'runs' / args.run_id
    out.mkdir(exist_ok=False)
    for folder in ['models', 'borders', 'predictions']:
        (out/folder).mkdir()
    budget = Budget(cfg['limits'])
    start = time.monotonic()
    inputs, stats, metrics, fits, reuses, design = {}, {}, [], [], [], []
    record = dict(status='preparing', created_utc=utcnow(), config=cfg, command_argv=sys.argv,
        input_hashes=inputs, empirical_fits=0, reused_monthly_slots=0, declared_new_fits=84,
        full_data_fits=0, historical_projections=0, scientific_gate_changes=0,
        occurrence_changes=0, background_changes=0, weight_vector_changes=0)

    def track(path, expected=None):
        rel = str(path.relative_to(ROOT))
        if rel not in inputs:
            inputs[rel] = sha256(path)
            stat = path.stat(); stats[rel] = (stat.st_size, stat.st_mtime_ns)
        if expected is not None and inputs[rel] != expected:
            raise ValueError('Changed frozen input ' + rel)
        return path

    def inherited(run, name):
        path = ROOT/'runs'/run/name
        rel = str(path.relative_to(ROOT))
        expected = prec['input_hashes'][rel]
        if expected != srec['input_hashes'][rel]:
            raise ValueError('Monthly parents disagree on input ' + rel)
        return track(path, expected)

    def source_asset(directory, manifest, name):
        return track(directory/name, manifest['output_hashes'][name])

    def flush():
        for name, rows in [('metrics.csv', metrics), ('fit_audit.csv', fits), ('monthly_reuse.csv', reuses)]:
            if rows:
                pd.DataFrame(rows).to_csv(out/name, index=False)
        record.update(cpu_seconds=time.process_time()-budget.start, elapsed_seconds=time.monotonic()-start,
                      peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write_json(out/'run.json', record)

    try:
        track(config_path)
        for p in [parent, soft, weights_dir]: track(p/'run.json')
        registry = json.loads(inherited(cfg['monthly_run'], 'registry.json').read_text())
        reader = MonthlyReader(ROOT, registry)
        for paths in reader.paths.values():
            for path in paths:
                rel = str(path.relative_to(ROOT))
                if prec['input_hashes'][rel] != srec['input_hashes'][rel]:
                    raise ValueError('Different physical monthly inputs')
                track(path, prec['input_hashes'][rel])
        genus = pd.read_csv(inherited(cfg['preparation_run'], 'cells.csv.gz'))
        ns = pd.read_csv(inherited(cfg['north_south_run'], 'species_cells.csv'))
        code = [Path(__file__), ROOT/'scripts/run_climate_statistics_comparison.py', ROOT/'src/spea/seasonal_model.py',
                ROOT/'src/spea/native_objective_comparison.py', ROOT/'src/spea/objective_measures.py',
                ROOT/'src/spea/native_benchmark.py', ROOT/'src/spea/spatial_thinning.py',
                ROOT/'src/spea/provenance.py', ROOT/'catalog/review_gates.json', ROOT/'uv.lock', ROOT/'pyproject.toml']
        for path in code: track(path, prec['input_hashes'].get(str(path.relative_to(ROOT))))
        original = pd.read_csv(source_asset(parent, prec, 'metrics.csv'))
        weighted = pd.read_csv(source_asset(soft, srec, 'metrics.csv'))
        old = {'equal_cells': original[(original.family=='monthly') & (original.sampling=='equal_cells')],
               'soft_geo_10km': weighted[weighted.family=='monthly']}
        if any(len(v)!=21 or v.duplicated(['species','fold']).any() for v in old.values()):
            raise ValueError('Missing or duplicate monthly reuse slots')
        # Complete identity checks for both old monthly arms before new fitting.
        for index, species in enumerate(cfg['species']):
            is_ns = species in ['hammondii_north','spea_south']
            for fold in range(3 if is_ns else 5):
                prefix = f'{species}_fold{fold}'
                source_run = cfg['north_south_run'] if is_ns else cfg['preparation_run']
                split_path = inherited(source_run, f'splits/{prefix}.npz')
                with np.load(split_path) as z:
                    p, q, held = z['train_native'], z['q_native'], z['test_native']
                    target = z['target_native']; held_rows = inherited_cell_rows(ns if is_ns else genus, z['test_cell_rows'])
                    if not np.array_equal(held_rows.native_flat, held) or not np.isin(p,q).all() or not np.isin(held,target).all():
                        raise ValueError('Inherited row identity changed')
                    row_counts = (len(p), len(held), len(np.unique(held_rows.block50)), len(target))
                wp = source_asset(weights_dir, wrec, f'weights/{prefix}.npz')
                with np.load(wp) as z:
                    w = z['presence_weights']
                    if not np.array_equal(z['eligible_native'],p) or np.any(w<=0) or not np.isclose(w.sum(),len(p)):
                        raise ValueError('Frozen weighting membership/mass changed')
                mask_dir = ROOT/'runs'/prec['config']['sampling_preparation_run']
                mask_rec = json.loads(track(mask_dir/'run.json').read_text())
                mask_path = source_asset(mask_dir,mask_rec,f'masks/{prefix}.npz')
                with np.load(mask_path) as z:
                    if not np.array_equal(z['eligible_native'],p) or not z['equal_cells_mask'].all():
                        raise ValueError('Equal-cell source removed original presences')
                border = source_asset(parent,prec,f'borders/{prefix}_monthly.tsv')
                source_asset(soft,srec,f'borders/{prefix}_monthly.tsv')
                if sha256(border)!=sha256(soft/f'borders/{prefix}_monthly.tsv'):
                    raise ValueError('Monthly cutoff identities differ')
                for family in cfg['new_fit_families']:
                    prefix_borders(border,out/f'borders/{prefix}_{family}.tsv',len(FAMILIES[family]))
                for arm, table in old.items():
                    row = table[(table.species==species)&(table.fold==fold)].iloc[0].to_dict()
                    directory, manifest = (parent,prec) if arm=='equal_cells' else (soft,srec)
                    mpath = source_asset(directory,manifest,'models/'+row['stem']+'_metadata.json')
                    meta = json.loads(mpath.read_text())
                    if (meta['feature_names']!=MONTHLY or meta['candidate']!=cfg['model']
                            or meta['seed']!=cfg['model_seed']+100*index+fold
                            or meta['source_split']!=str(split_path.relative_to(ROOT))
                            or meta['borders_sha256']!=sha256(border)
                            or meta['presence_mass']!=len(p) or meta['object_rows']!=len(p)+len(q)
                            or not np.isclose(meta['effective_leaf_damping'],.03*len(p))
                            or not np.isclose(meta['total_object_weight'],2*len(p))
                            or meta['tree_count']!=150 or not meta['no_held_data_in_fit']
                            or tuple(row[k] for k in ['train_cells','held_cells','held_blocks','target_cells'])!=row_counts):
                        raise ValueError('Monthly reuse is not matched: '+row['stem'])
                    if arm=='soft_geo_10km' and meta['frozen_weight_vector_sha256']!=sha256(wp):
                        raise ValueError('Different soft weight vector in monthly fit')
                    source_asset(directory,manifest,'models/'+row['stem']+'.cbm')
                    pp = source_asset(directory,manifest,'predictions/'+row['stem']+'.npz')
                    row.update(model_source_run=directory.name,reused=True)
                    metrics.append(row)
                    reuses.append(dict(species=species,fold=fold,sampling=arm,family='monthly',stem=row['stem'],
                        model_source_run=directory.name,source_split=str(split_path.relative_to(ROOT)),
                        metadata_sha256=sha256(mpath),predictions_path=str(pp.relative_to(ROOT)),
                        predictions_sha256=sha256(pp),borders_sha256=sha256(border),
                        weight_vector_sha256=sha256(wp) if arm=='soft_geo_10km' else None))
                design.append(dict(species=species,fold=fold,train_cells=row_counts[0],held_cells=row_counts[1],
                    held_blocks=row_counts[2],target_cells=row_counts[3],source_split=str(split_path.relative_to(ROOT)),
                    source_split_sha256=sha256(split_path),monthly_borders_sha256=sha256(border),weight_sha256=sha256(wp)))
        write_json(out/'frozen_design.json',dict(status='frozen_before_new_fits',slots=design,reused_monthly_slots=42,
                   new_fits=84,total_slots=126,config_sha256=sha256(config_path)))
        write_json(out/'effective_config.json',cfg)
        with zipfile.ZipFile(out/'source_snapshot.zip','w',zipfile.ZIP_DEFLATED) as z:
            for path in code+[config_path]: z.write(path,str(path.relative_to(ROOT)))
        record.update(status='fitting',design_frozen_utc=utcnow(),reused_monthly_slots=42)
        flush()
        for index,species in enumerate(cfg['species']):
            is_ns = species in ['hammondii_north','spea_south']
            for fold in range(3 if is_ns else 5):
                budget.check(); prefix=f'{species}_fold{fold}'
                split_path=ROOT/'runs'/(cfg['north_south_run'] if is_ns else cfg['preparation_run'])/f'splits/{prefix}.npz'
                with np.load(split_path) as z:
                    split={k:z[k] for k in ['train_native','q_native','q_weights_m2','test_native','test_cell_rows',
                                           'target_native','target_area_m2']}
                n=len(split['train_native']);encoding=encode_objective('bce',split['train_native'],split['q_native'],
                    split['q_weights_m2'],cfg['model']['leaf_damping_per_presence'])
                training=reader.read(encoding.native);target=reader.read(split['target_native'])
                positions=np.searchsorted(split['target_native'],split['test_native'])
                held_rows=inherited_cell_rows(ns if is_ns else genus,split['test_cell_rows']);blocks=held_rows.block50.to_numpy()
                wp=weights_dir/f'weights/{prefix}.npz'
                with np.load(wp) as z: positive=z['presence_weights'].copy()
                seed=cfg['model_seed']+100*index+fold
                audit_idx=np.unique(np.linspace(0,len(target)-1,min(256,len(target)),dtype=int))
                for arm in cfg['sampling_arms']:
                    current=encoding
                    if arm=='soft_geo_10km':
                        w=binary_class_weights(positive,split['q_weights_m2'],loss_mass=n)
                        empirical=np.zeros(len(split['q_native']));empirical[np.searchsorted(split['q_native'],split['train_native'])]=w[:n]
                        current=replace(encoding,weight=w,empirical_on_q=empirical,
                                        native_l2_leaf_reg=encoding.effective_leaf_damping/w.mean())
                    for family in cfg['new_fit_families']:
                        budget.check()
                        if record['empirical_fits']>=84: raise RuntimeError('Declared fit limit reached')
                        names=FAMILIES[family];borders=out/f'borders/{prefix}_{family}.tsv'
                        stem=f'{prefix}_{arm}_{family}';tick=time.monotonic();budget.reason=None
                        record['empirical_fits']+=1
                        model,details=fit_native_objective(training[names],current,'bce',cfg['model'],borders,seed,
                                                          callbacks=[budget],threads=2)
                        if budget.reason: raise RuntimeError(budget.reason)
                        model.save_model(str(out/f'models/{stem}.cbm'))
                        eta=np.empty(len(target))
                        for a in range(0,len(target),cfg['limits']['prediction_chunk_size']):
                            budget.check();b=a+cfg['limits']['prediction_chunk_size']
                            eta[a:b]=model.predict(target.iloc[a:b][names],prediction_type='RawFormulaVal',thread_count=2)
                        if not np.isfinite(eta).all():raise ValueError('Nonfinite target prediction')
                        novel=np.any((target.iloc[positions][names].to_numpy()<training[names].min().to_numpy())|
                                     (target.iloc[positions][names].to_numpy()>training[names].max().to_numpy()),axis=1)
                        result,scores,captures=held_diagnostics(eta[positions],eta,split['target_area_m2'],blocks,novel)
                        result.update(species=species,fold=fold,sampling=arm,family=family,features=len(names),
                            train_cells=n,held_cells=len(positions),held_blocks=len(np.unique(blocks)),
                            target_cells=len(target),stem=stem,model_source_run=args.run_id,reused=False)
                        metrics.append(result)
                        np.savez_compressed(out/f'predictions/{stem}.npz',target_raw_eta=eta,held_raw_eta=eta[positions],
                                            held_scores=scores,audit_indices=audit_idx,audit_raw_eta=eta[audit_idx])
                        details.update(species=species,fold=fold,sampling=arm,family=family,feature_names=names,
                            candidate=cfg['model'],seed=seed,source_split=str(split_path.relative_to(ROOT)),
                            source_split_sha256=inputs[str(split_path.relative_to(ROOT))],borders=f'borders/{prefix}_{family}.tsv',
                            borders_sha256=sha256(borders),parent_monthly_borders_sha256=sha256(parent/f'borders/{prefix}_monthly.tsv'),
                            target_log_normalizer=result['target_log_normalizer'],held_area_captures=captures,
                            background_cells_area_inclusion_weights_unchanged=True,all_original_presence_rows_retained=True,
                            frozen_weight_vector_sha256=sha256(wp) if arm=='soft_geo_10km' else None,
                            identity_only_weighting_change=True)
                        write_json(out/f'models/{stem}_metadata.json',details)
                        fits.append(dict(species=species,fold=fold,sampling=arm,family=family,stem=stem,status='completed',
                                         seconds=time.monotonic()-tick,training_metric_error=details['training_metric_error']))
                del training,target,eta
                flush()
                if sum(p.stat().st_size for p in out.rglob('*') if p.is_file())>cfg['limits']['output_gib']*1024**3:
                    raise RuntimeError('Declared output limit reached')
                print(json.dumps(dict(species=species,fold=fold,new_fits=record['empirical_fits'])),flush=True)
        if len(metrics)!=126 or len(fits)!=84 or len(reuses)!=42:
            raise ValueError('Incomplete matched design')
        for rel,value in stats.items():
            st=(ROOT/rel).stat()
            if (st.st_size,st.st_mtime_ns)!=value:raise ValueError('Input modified during run '+rel)
        record.update(status='completed',finished_utc=utcnow());flush()
        record['output_hashes']={str(p.relative_to(out)):sha256(p) for p in out.rglob('*') if p.is_file() and p.name!='run.json'}
        write_json(out/'run.json',record)
        print(json.dumps(dict(status='completed',new_fits=84,reused_monthly=42,cpu_seconds=record['cpu_seconds'])),flush=True)
    except Exception as error:
        record.update(status='failed_partial_results_preserved',error=str(error),traceback=traceback.format_exc());flush();raise


if __name__=='__main__':main()
