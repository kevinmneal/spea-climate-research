#!/usr/bin/env python3
"""Frozen reciprocal N/S transfer on existing spatial/source splits.

Reuse ordinary CatBoost fitting helpers. Neither group identity nor geography is
a feature. All target normalizers are separate by group and use complete native
held areas. Equal training group mass is explicit for both classifier classes.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import resource
import shutil
import sys
import time
import traceback
import zipfile
from pathlib import Path

from run_genus_environment_comparison import (
    Budget, BudgetExceeded, digest_array, fit_model, predict, save_model,
    training_reference,
)
import catboost
import numpy as np
import pandas as pd
from scipy.special import logsumexp
from spea.environment_benchmark import NativeEnvironment
from spea.provenance import ROOT, sha256, utcnow, write_json


def group_balanced_weights(pgroup, qgroup, area_weight):
    """Half total mass per class, equal group mass within either class."""
    pg = np.asarray(pgroup); qg = np.asarray(qgroup); w = np.asarray(area_weight, dtype=float)
    groups = np.unique(pg)
    if not len(pg) or not np.array_equal(groups, np.unique(qg)):
        raise ValueError('Every modeled group needs presence and background support')
    if len(qg) != len(w) or not np.isfinite(w).all() or np.any(w <= 0):
        raise ValueError('Invalid area/inclusion-probability weights')
    half = (len(pg)+len(qg))/2.; pm = np.empty(len(pg)); qm = np.empty(len(qg))
    for group in groups:
        pp, qq = pg == group, qg == group
        pm[pp] = half/len(groups)/pp.sum()
        scaled = w[qq]/w[qq].max(); scaled /= scaled.sum()
        qm[qq] = half/len(groups)*scaled
    return np.r_[pm, qm]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-id', required=True)
    parser.add_argument('--config', default='config/hammondii_reciprocal_transfer.json')
    parser.add_argument('--plan-only', action='store_true')
    args = parser.parse_args()
    if Path(args.run_id).name != args.run_id: raise ValueError('Run ID must be one path component')
    config_path = ROOT / args.config; cfg = json.loads(config_path.read_text())
    if args.plan_only: print(json.dumps(cfg, indent=2)); return
    if catboost.__version__ != '1.2.10': raise ValueError('Pinned CatBoost 1.2.10 required')
    slots_count = len(cfg['cases'])*len(cfg['folds'])*len(cfg['families'])*len(cfg['train_modes'])
    if slots_count != 36 or cfg['budget']['maximum_fits'] != 36:
        raise ValueError('Unexpected frozen comparison size')
    out = ROOT / 'runs' / args.run_id; out.mkdir(exist_ok=False)
    for folder in ['models', 'predictions', 'training_design', 'source_splits']: (out/folder).mkdir()
    previous, environment, preparation = [ROOT/'runs'/cfg[k] for k in
        ['previous_run', 'environment_run', 'preparation_run']]
    started = time.monotonic(); budget = Budget(cfg['budget'])
    inputs, checks, metrics, cellscores, audit, coverage, feature_novelty = {}, [], [], [], [], [], []
    slots, designs, splits = [], {}, {}; model_cache = {}
    record = dict(status='preparing', created_utc=utcnow(), command_argv=sys.argv, config=cfg,
        input_hashes=inputs, fit_attempts=0, completed_slots=0, reused_slots=0,
        full_data_fits=0, parameter_searches=0, formal_significance_tests=0, occurrence_changes=0)

    def track(path, expected=None):
        digest = sha256(path)
        if expected is not None and digest != expected: raise ValueError('Changed input: '+str(path))
        inputs[str(path.relative_to(ROOT))] = digest
        return path

    def check(name, passed):
        checks.append(dict(name=name, passed=bool(passed)))
        if not passed: raise ValueError(name)

    def disk_bytes(): return sum(p.stat().st_size for p in out.rglob('*') if p.is_file())

    def contrasts():
        frame = pd.DataFrame(metrics); rows = []
        if frame.empty: return rows
        keys = ['case', 'fold', 'family', 'target_cohort', 'group']
        for identity, part in frame.groupby(keys, sort=False):
            by_mode = {r.train_mode: r for r in part.itertuples()}
            own = 'north' if identity[-1] == 1 else 'south'; other = 'south' if own == 'north' else 'north'
            if own not in by_mode: continue
            for name, low in [('within_minus_cross', other), ('within_minus_pooled', 'pooled')]:
                if low not in by_mode: continue
                a, b = by_mode[own], by_mode[low]
                rows.append(dict(zip(keys, identity)) | dict(contrast=name, within_mode=own,
                    comparator_mode=low, within_slot_id=a.slot_id, comparator_slot_id=b.slot_id,
                    cells=a.cells, blocks=a.blocks, cell_gain=a.cell_score-b.cell_score,
                    block_gain=a.block_score-b.block_score))
        return rows

    def flush():
        for name, rows in [('metrics.csv', metrics), ('cell_scores.csv.gz', cellscores),
                ('fit_audit.csv', audit), ('coverage.csv', coverage),
                ('held_feature_novelty.csv', feature_novelty), ('paired_contrasts.csv', contrasts())]:
            if rows: pd.DataFrame(rows).to_csv(out/name, index=False)
        if metrics:
            pd.DataFrame(metrics).groupby(['case', 'family', 'train_mode', 'target_cohort', 'group'],
                sort=False).agg(folds=('fold', 'count'), cell_score=('cell_score', 'mean'),
                block_score=('block_score', 'mean')).to_csv(out/'mean_scores.csv')
        record.update(elapsed_seconds=time.monotonic()-started, cpu_seconds=budget.cpu(),
            peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, disk_bytes=disk_bytes())
        write_json(out/'run.json', record)

    write_json(out/'run.json', record)
    try:
        pm = json.loads(track(previous/'run.json').read_text())
        em = json.loads(track(environment/'run.json').read_text())
        prep = json.loads(track(preparation/'run.json').read_text())
        check('source_runs_completed', all(m['status'] == 'completed' for m in [pm, em, prep]))
        grid_path = track(previous/'analysis_grid.npz', pm['output_hashes']['analysis_grid.npz'])
        shutil.copyfile(grid_path, out/'analysis_grid.npz')
        with np.load(grid_path) as saved: grid = {k: saved[k] for k in saved.files}
        native, group, area, block = [grid[k] for k in ['native', 'group', 'area_m2', 'block50']]
        check('sorted_unique_native_ids', np.all(np.diff(native) > 0))
        check('two_valid_groups', np.array_equal(np.unique(group), [-1, 1]))
        check('positive_finite_areas', np.isfinite(area).all() and (area > 0).all())
        env = json.loads(track(environment/'registry.json', em['output_hashes']['registry.json']).read_text())
        base = json.loads(track(preparation/'feature_registry.json', prep['output_hashes']['feature_registry.json']).read_text())
        check('matching_native_grid', env['window_grid'] == base['window_grid'])
        check('native_ids_inside_grid', np.all(native >= 0) and np.all(native < np.prod(env['window_grid']['shape'])))
        check('broad50_exact', cfg['families']['broad50'] == env['families']['broad'])
        terrain = ['terrain_elevation_m', 'terrain_slope_degrees', 'terrain_relief_m', 'terrain_missing']
        check('climate46_exact', cfg['families']['climate46'] == [k for k in cfg['families']['broad50'] if k not in terrain])
        check('expected_feature_counts', len(cfg['families']['climate46']) == 46 and len(cfg['families']['broad50']) == 50)
        forbidden = {'group', 'latitude', 'longitude', 'distance_m', 'block50', 'cec_level2'}
        check('no_group_or_location_features', not forbidden.intersection(cfg['families']['broad50']))
        for field in env['fields']:
            if 'physical_path' in field: track(ROOT/field['physical_path'], em['input_hashes'][field['physical_path']])
        for path in [config_path, Path(__file__), ROOT/'scripts/run_genus_environment_comparison.py',
                ROOT/'src/spea/environment_benchmark.py', ROOT/'src/spea/native_benchmark.py',
                ROOT/'src/spea/genus_model.py', ROOT/'src/spea/genus_spatial.py', ROOT/'src/spea/provenance.py',
                ROOT/'pyproject.toml', ROOT/'uv.lock']:
            track(path)
        reader = NativeEnvironment(ROOT, env['window_grid']['shape'], env['fields'])
        # Prepare every design, alias and source link before any empirical fit.
        fingerprints = {}; evaluation_by_fold = {}
        for case in cfg['cases']:
            for fold in cfg['folds']:
                split_name = f'{case}_{fold}.npz'
                path = track(previous/'splits'/split_name, pm['output_hashes']['splits/'+split_name])
                shutil.copyfile(path, out/'source_splits'/split_name)
                with np.load(path) as saved: split = {k: saved[k] for k in saved.files}
                splits[(case, fold)] = split
                p, q, held, target = [split[k] for k in ['p', 'q', 'held', 'target']]
                check(case+str(fold)+':valid_analysis_indices', all(np.all(v >= 0) and np.all(v < len(native)) for v in [p, q, held, target]))
                check(case+str(fold)+':presence_in_quadrature', np.isin(p, q).all())
                check(case+str(fold)+':target_disjoint', not np.isin(q, target).any() and np.isin(held, target).all())
                check(case+str(fold)+':sorted_target', np.all(np.diff(target) > 0))
                if fold in evaluation_by_fold:
                    old = evaluation_by_fold[fold]
                    check(case+str(fold)+':identical_held_source_cases', np.array_equal(held, old['held']))
                    check(case+str(fold)+':identical_targets_source_cases', np.array_equal(target, old['target']))
                else: evaluation_by_fold[fold] = split
                for mode_index, mode in enumerate(cfg['train_modes']):
                    selected_group = {'north': 1, 'south': -1}.get(mode)
                    pp = p if selected_group is None else p[group[p] == selected_group]
                    qmask = np.ones(len(q), dtype=bool) if selected_group is None else group[q] == selected_group
                    qq, qw, pi = q[qmask], split['w'][qmask], split['pi'][qmask]
                    weight = group_balanced_weights(group[pp], group[qq], qw)
                    check(case+str(fold)+mode+':selected_presence_in_quadrature', np.isin(pp, qq).all())
                    check(case+str(fold)+mode+':positive_weights', np.isfinite(weight).all() and (weight > 0).all())
                    check(case+str(fold)+mode+':weight_mass', np.isclose(weight[:len(pp)].sum(), weight[len(pp):].sum()) and np.isclose(weight.mean(), 1.))
                    design_name = f'{case}_fold{fold}_{mode}.npz'
                    design = dict(train_native=native[pp], q_native=native[qq], train_analysis_index=pp,
                        q_analysis_index=qq, train_group=group[pp], q_group=group[qq], q_weights_m2=qw,
                        q_inclusion_probability=pi, class_object_weights=weight)
                    designs[(case, fold, mode)] = design
                    np.savez_compressed(out/'training_design'/design_name, **design)
                    seed = cfg['seed']+100*fold+mode_index
                    for family, keys in cfg['families'].items():
                        slot_id = f'{case}_fold{fold}_{family}_{mode}'
                        identity = dict(seed=seed, parameters=cfg['parameters'], feature_names=keys,
                            train_native=digest_array(native[pp]), q_native=digest_array(native[qq]),
                            class_object_weights=digest_array(weight),
                            held_native=digest_array(native[held]), target_native=digest_array(native[target]),
                            target_area_m2=digest_array(area[target]), target_group=digest_array(group[target]),
                            environment_registry_sha256=inputs[str((environment/'registry.json').relative_to(ROOT))])
                        fingerprint = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()
                        canonical = fingerprints.setdefault(fingerprint, slot_id)
                        slots.append(dict(slot_id=slot_id, case=case, fold=fold, family=family, train_mode=mode,
                            seed=seed, training_design_path='training_design/'+design_name, fingerprint=fingerprint,
                            model_id=canonical, reused_from=canonical if canonical != slot_id else '',
                            prediction_path='predictions/'+canonical+'.npz', model_path='models/'+canonical+'.cbm'))
                    for gg in [-1, 1]:
                        for role, ids in [('train', pp), ('quadrature', qq), ('held', held), ('target', target)]:
                            selected = ids[group[ids] == gg]
                            coverage.append(dict(case=case, fold=fold, train_mode=mode, group=gg, role=role,
                                cells=len(selected), blocks=len(np.unique(block[selected])),
                                without_S35_cells=int(grid['without_S35'][selected].sum()) if role in ['train', 'held'] else None))
        check('36_comparison_slots', len(slots) == 36)
        check('30_unique_fits', len(fingerprints) == cfg['budget']['expected_unique_fits'] == 30)
        check('six_exact_north_aliases', sum(bool(s['reused_from']) for s in slots) == 6 and all(s['train_mode'] == 'north' for s in slots if s['reused_from']))
        pd.DataFrame(slots).to_csv(out/'slot_registry.csv', index=False)
        write_json(out/'effective_config.json', cfg)
        write_json(out/'feature_registry.json', dict(window_grid=env['window_grid'], fields=env['fields'], families=cfg['families']))
        with zipfile.ZipFile(out/'source_snapshot.zip', 'w', zipfile.ZIP_DEFLATED) as archive:
            for relative in inputs:
                if relative.endswith(('.py', '.json', '.toml', '.lock')): archive.write(ROOT/relative, relative)
        flush()
        write_json(out/'prefit_receipt.json', dict(created_utc=utcnow(), input_hashes=dict(inputs),
            planned_slots=len(slots), unique_fits=len(fingerprints), aliases=6, empirical_fits=0,
            frozen_output_hashes={str(p.relative_to(out)): sha256(p) for p in out.rglob('*') if p.is_file() and p.name != 'run.json'}))
        record['status'] = 'running'; flush()
        for slot in slots:
            budget.check()
            if disk_bytes() >= cfg['budget']['max_output_bytes']: raise BudgetExceeded('Declared output disk budget exhausted')
            case, fold, family, mode = [slot[k] for k in ['case', 'fold', 'family', 'train_mode']]
            design = designs[(case, fold, mode)]; keys = cfg['families'][family]
            split = splits[(case, fold)]; held, target = split['held'], split['target']
            pids, qids, weight = [design[k] for k in ['train_native', 'q_native', 'class_object_weights']]
            entry = dict(slot, status='running', n_train=len(pids), n_quadrature=len(qids),
                empirical_fit=not bool(slot['reused_from']), class_object_weights_sha256=digest_array(weight))
            audit.append(entry); flush(); t = time.monotonic()
            if not slot['reused_from']:
                if record['fit_attempts'] >= cfg['budget']['expected_unique_fits']: raise BudgetExceeded('Declared unique-fit budget exhausted')
                record['fit_attempts'] += 1; flush()
                train = reader.read(np.r_[pids, qids], keys); reference = training_reference(train, len(pids), keys)
                model = fit_model(train, len(pids), weight, keys, cfg, slot['seed'], budget)
                save_model(model, slot['model_id'], out, keys, reference, weight, slot['seed'])
                fit_seconds = time.monotonic()-t
                prediction = predict(model, native[target], keys, reference, reader, cfg['budget']['prediction_chunk_size'], budget)
                raw = prediction['raw_eta']; normalized = np.empty(len(target)); normalizers = []
                for gg in [-1, 1]:
                    use = group[target] == gg; weights = area[target[use]]
                    norm = float(logsumexp(raw[use]+np.log(weights))-np.log(weights.sum()))
                    normalized[use] = raw[use]-norm; normalizers.append(norm)
                    check(slot['slot_id']+f':normalizer_{gg}', abs(logsumexp(normalized[use]+np.log(weights))-np.log(weights.sum())) < 1e-10)
                held_at = np.searchsorted(target, held)
                check(slot['slot_id']+':held_target_identity', np.array_equal(target[held_at], held))
                held_frame = reader.read(native[held], keys)
                held_raw = model.predict(held_frame, prediction_type='RawFormulaVal', thread_count=2)
                check(slot['slot_id']+':held_prediction_replay', np.array_equal(held_raw, raw[held_at]))
                prediction_data = dict(target_native=native[target], target_analysis_index=target, target_group=group[target],
                    target_area_m2=area[target], target_raw_eta=raw, target_log_relative_density=normalized,
                    held_native=native[held], held_analysis_index=held, held_group=group[held], held_block50=block[held],
                    held_without_S35=grid['without_S35'][held], held_raw_eta=held_raw,
                    held_log_relative_density=normalized[held_at], normalizer_groups=np.asarray([-1, 1]),
                    group_log_normalizers=np.asarray(normalizers), novel_target=prediction['novel_numeric'],
                    novel_held=prediction['novel_numeric'][held_at], target_numeric_novel_count=prediction['numeric_novel_count'],
                    held_numeric_novel_count=prediction['numeric_novel_count'][held_at],
                    target_numeric_missing_count=prediction['numeric_missing_count'],
                    held_numeric_missing_count=prediction['numeric_missing_count'][held_at],
                    training_feature_min=reference['minimum'], training_feature_max=reference['maximum'],
                    training_numeric_feature_names=np.asarray(reference['numeric_keys']))
                np.savez_compressed(out/slot['prediction_path'], **prediction_data)
                model_cache[slot['model_id']] = dict(fit_seconds=fit_seconds)
                entry.update(fit_seconds=fit_seconds, tree_count=int(model.tree_count_))
                del train, model, prediction, raw, normalized, held_frame
            else:
                check(slot['slot_id']+':canonical_already_complete', slot['model_id'] in model_cache)
                with np.load(out/slot['prediction_path']) as saved: prediction_data = {k: saved[k] for k in saved.files}
                check(slot['slot_id']+':reuse_held_identity', np.array_equal(prediction_data['held_native'], native[held]))
                check(slot['slot_id']+':reuse_target_identity', np.array_equal(prediction_data['target_native'], native[target]))
                check(slot['slot_id']+':reuse_target_areas', np.array_equal(prediction_data['target_area_m2'], area[target]))
                entry.update(fit_seconds=0., tree_count=cfg['parameters']['iterations']); record['reused_slots'] += 1
            y = prediction_data['held_log_relative_density']
            held_values = reader.read(native[held], keys).to_numpy(float)
            lower, upper = prediction_data['training_feature_min'], prediction_data['training_feature_max']
            finite = np.isfinite(held_values); supported = np.isfinite(lower)&np.isfinite(upper)
            outside = finite&((held_values < lower)|(held_values > upper)|~supported)
            common = dict(slot_id=slot['slot_id'], model_id=slot['model_id'], case=case, fold=fold, family=family, train_mode=mode)
            for cohort, keep in [('all', np.ones(len(held), dtype=bool)), ('without_S35', grid['without_S35'][held])]:
                for gg in [-1, 1]:
                    use = keep&(group[held] == gg); yy, bb = y[use], block[held[use]]
                    check(slot['slot_id']+cohort+str(gg)+':nonempty_evaluation', len(yy) >= 2)
                    block_mean = pd.DataFrame(dict(score=yy, block=bb)).groupby('block').score.mean().mean()
                    target_use = group[target] == gg
                    metrics.append(common | dict(target_cohort=cohort, group=gg, cells=len(yy), blocks=len(np.unique(bb)),
                        cell_score=float(yy.mean()), block_score=float(block_mean),
                        target_area_km2=float(area[target[target_use]].sum()/1e6),
                        novel_cell_fraction=float(prediction_data['novel_held'][use].mean()),
                        novel_target_area_fraction=float(np.average(prediction_data['novel_target'][target_use], weights=area[target[target_use]])),
                        missing_target_area_fraction=float(np.average(prediction_data['target_numeric_missing_count'][target_use] > 0, weights=area[target[target_use]])),
                        missing_cell_fraction=float((prediction_data['held_numeric_missing_count'][use] > 0).mean())))
                    for j, key in enumerate(keys):
                        feature_novelty.append(common | dict(target_cohort=cohort, group=gg, feature=key,
                            held_cells=int(use.sum()), out_of_range_cells=int(outside[use, j].sum()),
                            missing_cells=int((~finite[use, j]).sum()),
                            training_min=float(lower[j]) if np.isfinite(lower[j]) else None,
                            training_max=float(upper[j]) if np.isfinite(upper[j]) else None))
            for i, ii in enumerate(held):
                cellscores.append(common | dict(native_flat=int(native[ii]), analysis_index=int(ii), group=int(group[ii]),
                    block50=int(block[ii]), without_S35=bool(grid['without_S35'][ii]), score=float(y[i]),
                    raw_eta=float(prediction_data['held_raw_eta'][i]), novel_numeric=bool(prediction_data['novel_held'][i]),
                    numeric_novel_count=int(prediction_data['held_numeric_novel_count'][i]),
                    numeric_missing_count=int(prediction_data['held_numeric_missing_count'][i])))
            entry.update(status='passed', total_seconds=time.monotonic()-t)
            record['completed_slots'] += 1; flush()
            print(json.dumps({k: entry[k] for k in ['slot_id', 'model_id', 'empirical_fit', 'status', 'fit_seconds', 'total_seconds']}), flush=True)
            del prediction_data, held_values
        check('30_completed_empirical_fits', record['fit_attempts'] == 30)
        check('36_completed_slots', record['completed_slots'] == 36 and record['reused_slots'] == 6)
        check('144_metrics', len(metrics) == 144); check('96_matched_contrasts', len(contrasts()) == 96)
        check('within_disk_budget', disk_bytes() < cfg['budget']['max_output_bytes'])
        for relative, digest in inputs.items(): check('unchanged:'+relative, sha256(ROOT/relative) == digest)
        record.update(status='completed', finished_utc=utcnow())
    except BaseException as error:
        if audit and audit[-1]['status'] == 'running': audit[-1].update(status='failed', error=str(error))
        record.update(status='budget_exhausted' if isinstance(error, BudgetExceeded) else 'interrupted_or_failed',
            error=str(error), traceback=traceback.format_exc(), finished_utc=utcnow())
        raise
    finally:
        flush(); write_json(out/'checks.json', dict(passed=sum(c['passed'] for c in checks),
            failed=sum(not c['passed'] for c in checks), checks=checks))
        record['output_hashes'] = {str(p.relative_to(out)): sha256(p) for p in out.rglob('*') if p.is_file() and p.name != 'run.json'}
        write_json(out/'run.json', record)
    print(json.dumps({k: record[k] for k in ['status', 'fit_attempts', 'completed_slots', 'reused_slots', 'cpu_seconds', 'disk_bytes']}), flush=True)


if __name__ == '__main__': main()
