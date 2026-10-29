#!/usr/bin/env python3
"""Freeze separate-species seasonal splits; inspect feasibility without fitting."""
from __future__ import annotations
import argparse
import json
import time
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from spea.provenance import ROOT, sha256, utcnow, write_json
from spea.genus_model import geographic_blocks, stratified_quadrature
from spea.genus_spatial import exact_pair_min_clearance


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config', default='config/hammondii_species_seasonal_comparison.json')
    args = ap.parse_args()
    cfgpath = ROOT / args.config
    cfg = json.loads(cfgpath.read_text())
    out = ROOT / 'runs' / cfg['species_split_run']
    out.mkdir(exist_ok=False)
    (out / 'splits').mkdir()
    start = time.monotonic()
    inputs, checks, designs, outer_designs, manifests = {}, [], [], [], {}
    record = dict(status='preparing', created_utc=utcnow(), config=cfg, model_fits=0,
                  canonical_occurrence_changes=0, scientific_gate_changes=0)

    def check(name, value):
        checks.append(dict(name=name, passed=bool(value)))
        if not value:
            raise ValueError(name)

    def track(path, expected=None):
        path = Path(path)
        actual = sha256(path)
        check('hash:' + str(path.relative_to(ROOT)), expected is None or actual == expected)
        inputs[str(path.relative_to(ROOT))] = actual
        return path

    def asset(run, name):
        directory = ROOT / 'runs' / run
        if run not in manifests:
            manifests[run] = json.loads(track(directory / 'run.json').read_text())
            check('source_completed:' + run, manifests[run]['status'] == 'completed')
        return track(directory / name, manifests[run]['output_hashes'][name])

    try:
        track(cfgpath)
        for path in [Path(__file__).resolve(), ROOT / 'src/spea/genus_model.py',
                     ROOT / 'src/spea/genus_spatial.py', ROOT / 'src/spea/provenance.py']:
            track(path)
        parent = cfg['parent_spatial_run']
        with np.load(asset(parent, 'analysis_grid.npz')) as source:
            grid = {key: source[key] for key in source.files}
        native = grid['native'].astype('int64')
        group = grid['group']
        check('sorted_native', np.all(np.diff(native) > 0))
        check('two_assigned_species', set(np.unique(group)) == {-1, 1})
        cells = pd.read_csv(asset(cfg['preparation_run'], 'cells.csv.gz'))
        registry = json.loads(asset(cfg['preparation_run'], 'feature_registry.json').read_text())
        rowarea = np.load(asset(cfg['preparation_run'], 'row_area_m2.npy'))
        shape = registry['window_grid']['shape']
        check('physical_grid_areas', np.array_equal(grid['area_m2'], rowarea[native // shape[1]]))
        h = cells[cells.working_taxon.eq('Spea hammondii')].copy()
        by_native = h.reset_index().set_index('native_flat')['index']
        cp = grid['presence_index']
        check('same_591_occurrences', len(h) == 591 and
              np.array_equal(np.sort(h.native_flat), np.sort(native[cp])))
        parent_cfg = manifests[parent]['config']
        assignments = pd.read_csv(asset(parent_cfg['support_run'], 'cell_geographic_assignments.csv'))
        assignment_idx = np.searchsorted(native, assignments.native_flat)
        check('inherited_assignment_identity', np.array_equal(native[assignment_idx], assignments.native_flat)
              and np.array_equal(group[assignment_idx], assignments.provisional_group.map({'North': 1, 'South': -1})))
        h['source_cell_row'] = h.index
        h['analysis_species'] = np.where(group[np.searchsorted(native, h.native_flat)] == 1,
                                         'hammondii_north', 'spea_south')
        h['assignment_basis'] = 'Inherited geographic N/S assignment; not a new genotype determination'
        h['assignment_source_run'] = parent_cfg['support_run']
        h.to_csv(out / 'species_cells.csv', index=False)
        blocks = geographic_blocks(native, shape, registry['window_grid']['transform'],
                                   cfg['quadrature']['stratum_width_m'])
        splits = {}
        for fold in cfg['folds']:
            with np.load(asset(parent, f'splits/{cfg["parent_case"]}_{fold}.npz')) as saved:
                splits[fold] = {key: saved[key] for key in saved.files}

        def validate_and_save(prefix, species, p, held, domain, target, q, w, pi, inner):
            rows = by_native.loc[native[p]].to_numpy('int64')
            heldrows = by_native.loc[native[held]].to_numpy('int64')
            check(prefix + ':single_species', all(np.all(group[a] == cfg['species_group_codes'][species])
                                                 for a in [p, held, domain, target, q]))
            check(prefix + ':membership', np.isin(p, domain).all() and np.isin(q, domain).all()
                  and np.isin(held, target).all() and not np.isin(domain, target).any())
            check(prefix + ':event_inclusive', np.all(np.diff(native[q]) > 0) and np.isin(p, q).all())
            check(prefix + ':area_probability', np.array_equal(w, grid['area_m2'][q] / pi)
                  and np.all((pi > 0) & (pi <= 1)))
            gap = exact_pair_min_clearance(cells.loc[rows, 'longitude'].to_numpy(),
                cells.loc[rows, 'latitude'].to_numpy(), cells.loc[rows, 'radius_unknown1km_m'].to_numpy(),
                cells.loc[heldrows, 'longitude'].to_numpy(), cells.loc[heldrows, 'latitude'].to_numpy(),
                cells.loc[heldrows, 'radius_unknown1km_m'].to_numpy())
            check(prefix + ':support_gap30km', len(gap) > 0 and np.all(gap >= 30000 - 1e-6))
            check(prefix + ':block_identity', np.array_equal(grid['block50'][held], cells.loc[heldrows, 'block50']))
            key = 'validation' if inner else 'test'
            payload = dict(train_native=native[p], train_cell_rows=rows, q_native=native[q],
                q_weights_m2=w, q_inclusion_probability=pi, target_native=native[target],
                target_area_m2=grid['area_m2'][target], eligible_native=native[domain],
                **{key + '_native': native[held], key + '_cell_rows': heldrows})
            np.savez_compressed(out / 'splits' / (prefix + '.npz'), **payload)
            return float(gap.min())

        actual = {}
        for species in cfg['species']:
            code = cfg['species_group_codes'][species]
            keep = lambda indices: np.asarray(indices)[group[indices] == code]
            for outer in cfg['folds']:
                source = splits[outer]
                p, held, domain, target = [keep(source[k]) for k in ['p', 'held', 'domain', 'target']]
                qkeep = group[source['q']] == code
                q, w, pi = source['q'][qkeep], source['w'][qkeep], source['pi'][qkeep]
                prefix = f'{species}_fold{outer}'
                gap = validate_and_save(prefix, species, p, held, domain, target, q, w, pi, False)
                outer_designs.append(dict(species=species, fold=outer, train_cells=len(p), held_cells=len(held),
                    held_blocks=len(np.unique(grid['block50'][held])), quadrature_cells=len(q),
                    target_cells=len(target), target_area_km2=grid['area_m2'][target].sum()/1e6,
                    min_support_gap_m=gap))
                for inner in cfg['folds']:
                    if inner == outer:
                        continue
                    other = splits[inner]
                    ip = keep(np.intersect1d(source['p'], other['p']))
                    vp = keep(np.intersect1d(source['p'], other['held']))
                    eligible = keep(np.intersect1d(source['domain'], other['domain']))
                    target_i = keep(np.intersect1d(source['domain'], other['target']))
                    prefix = f'{species}_outer{outer}_inner{inner}'
                    check(prefix + ':outer_holdout_excluded', not np.isin(np.r_[eligible, target_i], source['target']).any()
                          and not np.isin(np.r_[ip, vp], source['held']).any())
                    check(prefix + ':nested_domain_membership', np.isin(ip, eligible).all() and np.isin(vp, target_i).all())
                    nblocks = len(np.unique(grid['block50'][vp]))
                    use = len(ip) >= cfg['minimum_inner_training_cells'] and nblocks >= cfg['minimum_inner_validation_blocks']
                    design = dict(species=species, prefix=prefix, outer_fold=outer, inner_fold=inner,
                        train_cells=len(ip), validation_cells=len(vp), validation_blocks=nblocks, eligible=bool(use),
                        exclusion_reason='' if use else 'fewer_than_20_training_cells',
                        target_cells=len(target_i), min_support_gap_m=None)
                    if use:
                        seed = cfg['quadrature']['seed'] + 1000*cfg['species'].index(species) + 10*outer + inner
                        iq, iw, ipi, _, _ = stratified_quadrature(native[eligible], native[ip], blocks[eligible],
                            grid['area_m2'][eligible], cfg['quadrature']['nonpresence_per_stratum'], seed)
                        qi = np.searchsorted(native, iq)
                        design['min_support_gap_m'] = validate_and_save(prefix, species, ip, vp, eligible,
                            target_i, qi, iw, ipi, True)
                        design['quadrature_seed'] = seed
                    designs.append(design)
            actual[species] = [sum(d['eligible'] and d['outer_fold'] == f and d['species'] == species
                                   for d in designs) for f in cfg['folds']]
            check(species + ':expected_feasibility', actual[species] == cfg['expected_inner_folds'][species])
        check('120_inner_slots', sum(map(sum, actual.values()))*len(cfg['families'])*len(cfg['candidates']) == 120)
        pd.DataFrame(designs).to_csv(out / 'inner_design.csv', index=False)
        pd.DataFrame(outer_designs).to_csv(out / 'outer_design.csv', index=False)
        write_json(out / 'effective_config.json', cfg)
        with zipfile.ZipFile(out / 'source_snapshot.zip', 'w', zipfile.ZIP_DEFLATED) as z:
            for rel in inputs:
                if rel.endswith(('.py', '.json')):
                    z.write(ROOT / rel, rel)
        for rel, expected in inputs.items():
            check('unchanged:' + rel, sha256(ROOT / rel) == expected)
        record.update(status='completed', finished_utc=utcnow(), eligible_inner_folds=actual,
                      held_cells=sum(d['held_cells'] for d in outer_designs))
    except BaseException as error:
        record.update(status='failed', error=str(error))
        raise
    finally:
        write_json(out / 'checks.json', dict(passed=sum(c['passed'] for c in checks),
                   failed=sum(not c['passed'] for c in checks), checks=checks))
        record.update(input_hashes=inputs, elapsed_seconds=time.monotonic()-start,
                      output_hashes={str(p.relative_to(out)): sha256(p) for p in out.rglob('*')
                                     if p.is_file() and p.name != 'run.json'})
        write_json(out / 'run.json', record)
    print(json.dumps(dict(status=record['status'], checks=len(checks), eligible_inner_folds=actual)))


if __name__ == '__main__':
    main()
