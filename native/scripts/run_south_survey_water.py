#!/usr/bin/env python3
"""Fixed terrestrial / mapped-water survey-history comparison; no locations inferred."""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(name, '1')

import argparse
import json
import resource
import time
import traceback
from xml.etree import ElementTree as ET

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import logit

from spea.provenance import ROOT, sha256, utcnow, write_json
from spea.survey_history import annual_usgs_history
from spea.survey_water import (
    HabitatTransform, history_terms, numerical_hessian, objective, observation_design,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', required=True)
    parser.add_argument('--run-id', required=True)
    args = parser.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())
    out = ROOT / 'runs' / args.run_id
    out.mkdir(exist_ok=False); (out / 'fits').mkdir()
    write_json(out / 'effective_config.json', cfg)
    cpu = time.process_time()
    checks = []
    state = dict(status='preflight', created_utc=utcnow(), fits=0, optimizer_calls=0,
                 occurrence_changes=0, spatial_crosswalks_inferred=0, hindcast_fits=0)

    def check(name, passed):
        checks.append(dict(name=name, passed=bool(passed)))
        if not passed:
            raise ValueError(name)

    def flush():
        state.update(cpu_seconds=time.process_time() - cpu,
                     peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write_json(out / 'run.json', state)
        write_json(out / 'checks.json', dict(checks=checks, passed=sum(c['passed'] for c in checks),
                                           failed=sum(not c['passed'] for c in checks)))

    try:
        for path, digest in cfg['input_hashes'].items():
            check('frozen:' + path, sha256(ROOT / path) == digest)
        survey = pd.read_csv(ROOT / cfg['survey_table']).sort_values('cellID').reset_index(drop=True)
        effort = pd.read_csv(ROOT / cfg['effort_table'])
        metadata = ET.parse(ROOT / cfg['metadata']).getroot()
        for year in cfg['years']:
            name = f'eff.{year % 100:02}'
            definitions = [a for a in metadata.findall('.//attr') if a.findtext('attrlabl') == name]
            check(name + ':source_NA_no_trapping', len(definitions) == 1 and
                  any(e.findtext('edomv') == 'NA' and
                      e.findtext('edomvd') == 'No trapping in that cell in this year.'
                      for e in definitions[0].findall('.//edom')))
        history = annual_usgs_history(survey, effort, cfg['years'])
        ids = survey.cellID.to_numpy()
        years = cfg['years']

        def matrix(field):
            return history.pivot(index='source_cell_id', columns='year', values=field).loc[ids, years].to_numpy()

        observed = matrix('observed').astype(bool)
        y = matrix('encounter_outcome')
        obs_values = np.stack([matrix(c) for c in cfg['observation_features']], axis=2)
        fields = cfg['terrestrial_features'] + cfg['water_features']
        values = survey[fields].to_numpy(float)
        check('all_habitat_complete', np.isfinite(values).all())
        check('water_nonnegative', (values[:, -2:] >= 0).all())
        check('complete376_cells', len(ids) == 376 and np.all(observed.sum(axis=1) > 0))
        check('source8272_years_2133_surveyed_150_detections', len(history) == 8272 and
              observed.sum() == 2133 and np.nansum(y) == 150)
        history.to_csv(out / 'local_annual_histories.csv.gz', index=False)
        inventory = survey[['cellID', 'block', 'n']].copy()
        inventory['surveyed_years'] = observed.sum(axis=1)
        inventory['detection_years'] = np.nansum(y, axis=1).astype(int)
        inventory['detected'] = inventory.detection_years > 0
        inventory.to_csv(out / 'local_cell_inventory.csv', index=False)
        coverage = inventory.groupby('block').agg(
            cells=('cellID', 'size'), detected_cells=('detected', 'sum'),
            surveyed_years=('surveyed_years', 'sum'), detection_years=('detection_years', 'sum'))
        coverage.to_csv(out / 'block_coverage.csv')
        scores, cell_scores, summaries = [], [], []
        for block in [*cfg['blocks'], 'full']:
            train = np.arange(len(ids)) if block == 'full' else np.flatnonzero(survey.block != block)
            held = np.arange(len(ids)) if block == 'full' else np.flatnonzero(survey.block == block)
            if block != 'full':
                check(str(block) + ':disjoint_source_cells', not set(ids[train]) & set(ids[held]))
            w, obs_mean, obs_scale = observation_design(obs_values, observed, train)
            for arm, water in cfg['arms'].items():
                raw = values.copy() if water else values[:, :len(cfg['terrestrial_features'])].copy()
                if water:
                    raw[:, -2:] = np.log1p(raw[:, -2:])
                transform = HabitatTransform.fit(raw[train], len(cfg['terrestrial_features']))
                x = transform.apply(raw)
                feature_names = ['intercept'] + cfg['terrestrial_features'] + (cfg['water_features'] if water else [])
                feature_names += [c + ':quadratic' for c in cfg['terrestrial_features']]
                fun = lambda theta: objective(theta, x[train], w[train], y[train], observed[train], cfg['ridge_precision'])
                check('model_fit_cap', state['fits'] < cfg['maximum_fits'])
                state['fits'] += 1
                starts = []
                for start in cfg['starts']:
                    check('optimizer_cap', state['optimizer_calls'] < cfg['maximum_optimizer_calls'])
                    state['optimizer_calls'] += 1; flush()
                    theta0 = np.zeros(x.shape[1] + w.shape[2])
                    theta0[0] = logit(start['psi']); theta0[x.shape[1]] = logit(start['p'])
                    fit = minimize(fun, theta0, jac=True, method='L-BFGS-B', options=cfg['optimizer_options'])
                    starts.append(dict(theta=fit.x.tolist(), value=float(fit.fun), success=bool(fit.success),
                                       message=str(fit.message), iterations=int(fit.nit),
                                       gradient_inf=float(np.max(np.abs(fit.jac)))))
                    check('CPU_budget', time.process_time() - cpu < cfg['maximum_cpu_seconds'])
                eligible = [i for i, s in enumerate(starts) if s['success'] and np.isfinite(s['value'])]
                check(f'{arm}:{block}:converged_start', bool(eligible))
                chosen = min(eligible, key=lambda i: starts[i]['value'])
                theta = np.asarray(starts[chosen]['theta'])
                check(f'{arm}:{block}:stationary', starts[chosen]['gradient_inf'] < cfg['maximum_gradient_inf'])
                eigen = np.linalg.eigvalsh(numerical_hessian(theta, fun))
                train_terms = history_terms(theta, x[train], w[train], y[train], observed[train])
                record = dict(arm=arm, held_block=block, training_cells=ids[train].tolist(),
                              evaluated_cells=ids[held].tolist(), habitat_features=feature_names,
                              observation_features=['intercept'] + cfg['observation_features'],
                              habitat_transform=transform.as_dict(), observation_mean=obs_mean.tolist(),
                              observation_scale=obs_scale.tolist(), theta=theta.tolist(),
                              selected_start=chosen, starts=starts,
                              start_objective_spread=float(np.ptp([s['value'] for s in starts])),
                              penalized_hessian_eigenvalues=eigen.tolist(),
                              training_unpenalized_nll=float(-train_terms['log_history'].sum()))
                write_json(out / 'fits' / f'{arm}_{block}.json', record)
                terms = history_terms(theta, x[held], w[held], y[held], observed[held])
                check(f'{arm}:{block}:finite_predictions', all(np.isfinite(terms[key]).all()
                      for key in ['log_history', 'log_any', 'log_pattern_given_any', 'q', 'psi']))
                role = 'training_diagnostic' if block == 'full' else 'held_block'
                for j, i in enumerate(held):
                    cell_scores.append(dict(arm=arm, held_block=block, role=role, source_cell_id=int(ids[i]),
                        surveyed_years=int(observed[i].sum()), detection_years=int(np.nansum(y[i])),
                        any_detection=bool(terms['positive'][j]),
                        joint_log_score=float(terms['log_history'][j]), any_log_score=float(terms['log_any'][j]),
                        conditional_pattern_log_score=float(terms['log_pattern_given_any'][j]),
                        predicted_any_encounter=float(terms['q'][j]), modeled_static_state=float(terms['psi'][j])))
                if block != 'full':
                    scores.append(dict(arm=arm, held_block=block, training_cells=len(train), cells=len(held),
                        surveyed_years=int(observed[held].sum()), detection_years=int(np.nansum(y[held])),
                        detected_cells=int(terms['positive'].sum()),
                        mean_joint_log_score=float(terms['log_history'].mean()),
                        mean_any_log_score=float(terms['log_any'].mean()),
                        mean_conditional_pattern_log_score=float(terms['log_pattern_given_any'].mean()),
                        predicted_detected_cells=float(terms['q'].sum())))
                summaries.append({k: v for k, v in record.items() if k not in
                                  ['training_cells', 'evaluated_cells', 'theta', 'habitat_transform', 'starts']})
                print(json.dumps(dict(completed=f'{arm}_{block}', fits=state['fits'],
                                      optimizer_calls=state['optimizer_calls'])), flush=True)
        pd.DataFrame(scores).to_csv(out / 'scores.csv', index=False)
        pd.DataFrame(cell_scores).to_csv(out / 'local_cell_scores.csv', index=False)
        write_json(out / 'fit_summary.json', summaries)
        for path, digest in cfg['input_hashes'].items():
            check('preserved:' + path, sha256(ROOT / path) == digest)
        check('exact_fit_budget', state['fits'] == cfg['maximum_fits'] and
              state['optimizer_calls'] == cfg['maximum_optimizer_calls'])
        state.update(status='completed_developmental_static_survey_history_comparison', completed_utc=utcnow(),
                     source_cells=376, surveyed_years=2133, no_trapping_years=6139, detection_years=150,
                     held_block_fits=10, descriptive_full_fits=2, input_hashes=cfg['input_hashes'])
        flush()
        state['output_hashes'] = {str(p.relative_to(out)): sha256(p) for p in sorted(out.rglob('*'))
                                  if p.is_file() and p.name != 'run.json'}
        write_json(out / 'run.json', state)
        print(json.dumps({k: v for k, v in state.items() if k not in ['input_hashes', 'output_hashes']}))
    except Exception as exc:
        state.update(status='failed', error=str(exc), traceback=traceback.format_exc()); flush(); raise


if __name__ == '__main__':
    main()
