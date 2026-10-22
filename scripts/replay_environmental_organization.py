#!/usr/bin/env python3
"""Audit the saved environmental-organization arithmetic using only Python's standard library."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
TOLERANCE = 1e-10
GENERATORS = ('uniform', 'shared', 'smooth')
MODELS = ('shared', 'smooth')
FEATURES = ('mean_temperature_c', 'log1p_annual_precipitation_mm', 'temperature_range_c',
            'precipitation_cv', 'precipitation_harmonic_cos', 'precipitation_harmonic_sin')
DATA_KEYS = ('cohort', 'generator', 'replicate', 'selection')
FIELDS = {
    'organization_statistics.csv': DATA_KEYS + ('statistic', 'held_blocks', 'small_component_folds'),
    'mixture_fold_scores.csv': DATA_KEYS + ('fold', 'train_blocks', 'held_blocks',
        'single_mean_log_density', 'mixture_mean_log_density', 'gain_sum', 'gain_mean',
        'smallest_component_mass', 'small_component', 'single_iterations', 'mixture_iterations'),
    'envelope_summary.csv': ('generator', 'selection', 'observed', 'simulated_05', 'simulated_median',
        'simulated_95', 'simulations_at_or_above_observed', 'simulations', 'observed_above_95'),
    'association_fold_scores.csv': ('fold', 'model', 'held_blocks', 'held_cells', 'mean_log_density'),
    'association_summary.csv': ('model', 'mean_held_log_density', 'held_blocks', 'held_cells'),
    'feature_variance.csv': ('feature', 'within_block_variance_fraction', 'within_block_variance_percent'),
}


def read_csv(path):
    with path.open(newline='', encoding='utf-8') as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != FIELDS[path.name]:
            raise ValueError(f'Unexpected columns in {path.name}')
        rows = list(reader)
    if not rows or any(None in row or any(v is None or v == '' for v in row.values()) for row in rows):
        raise ValueError(f'Empty or incomplete rows in {path.name}')
    return rows


def index(rows, keys):
    result = {tuple(str(row[k]) for k in keys): row for row in rows}
    if len(result) != len(rows):
        raise ValueError(f'Duplicate keys: {keys}')
    return result


def number(value):
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f'Invalid number: {value}') from error
    if not math.isfinite(result):
        raise ValueError(f'Nonfinite number: {value}')
    return result


def integer(value, minimum=0, maximum=None):
    result = number(value)
    if not result.is_integer() or result < minimum or (maximum is not None and result > maximum):
        raise ValueError(f'Invalid integer count: {value}')
    return int(result)


def boolean(value):
    if value not in ('True', 'False', True, False):
        raise ValueError(f'Invalid Boolean: {value}')
    return value in ('True', True)


def near(left, right):
    if not math.isclose(number(left), number(right), rel_tol=0, abs_tol=TOLERANCE):
        raise ValueError(f'Aggregate mismatch: {left} versus {right}')


def same_keys(actual, expected, label):
    if set(actual) != set(expected):
        raise ValueError(f'Incomplete or unexpected {label}')


def quantile(values, probability):
    """Linear interpolation at (n - 1) * probability, matching the saved envelopes."""
    ordered = sorted(number(v) for v in values)
    if not ordered or not 0 <= probability <= 1:
        raise ValueError('Quantile needs values and a probability in [0, 1]')
    position = (len(ordered) - 1) * probability
    lower, upper = math.floor(position), math.ceil(position)
    return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])


def expected_datasets():
    observed = {(cohort, 'observed', '-1', str(selection))
                for cohort in ('primary', 'broader_centers') for selection in range(3)}
    simulated = {('primary', generator, str(replicate), str(selection))
                 for generator in GENERATORS for replicate in range(99) for selection in range(3)}
    return observed | simulated


def reconstruct(data):
    """Validate complete saved evidence and recompute aggregates; never fit a model."""
    data = Path(data)
    cohort = json.loads((data / 'cohort_summary.json').read_text())
    cohort_fields = {'original_cells', 'original_blocks', 'supported_cells', 'supported_blocks',
        'primary_cells', 'primary_blocks', 'outside_mask_cells', 'supported_precision_exclusions',
        'available_native_cells', 'primary_available_native_cells'}
    same_keys(cohort, cohort_fields, 'cohort fields')
    if any(type(v) is not int or v < 0 for v in cohort.values()):
        raise ValueError('Invalid cohort counts')
    near(cohort['original_cells'], cohort['supported_cells'] + cohort['outside_mask_cells'])
    near(cohort['supported_cells'], cohort['primary_cells'] + cohort['supported_precision_exclusions'])
    if not (0 < cohort['primary_blocks'] <= cohort['supported_blocks'] <= cohort['original_blocks']
            and cohort['primary_cells'] >= cohort['primary_blocks']
            and cohort['available_native_cells'] >= cohort['primary_available_native_cells'] > 0):
        raise ValueError('Invalid cohort support ordering')

    fold_rows = read_csv(data / 'mixture_fold_scores.csv')
    folds = index(fold_rows, DATA_KEYS + ('fold',))
    expected = expected_datasets()
    same_keys(folds, {key + (str(fold),) for key in expected for fold in range(5)},
              'mixture fold slots (expected 4,485)')
    stats = []
    for key in sorted(expected):
        part = [folds[key + (str(fold),)] for fold in range(5)]
        for row in part:
            train = integer(row['train_blocks'], 20)
            held = integer(row['held_blocks'], 5)
            available_blocks = cohort['primary_blocks' if key[0] == 'primary' else 'supported_blocks']
            if train + held > available_blocks:
                raise ValueError('Training and held counts exceed available blocks')
            single, mixture = number(row['single_mean_log_density']), number(row['mixture_mean_log_density'])
            near(row['gain_mean'], mixture - single)
            near(row['gain_sum'], held * number(row['gain_mean']))
            mass = number(row['smallest_component_mass'])
            if not 0 <= mass <= train or boolean(row['small_component']) != (mass < 5):
                raise ValueError('Inconsistent component diagnostic')
            integer(row['single_iterations'], 1, 500)
            integer(row['mixture_iterations'], 1, 500)
            if key[1] != 'observed':
                reference = folds[('primary', 'observed', '-1', key[3], row['fold'])]
                near(train, reference['train_blocks'])
                near(held, reference['held_blocks'])
        held = sum(integer(row['held_blocks'], 5) for row in part)
        stats.append(dict(zip(DATA_KEYS, key),
                          statistic=math.fsum(number(row['gain_sum']) for row in part) / held,
                          held_blocks=held,
                          small_component_folds=sum(boolean(row['small_component']) for row in part)))
    reconstructed = index(stats, DATA_KEYS)
    saved = index(read_csv(data / 'organization_statistics.csv'), DATA_KEYS)
    same_keys(saved, expected, 'organization datasets (expected 897)')
    for key, row in reconstructed.items():
        near(row['statistic'], saved[key]['statistic'])
        near(row['held_blocks'], integer(saved[key]['held_blocks'], 1))
        near(row['small_component_folds'], integer(saved[key]['small_component_folds'], 0, 5))

    envelopes = []
    saved_envelopes = index(read_csv(data / 'envelope_summary.csv'), ('generator', 'selection'))
    same_keys(saved_envelopes, {(g, str(s)) for g in GENERATORS for s in range(3)}, 'nine envelopes')
    for generator in GENERATORS:
        for selection in range(3):
            observed = reconstructed[('primary', 'observed', '-1', str(selection))]['statistic']
            values = [reconstructed[('primary', generator, str(rep), str(selection))]['statistic']
                      for rep in range(99)]
            row = dict(generator=generator, selection=selection, observed=observed,
                       simulated_05=quantile(values, .05), simulated_median=quantile(values, .5),
                       simulated_95=quantile(values, .95),
                       simulations_at_or_above_observed=sum(value >= observed for value in values),
                       simulations=len(values))
            row['observed_above_95'] = observed > row['simulated_95']
            reference = saved_envelopes[(generator, str(selection))]
            for field in ('observed', 'simulated_05', 'simulated_median', 'simulated_95'):
                near(row[field], reference[field])
            for field in ('simulations', 'simulations_at_or_above_observed'):
                near(row[field], integer(reference[field], 0, 99))
            if row['observed_above_95'] != boolean(reference['observed_above_95']):
                raise ValueError('Changed upper-envelope flag')
            envelopes.append(row)

    association_rows = read_csv(data / 'association_fold_scores.csv')
    association = index(association_rows, ('model', 'fold'))
    same_keys(association, {(m, str(f)) for m in MODELS for f in range(5)}, 'association folds')
    association_summary = []
    reference_summary = index(read_csv(data / 'association_summary.csv'), ('model',))
    same_keys(reference_summary, {(m,) for m in MODELS}, 'association models')
    for fold in range(5):
        left, right = (association[(model, str(fold))] for model in MODELS)
        for field in ('held_blocks', 'held_cells'):
            near(integer(left[field], 1), integer(right[field], 1))
        if integer(left['held_cells'], 1) < integer(left['held_blocks'], 1):
            raise ValueError('Fewer cells than held blocks')
    for model in MODELS:
        part = [association[(model, str(fold))] for fold in range(5)]
        blocks = sum(integer(row['held_blocks'], 1) for row in part)
        cells = sum(integer(row['held_cells'], 1) for row in part)
        if blocks > cohort['primary_blocks'] or cells > cohort['primary_cells']:
            raise ValueError('Held support exceeds the primary cohort')
        row = dict(model=model,
                   mean_held_log_density=math.fsum(number(r['mean_log_density']) * integer(r['held_blocks'], 1)
                                                  for r in part) / blocks,
                   held_blocks=blocks, held_cells=cells)
        for field in ('mean_held_log_density', 'held_blocks', 'held_cells'):
            near(row[field], reference_summary[(model,)][field])
        association_summary.append(row)

    variance = read_csv(data / 'feature_variance.csv')
    same_keys(index(variance, ('feature',)), {(f,) for f in FEATURES}, 'six climate features')
    for row in variance:
        fraction = number(row['within_block_variance_fraction'])
        if not 0 <= fraction <= 1:
            raise ValueError('Variance fraction outside [0, 1]')
        near(100 * fraction, row['within_block_variance_percent'])
    return dict(statistics=stats, envelopes=envelopes, association_summary=association_summary,
                association_folds=association_rows, feature_variance=variance, cohort=cohort)


def write_csv(path, rows):
    with path.open('w', newline='', encoding='utf-8') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def observed_table(result):
    lookup = index(result['statistics'], DATA_KEYS)
    envelopes = index(result['envelopes'], ('generator', 'selection'))
    text = ('| Cohort | Selection | Held blocks | Two-minus-one gain (nats/block) | Shared at or above | Gradual at or above |\n'
            '|---|---|---:|---:|---:|---:|\n')
    for cohort, label in [('primary', 'Known uncertainty at most 5 km'), ('broader_centers', 'Broader centers')]:
        for selection, selection_label in enumerate(['Primary', 'Sensitivity 1', 'Sensitivity 2']):
            row = lookup[(cohort, 'observed', '-1', str(selection))]
            counts = [f"{envelopes[(g, str(selection))]['simulations_at_or_above_observed']}/99"
                      if cohort == 'primary' else 'Not simulated' for g in MODELS]
            text += f"| {label} | {selection_label} | {row['held_blocks']} | {row['statistic']:+.6f} | {' | '.join(counts)} |\n"
    return text + '\nSimulation counts are fitted-model diagnostics, not calibrated p-values.\n'


def association_table(result):
    lookup = index(result['association_folds'], ('model', 'fold'))
    text = ('| Stored fold | Held blocks | Shared versus area (nats) | Gradual versus area (nats) |\n'
            '|---|---:|---:|---:|\n')
    for fold in range(5):
        rows = [lookup[(model, str(fold))] for model in MODELS]
        text += f"| {fold} | {rows[0]['held_blocks']} | {float(rows[0]['mean_log_density']):+.6f} | {float(rows[1]['mean_log_density']):+.6f} |\n"
    rows = result['association_summary']
    text += f"| All held blocks | {rows[0]['held_blocks']} | {rows[0]['mean_held_log_density']:+.6f} | {rows[1]['mean_held_log_density']:+.6f} |\n"
    return text + '\nThe aggregate weights held blocks equally, not folds equally.\n'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=ROOT / 'data/environmental_organization')
    parser.add_argument('--output', type=Path, default=ROOT / 'environmental_replay_output')
    args = parser.parse_args()
    data, out = args.data.resolve(), args.output.resolve()
    if out.exists():
        raise FileExistsError('Choose a new output directory; existing work is never overwritten')
    result = reconstruct(data)
    out.mkdir(parents=True)
    for name in ['statistics', 'envelopes', 'association_summary']:
        write_csv(out / f'{name}.csv', result[name])
    (out / 'observed_table.md').write_text(observed_table(result), encoding='utf-8')
    (out / 'association_table.md').write_text(association_table(result), encoding='utf-8')
    receipt = dict(status='completed_environmental_organization_aggregate_replay',
        created_utc=datetime.now(timezone.utc).isoformat(), python=sys.version.split()[0],
        organization_fold_rows=4485, organization_statistics=len(result['statistics']),
        envelopes=len(result['envelopes']), association_fold_rows=len(result['association_folds']),
        observed_representative_datasets=6, simulated_representative_datasets=297,
        representative_role_layouts=3, tolerance=TOLERANCE,
        inputs={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(data.iterdir()) if p.suffix in ['.csv', '.json']},
        new_fits=0, simulations_rerun=False, native_predictions_reproduced=False,
        source_records_rebuilt=False, native_variance_recomputed=False,
        simulation_counts_are_calibrated_pvalues=False,
        scope='Saved aggregate arithmetic and internal consistency only; not independent biological validation.')
    (out / 'receipt.json').write_text(json.dumps(receipt, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: receipt[k] for k in ['status', 'organization_statistics', 'envelopes']}))


if __name__ == '__main__':
    main()
