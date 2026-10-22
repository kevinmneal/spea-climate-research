"""Checks of the released arithmetic and deliberate corruption of its evidence."""
import importlib.util
from pathlib import Path
import shutil
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('environmental_replay', ROOT / 'scripts/replay_environmental_organization.py')
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)
DATA = ROOT / 'data/environmental_organization'


class EnvironmentalReplayTests(unittest.TestCase):
    def mutate(self, filename, change):
        with tempfile.TemporaryDirectory() as tmp:
            dst = Path(tmp) / 'data'
            shutil.copytree(DATA, dst)
            rows = replay.read_csv(dst / filename)
            change(rows)
            replay.write_csv(dst / filename, rows)
            with self.assertRaises(ValueError):
                replay.reconstruct(dst)

    def test_complete_scores_preserve_negative_result_and_sensitivity_mismatch(self):
        result = replay.reconstruct(DATA)
        self.assertEqual((len(result['statistics']), len(result['envelopes'])), (897, 9))
        observed = [r for r in result['statistics'] if r['cohort'] == 'primary' and r['generator'] == 'observed']
        for actual, expected in zip(observed, [-.14328483221264962, -.19165086590595448, -.008718556411569916]):
            self.assertAlmostEqual(actual['statistic'], expected, places=12)
        tail = [r for r in result['envelopes'] if r['selection'] == 2 and r['generator'] in ('shared', 'smooth')]
        self.assertEqual(len(tail), 2)
        self.assertTrue(all(r['observed_above_95'] and r['simulations_at_or_above_observed'] == 4 for r in tail))

    def test_linear_quantile_known_values(self):
        self.assertAlmostEqual(replay.quantile([30, 0, 20, 10], .05), 1.5, places=12)
        self.assertEqual(replay.quantile([30, 0, 20, 10], .5), 15)
        self.assertAlmostEqual(replay.quantile([30, 0, 20, 10], .95), 28.5, places=12)

    def test_unequal_fold_sizes_use_block_weighted_aggregate(self):
        result = replay.reconstruct(DATA)
        expected = {'shared': .17913054579848128, 'smooth': .11375851749893262}
        for row in result['association_summary']:
            self.assertEqual((row['held_blocks'], row['held_cells']), (134, 358))
            self.assertAlmostEqual(row['mean_held_log_density'], expected[row['model']], places=12)
            unweighted = sum(float(r['mean_log_density']) for r in result['association_folds'] if r['model'] == row['model']) / 5
            self.assertGreater(abs(row['mean_held_log_density'] - unweighted), .001)

    def test_missing_or_duplicate_mixture_fold_is_rejected(self):
        for change in [lambda rows: rows.pop(), lambda rows: rows.append(rows[0].copy())]:
            with self.subTest(change=change):
                self.mutate('mixture_fold_scores.csv', change)

    def test_nonfinite_or_inconsistent_fold_score_is_rejected(self):
        for field, value in [('gain_sum', 'nan'), ('single_mean_log_density', 'inf'), ('gain_mean', '99')]:
            with self.subTest(field=field):
                self.mutate('mixture_fold_scores.csv', lambda rows: rows[0].update({field: value}))

    def test_missing_dataset_and_changed_statistic_are_rejected(self):
        for change in [lambda rows: rows.pop(), lambda rows: rows[0].update(statistic='123')]:
            with self.subTest(change=change):
                self.mutate('organization_statistics.csv', change)

    def test_envelope_count_and_tail_flag_are_checked(self):
        for field, value in [('simulations_at_or_above_observed', '0'), ('simulations', '98'), ('observed_above_95', 'True')]:
            with self.subTest(field=field):
                self.mutate('envelope_summary.csv', lambda rows: rows[0].update({field: value}))

    def test_association_support_score_and_variance_are_checked(self):
        for filename, field, value in [('association_fold_scores.csv', 'held_blocks', '1'),
                                       ('association_fold_scores.csv', 'mean_log_density', 'nan'),
                                       ('association_summary.csv', 'mean_held_log_density', '99'),
                                       ('feature_variance.csv', 'within_block_variance_percent', '100')]:
            with self.subTest(filename=filename, field=field):
                self.mutate(filename, lambda rows: rows[0].update({field: value}))


if __name__ == '__main__':
    unittest.main()
