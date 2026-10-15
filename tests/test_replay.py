"""Portable replay checks, including rejection of changed or missing evidence."""
import csv
import importlib.util
from pathlib import Path
import shutil
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("replay", ROOT / "scripts/replay.py")
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)


class ReplayTests(unittest.TestCase):
    def test_complete_saved_comparison_and_adverse_reference(self):
        families, folds, summaries = replay.matched(ROOT / "data/matched")
        self.assertEqual((len(families), len(folds), len(summaries)), (40, 168, 40))
        # Preserve the adverse result as well as the positive simple-control result.
        reference = [r for r in summaries if r["contrast"] == "monthly_minus_bioclim"]
        self.assertEqual(len(reference), 10)
        self.assertTrue(all(r["mean_block_log_score"] < 0 for r in reference))

    def test_corrupted_or_missing_fold_is_rejected(self):
        for change in ["missing", "duplicate", "changed_score"]:
            with self.subTest(change=change), tempfile.TemporaryDirectory() as tmp:
                dst = Path(tmp) / "matched"
                shutil.copytree(ROOT / "data/matched", dst)
                rows = replay.read(dst / "fold_metrics.csv")
                if change == "missing":
                    rows.pop()
                elif change == "duplicate":
                    rows.append(rows[0])
                else:
                    rows[0]["mean_block_log_score"] = str(float(rows[0]["mean_block_log_score"]) + .1)
                replay.write_csv(dst / "fold_metrics.csv", rows)
                with self.assertRaises(ValueError):
                    replay.matched(dst)


if __name__ == "__main__":
    unittest.main()
