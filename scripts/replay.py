#!/usr/bin/env python3
"""Recompute published aggregate scores; no fitted models or localities needed."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
UNITS = {"intermontana": ("Intermontana", 5), "bombifrons": ("Bombifrons", 5),
         "multiplicata": ("Multiplicata", 5), "hammondii_north": ("Northern hammondii", 3),
         "spea_south": ("Southern hammondii", 3)}
ARMS = ("equal_cells", "soft_geo_10km")
FAMILIES = ("annual", "variability", "monthly", "bioclim")
METRICS = ("mean_block_log_score", "mean_cell_log_score", "area_weighted_source_auc",
           "source_average_precision_prior_half")
CONTRASTS = {"monthly_minus_annual": ("monthly", "annual"),
             "monthly_minus_variability": ("monthly", "variability"),
             "monthly_minus_bioclim": ("monthly", "bioclim"),
             "bioclim_minus_variability": ("bioclim", "variability")}


def read(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def index(rows, keys):
    result = {tuple(str(row[k]) for k in keys): row for row in rows}
    if len(result) != len(rows):
        raise ValueError(f"Duplicate keys: {keys}")
    return result


def near(left, right):
    if not math.isfinite(float(left)) or not math.isfinite(float(right)) or not math.isclose(
            float(left), float(right), rel_tol=0, abs_tol=1e-10):
        raise ValueError(f"Aggregate mismatch: {left} versus {right}")


def compare(actual, saved, keys):
    a, b = index(actual, keys), index(saved, keys)
    if set(a) != set(b):
        raise ValueError(f"Incomplete summary keys: {keys}")
    for key, row in a.items():
        for field, value in row.items():
            if field not in keys:
                near(value, b[key][field])


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def mean(values):
    return math.fsum(values) / len(values)


def matched(data):
    rows = read(data / "fold_metrics.csv")
    keys = ("species", "sampling", "family", "fold")
    slots = index(rows, keys)
    expected = {(species, arm, family, str(fold)) for species, (_, n) in UNITS.items()
                for arm in ARMS for family in FAMILIES for fold in range(n)}
    if set(slots) != expected:
        raise ValueError("Expected all 168 matched family/arm/fold slots")
    for species, (_, n) in UNITS.items():
        for fold in range(n):
            reference = slots[(species, ARMS[0], FAMILIES[0], str(fold))]
            for arm in ARMS:
                for family in FAMILIES:
                    row = slots[(species, arm, family, str(fold))]
                    for count in ["train_cells", "held_cells", "held_blocks", "target_cells"]:
                        near(row[count], reference[count])
                    for field in METRICS:
                        if not math.isfinite(float(row[field])):
                            raise ValueError("Nonfinite fold score")

    families, folds, summaries = [], [], []
    for species, (_, n) in UNITS.items():
        for arm in ARMS:
            for family in FAMILIES:
                row = dict(species=species, sampling=arm, family=family)
                row.update({m: mean([float(slots[(species, arm, family, str(k))][m])
                                     for k in range(n)]) for m in METRICS})
                families.append(row)
            for contrast, (left, right) in CONTRASTS.items():
                part = []
                for k in range(n):
                    row = dict(species=species, fold=k, sampling=arm, contrast=contrast)
                    row.update({m: float(slots[(species, arm, left, str(k))][m]) -
                                   float(slots[(species, arm, right, str(k))][m]) for m in METRICS})
                    folds.append(row)
                    part.append(row)
                row = dict(species=species, sampling=arm, contrast=contrast, folds=n)
                for m in METRICS:
                    values = [r[m] for r in part]
                    row.update({m: mean(values), m + "_min": min(values), m + "_max": max(values),
                                m + "_positive_folds": sum(v > 0 for v in values)})
                summaries.append(row)
    compare(families, read(data / "family_summary.csv"), ("species", "sampling", "family"))
    compare(folds, read(data / "contrast_folds.csv"), ("species", "sampling", "contrast", "fold"))
    compare(summaries, read(data / "contrast_summary.csv"), ("species", "sampling", "contrast"))
    return families, folds, summaries


def table(summaries):
    lookup = index(summaries, ("species", "sampling", "contrast"))
    text = ("| Unit | Training weights | Monthly − variability | Monthly − full BIOCLIM |\n"
            "|---|---|---:|---:|\n")
    for species, (label, n) in UNITS.items():
        for arm in ARMS:
            cells = []
            for contrast in ["monthly_minus_variability", "monthly_minus_bioclim"]:
                row = lookup[(species, arm, contrast)]
                cells.append(f'{row[METRICS[0]]:+.4f} ({row[METRICS[0]+"_positive_folds"]}/{n})')
            arm_label = "Equal" if arm == ARMS[0] else "Soft 10 km"
            text += f"| {label} | {arm_label} | {' | '.join(cells)} |\n"
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT / "data")
    parser.add_argument("--output", type=Path, default=ROOT / "replay_output")
    args = parser.parse_args()
    data, out = args.data.resolve(), args.output.resolve()
    if out.exists():
        raise FileExistsError("Choose a new output directory; existing work is never overwritten")
    families, folds, summaries = matched(data / "matched")
    out.mkdir(parents=True)
    subprocess.run([sys.executable, str(ROOT / "scripts/replay_original_tables.py"),
                    "--input", str(data / "original"), "--output", str(out / "original")], check=True)
    subprocess.run([sys.executable, str(ROOT / "scripts/replay_environmental_organization.py"),
                    "--data", str(data / "environmental_organization"),
                    "--output", str(out / "environmental_organization")], check=True)
    write_csv(out / "matched_family_means.csv", families)
    write_csv(out / "matched_fold_contrasts.csv", folds)
    write_csv(out / "matched_contrast_summary.csv", summaries)
    (out / "matched_primary_table.md").write_text(table(summaries), encoding="utf-8")
    hashes = {str(p.relative_to(data)): hashlib.sha256(p.read_bytes()).hexdigest()
              for directory in ["original", "matched", "environmental_organization"]
              for p in sorted((data / directory).iterdir()) if p.suffix in (".csv", ".json")}
    receipt = dict(status="completed_aggregate_replay", created_utc=datetime.now(timezone.utc).isoformat(),
                   python=sys.version.split()[0], original_slots=126, matched_slots=168,
                   matched_family_means=len(families), matched_fold_contrasts=len(folds),
                   matched_contrast_summaries=len(summaries), inputs=hashes, tolerance=1e-10,
                   environmental_organization_statistics=897, environmental_envelopes=9,
                   new_fits=0, native_predictions_reproduced=False, source_records_rebuilt=False,
                   historical_or_genetic_claims_tested=False)
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: receipt[k] for k in ["status", "original_slots", "matched_slots"]}))


if __name__ == "__main__":
    main()
