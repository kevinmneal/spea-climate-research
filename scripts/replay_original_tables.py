#!/usr/bin/env python3
"""Rebuild seasonal manuscript Tables 2/3 using only five aggregate CSV files.

Python standard library only. This replays saved score arithmetic, not models,
native normalizers, occurrence preparation, folds, maps or biological evidence.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path

UNITS = {
    "intermontana": ("Intermontana", 5),
    "bombifrons": ("Bombifrons", 5),
    "multiplicata": ("Multiplicata", 5),
    "hammondii_north": ("Northern hammondii", 3),
    "spea_south": ("Southern species", 3),
}
FAMILIES = ["annual", "variability", "monthly"]
METRICS = ["mean_block_log_score", "mean_cell_log_score",
           "area_weighted_source_auc", "source_average_precision_prior_half"]
CONTRASTS = {"monthly_minus_annual": ("monthly", "annual"),
             "monthly_minus_variability": ("monthly", "variability"),
             "variability_minus_annual": ("variability", "annual")}
FILES = ["primary_folds.csv", "primary_summary.csv", "primary_contrasts.csv",
         "procedure_contrasts.csv", "support_by_fold.csv"]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def near(a, b):
    if not math.isfinite(float(a)) or not math.isclose(float(a), float(b), rel_tol=0, abs_tol=1e-10):
        raise ValueError(f"Saved arithmetic differs: {a} versus {b}")


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def markdown(headers, rows):
    return ("| " + " | ".join(headers) + " |\n| " + " | ".join(["---"] * len(headers))
            + " |\n" + "".join("| " + " | ".join(row) + " |\n" for row in rows))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, type=Path, help="Directory containing the five aggregate CSVs")
    ap.add_argument("--output", required=True, type=Path, help="Fresh output directory; existing output is refused")
    args = ap.parse_args()
    source, out = args.input.resolve(), args.output.resolve()
    if out.exists():
        raise FileExistsError("Use a fresh output directory")
    data = {name: read(source / name) for name in FILES}
    slots = {}
    for row in data["primary_folds.csv"]:
        key = (row["species"], row["variant"], row["family"], int(row["fold"]))
        if key in slots:
            raise ValueError("Duplicate family/fold slot")
        slots[key] = row
    expected = {(s, v, f, k) for s, (_, n) in UNITS.items()
                for v in ["fixed", "selected"] for f in FAMILIES for k in range(n)}
    if set(slots) != expected:
        raise ValueError("Expected the complete original five-unit, fixed/selected design")

    means, paired, changes = [], {}, []
    for species, (_, count) in UNITS.items():
        for variant in ["fixed", "selected"]:
            for family in FAMILIES:
                part = [slots[(species, variant, family, fold)] for fold in range(count)]
                row = dict(species=species, variant=variant, family=family, folds=count)
                row.update({metric: math.fsum(float(r[metric]) for r in part) / count for metric in METRICS})
                means.append(row)
            for contrast, (a, b) in CONTRASTS.items():
                values = []
                for fold in range(count):
                    left, right = slots[(species, variant, a, fold)], slots[(species, variant, b, fold)]
                    value = {metric: float(left[metric]) - float(right[metric]) for metric in METRICS[:2]}
                    paired[(species, variant, contrast, fold)] = value
                    values.append(value)
                changes.append(dict(species=species, procedure=variant, contrast=contrast, folds=count,
                    block_mean=math.fsum(r[METRICS[0]] for r in values) / count,
                    cell_mean=math.fsum(r[METRICS[1]] for r in values) / count,
                    positive_block_folds=sum(r[METRICS[0]] > 0 for r in values),
                    positive_cell_folds=sum(r[METRICS[1]] > 0 for r in values)))
    by_mean = {(r["species"], r["variant"], r["family"]): r for r in means}
    by_change = {(r["species"], r["procedure"], r["contrast"]): r for r in changes}
    if len(data["primary_summary.csv"]) != len(means) or len(data["primary_contrasts.csv"]) != len(paired):
        raise ValueError("Incomplete saved summary or contrast table")
    if {(r["species"], r["variant"], r["family"]) for r in data["primary_summary.csv"]} != set(by_mean):
        raise ValueError("Incomplete or duplicate family summaries")
    if {(r["species"], r["variant"], r["contrast"], int(r["fold"]))
            for r in data["primary_contrasts.csv"]} != set(paired):
        raise ValueError("Incomplete or duplicate paired contrasts")
    for row in data["primary_summary.csv"]:
        expected_row = by_mean[(row["species"], row["variant"], row["family"])]
        near(row["folds"], expected_row["folds"])
        for metric in METRICS:
            near(row[metric], expected_row[metric])
    for row in data["primary_contrasts.csv"]:
        expected_row = paired[(row["species"], row["variant"], row["contrast"], int(row["fold"]))]
        for metric in METRICS[:2]:
            near(row[metric], expected_row[metric])
    if len(data["procedure_contrasts.csv"]) != len(changes):
        raise ValueError("Incomplete saved procedure table")
    if {(r["species"], r["procedure"], r["contrast"]) for r in data["procedure_contrasts.csv"]} != set(by_change):
        raise ValueError("Incomplete or duplicate procedure contrasts")
    for row in data["procedure_contrasts.csv"]:
        expected_row = by_change[(row["species"], row["procedure"], row["contrast"])]
        for metric in ["block_mean", "cell_mean", "positive_block_folds", "positive_cell_folds", "folds"]:
            near(row[metric], expected_row[metric])

    support = data["support_by_fold.csv"]
    support_keys = [(r["species"], int(r["fold"])) for r in support]
    if len(set(support_keys)) != 21 or len(support) != 21 or sum(int(r["n_test"]) for r in support) != 3914:
        raise ValueError("Incomplete saved fold support")
    for row in support:
        r = slots[(row["species"], "selected", "annual", int(row["fold"]))]
        for key in ["n_train", "n_test", "n_test_blocks", "target_native_cells", "target_area_km2"]:
            near(row[key], r[key])
    table2, table3 = [], []
    for species, (label, count) in UNITS.items():
        fam = {f: by_mean[(species, "selected", f)] for f in FAMILIES}
        ma = by_change[(species, "selected", "monthly_minus_annual")]
        mv = by_change[(species, "selected", "monthly_minus_variability")]
        table2.append([label] + [f'{fam[f][METRICS[0]]:+.3f}' for f in FAMILIES] +
                      [f'{ma["block_mean"]:+.3f}', f'{mv["block_mean"]:+.3f}; {mv["positive_block_folds"]}/{count}',
                       f'{fam["monthly"][METRICS[2]]:.3f}', f'{fam["monthly"][METRICS[3]]:.3f}'])
        row = [label]
        for contrast, variant in [("monthly_minus_annual", "selected"), ("monthly_minus_annual", "fixed"),
                                  ("monthly_minus_variability", "selected"), ("monthly_minus_variability", "fixed")]:
            r = by_change[(species, variant, contrast)]
            row.append(f'{r["block_mean"]:+.3f} ({r["positive_block_folds"]}/{count})')
        table3.append(row)
    out.mkdir(parents=True)
    write_csv(out / "family_means.csv", means)
    write_csv(out / "procedure_contrasts.csv", changes)
    (out / "table2.md").write_text(markdown(
        ["Unit", "Annual", "Variability", "Monthly", "Monthly − annual",
         "Monthly − variability; positive folds", "Monthly ROC AUC", "Monthly AP, 1:1"], table2))
    (out / "table3.md").write_text(markdown(
        ["Unit", "Selected/fallback monthly − annual", "Fixed baseline monthly − annual",
         "Selected/fallback monthly − variability", "Fixed baseline monthly − variability"], table3))
    receipt = dict(status="completed_aggregate_table_replay", created_utc=datetime.now(timezone.utc).isoformat(),
                   inputs={name: sha(source / name) for name in FILES}, script_sha256=sha(Path(__file__)),
                   family_fold_slots=len(slots), family_mean_groups=len(means), paired_fold_contrasts=len(paired),
                   procedure_contrast_groups=len(changes), fold_support_rows=21, held_cells=3914,
                   tolerance=1e-10, full_native_normalizers_replayed=False, model_predictions_replayed=False,
                   source_records_rebuilt=False, new_fits=0, new_predictions=0,
                   outputs={p.name: sha(p) for p in out.iterdir()})
    (out / "replay_receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps({"status": receipt["status"], "family_fold_slots": len(slots), "tables": 2}))


if __name__ == "__main__":
    main()
