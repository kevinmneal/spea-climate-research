#!/usr/bin/env python3
"""Portable arithmetic replay of packaged, coordinate-free saved results only."""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean

PACKAGE = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--source-dir", type=Path, default=PACKAGE / "source_data")
parser.add_argument("--output", type=Path, default=PACKAGE / "rebuilt/aggregate_replay_receipt.json")
args = parser.parse_args()
SOURCE = args.source_dir.resolve()
checks = Counter()


def read(name):
    with (SOURCE / name).open(newline="") as f:
        return list(csv.DictReader(f))


def obj(name):
    return json.loads((SOURCE / name).read_text())


def close(actual, expected, category):
    assert math.isclose(float(actual), float(expected), rel_tol=1e-12, abs_tol=1e-12), (category, actual, expected)
    checks[category] += 1


def key(r, fields):
    return tuple(r[k] for k in fields)


manifest = obj("source_manifest.json")
for name, identity in manifest["packaged_files"].items():
    assert hashlib.sha256((SOURCE / name).read_bytes()).hexdigest() == identity["sha256"], name
    checks["packaged_file_sha256"] += 1

structured = read("structured_all_scores.csv")
settings = ["case", "geography", "ridge", "target_cohort"]
grouped = defaultdict(list)
index = {}
for r in structured:
    grouped[key(r, settings + ["structure"])].append(r)
    index[key(r, settings + ["fold", "group", "structure"])] = r
assert len(structured) == 384
for r in read("structured_all_means.csv"):
    for metric in ["cell_score", "block_score"]:
        close(mean(float(x[metric]) for x in grouped[key(r, settings + ["structure"])]), r[metric], "structured_all_means")
for r in read("structured_all_paired_gains.csv"):
    a, b = r["contrast"].split("-")
    prefix = key(r, settings + ["fold", "group"])
    close(float(index[prefix + (a,)][r["metric"]]) - float(index[prefix + (b,)][r["metric"]]), r["gain"], "structured_all_paired_gains")
primary = [r for r in structured if r["case"] == "all_1km" and r["geography"] == "spline" and float(r["ridge"]) == 1 and r["target_cohort"] == "all"]
for r in read("table1_controlled_responses.csv"):
    close(mean(float(x["cell_score"]) for x in primary if x["structure"] == r["structure"]), r["score_nats_relative_recording_density"], "figure1_saved_means")

season = read("seasonal_all_fold_metrics.csv")
metrics = ["mean_block_log_score", "mean_cell_log_score", "area_weighted_source_auc", "source_average_precision_prior_half"]
season_index = {key(r, ["species", "sampling", "family", "fold"]): r for r in season}
assert len(season_index) == 168
season_groups = defaultdict(list)
for r in season:
    season_groups[key(r, ["species", "sampling", "family"])].append(r)
for r in read("seasonal_all_family_means.csv"):
    rr = season_groups[key(r, ["species", "sampling", "family"])]
    for metric in metrics:
        close(mean(float(x[metric]) for x in rr), r[metric], "seasonal_family_means")
contrast_groups = defaultdict(list)
for r in read("seasonal_all_fold_contrasts.csv"):
    a, b = r["contrast"].split("_minus_")
    for metric in metrics:
        aa = season_index[key(r, ["species", "sampling"]) + (a, r["fold"])]
        bb = season_index[key(r, ["species", "sampling"]) + (b, r["fold"])]
        close(float(aa[metric]) - float(bb[metric]), r[metric], "seasonal_fold_contrasts")
    contrast_groups[key(r, ["species", "sampling", "contrast"])].append(r)
for r in read("seasonal_all_contrast_summaries.csv"):
    rr = contrast_groups[key(r, ["species", "sampling", "contrast"])]
    close(len(rr), r["folds"], "seasonal_contrast_fold_counts")
    for metric in metrics:
        dd = [float(x[metric]) for x in rr]
        for suffix, value in [("", mean(dd)), ("_min", min(dd)), ("_max", max(dd)), ("_positive_folds", sum(x > 0 for x in dd))]:
            close(value, r[metric + suffix], "seasonal_contrast_summary_statistics")
for r in read("seasonal_source_alignment.csv"):
    rr = [x for x in season if key(x, ["species", "fold", "sampling"]) == key(r, ["species", "fold", "sampling"])]
    assert len(rr) == 4
    for metric in ["train_cells", "held_cells", "held_blocks", "target_cells"]:
        for x in rr:
            close(x[metric], r[metric], "seasonal_saved_alignment_counts")

pool = read("pooling_transfer_scores.csv")
pool_fields = ["case", "family", "target_cohort", "train_mode", "group"]
pool_groups = defaultdict(list)
pool_index = {}
for r in pool:
    pool_groups[key(r, pool_fields)].append(r)
    pool_index[key(r, ["case", "fold", "family", "target_cohort", "group", "train_mode"])] = r
for r in read("pooling_transfer_means.csv"):
    for metric in ["cell_score", "block_score"]:
        close(mean(float(x[metric]) for x in pool_groups[key(r, pool_fields)]), r[metric], "pooling_transfer_means")
for r in read("pooling_transfer_paired_gains.csv"):
    prefix = key(r, ["case", "fold", "family", "target_cohort", "group"])
    a, b = pool_index[prefix + (r["within_mode"],)], pool_index[prefix + (r["comparator_mode"],)]
    for metric in ["cell", "block"]:
        close(float(a[metric + "_score"]) - float(b[metric + "_score"]), r[metric + "_gain"], "pooling_transfer_paired_gains")

water_metrics = read("occurrence_water_block_metrics.csv")
water_gains = read("occurrence_water_paired_gains.csv")
for r in water_gains:
    mm = {x["water"]: x for x in water_metrics if x["split"] == r["split"] and x["family"] == r["family"]}
    source_metric = "mean_source_block_score" if r["kind"] == "spatial" else "mean_record_score"
    for field, arm in [("climate_score", "False"), ("water_score", "True")]:
        close(mm[arm][source_metric], r[field], "occurrence_water_source_scores")
    close(float(r["water_score"]) - float(r["climate_score"]), r["water_gain"], "occurrence_water_paired_gains")
for family, summary in obj("occurrence_water_summary.json")["primary_equal_five_source_block_means"].items():
    rr = [r for r in water_gains if r["kind"] == "spatial" and r["family"] == family]
    for metric in ["climate_score", "water_score", "water_gain"]:
        close(mean(float(x[metric]) for x in rr), summary[metric], "occurrence_water_primary_means")
    close(sum(float(x["water_gain"]) > 0 for x in rr), summary["positive_water_gains"], "occurrence_water_positive_blocks")

survey_metrics = read("survey_water_block_metrics.csv")
survey_gains = read("survey_water_paired_gains.csv")
for r in survey_gains:
    mm = {x["arm"]: x for x in survey_metrics if x["held_block"] == r["held_block"]}
    for field, metric in [("joint_gain", "mean_joint_log_score"), ("any_encounter_gain", "mean_any_log_score"), ("conditional_pattern_gain", "mean_conditional_pattern_log_score")]:
        close(float(mm["terrestrial_water"][metric]) - float(mm["terrestrial"][metric]), r[field], "survey_water_paired_gains")
survey_summary = obj("survey_water_summary.json")
for field, saved in [("joint_gain", "primary_equal_block_joint_gain"), ("any_encounter_gain", "equal_block_any_encounter_gain"), ("conditional_pattern_gain", "equal_block_conditional_pattern_gain")]:
    close(mean(float(r[field]) for r in survey_gains), survey_summary[saved], "survey_water_equal_block_means")
close(sum(int(r["cells"]) * float(r["joint_gain"]) for r in survey_gains) / sum(int(r["cells"]) for r in survey_gains), survey_summary["secondary_pooled_cell_joint_gain"], "survey_water_cell_weighted_mean")
close(mean(float(r["joint_gain"]) for r in survey_gains if r["held_block"] != "5"), survey_summary["descriptive_equal_block_joint_gain_without_block5"], "survey_water_sparse_block_sensitivity")

original = read("seasonal_original_procedure_folds.csv")
original_index = {key(r, ["species", "variant", "family", "fold"]): r for r in original}
for r in read("seasonal_original_procedure_means.csv"):
    rr = [x for x in original if key(x, ["species", "variant", "family"]) == key(r, ["species", "variant", "family"])]
    for metric in ["mean_block_log_score", "mean_cell_log_score"]:
        close(mean(float(x[metric]) for x in rr), r[metric], "original_seasonal_procedure_means")
for r in read("seasonal_original_procedure_contrasts.csv"):
    a, b = r["contrast"].split("_minus_")
    dd = []
    for fold in range(int(r["folds"])):
        prefix = (r["species"], r["procedure"])
        aa, bb = original_index[prefix + (a, str(fold))], original_index[prefix + (b, str(fold))]
        dd.append({m: float(aa[m]) - float(bb[m]) for m in ["mean_block_log_score", "mean_cell_log_score"]})
    for metric, field in [("mean_block_log_score", "block"), ("mean_cell_log_score", "cell")]:
        close(mean(d[metric] for d in dd), r[field + "_mean"], "original_seasonal_procedure_contrasts")
        close(sum(d[metric] > 0 for d in dd), r["positive_" + field + "_folds"], "original_seasonal_procedure_positive_counts")

support = read("seasonal_support_by_fold.csv")
assert len(support) == 21
close(sum(int(r["n_test"]) for r in support), 3914, "saved_coverage_display_total")
coverage = obj("seasonal_coverage_map_provenance.json")
assert len(coverage["maps"]) == 5
for r in coverage["maps"]:
    close(sum(int(x["n_test"]) for x in support if x["species"] == r["species"]), r["displayed_held_cells"], "saved_coverage_per_unit_counts")

receipt = {"status": "passed", "created_utc": datetime.now(timezone.utc).isoformat(),
           "check_count": sum(checks.values()), "checks_by_category": dict(checks),
           "tolerance_absolute_and_relative": 1e-12,
           "source_manifest_sha256": hashlib.sha256((SOURCE / "source_manifest.json").read_bytes()).hexdigest(),
           "replay_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
           "scope": "Independent standard-library arithmetic from packaged saved aggregate files; source identity checks establish copy fidelity, not biological validity.",
           "score_units": {"structured_and_seasonal_and_pooling_and_occurrence_water": "relative recording-density nats under each design's declared normalization and evaluation weights; never pooled across designs", "survey_water": "nats per complete surveyed cell history; no common effect with density scores"},
           "new_fits": 0, "new_predictions": 0, "native_arrays_or_models_loaded": False,
           "native_audit": "Earlier 2510-check audit receipt is packaged separately as verification/aggregate_native_source_audit_receipt.json. No native audit repeated.",
           "map_rebuild": "Not included: S1 requires original saved native support and coarse point-grouping archive; this replay checks only packaged count identities."}
args.output.parent.mkdir(parents=True, exist_ok=True)
args.output.write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps({"status": receipt["status"], "check_count": receipt["check_count"], "receipt": str(args.output)}, indent=2))
