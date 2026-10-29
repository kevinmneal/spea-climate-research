#!/usr/bin/env python3
"""Fixed climate-summary by spatial-declustering comparison; no tuning."""
import os
for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ.setdefault(key, "1")

import argparse
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

from spea.climate_statistics import FAMILIES, UNION, StatisticsReader, subset_borders
from spea.spatial_thinning import inherited_cell_rows
from spea.native_objective_comparison import encode_objective, save_common_borders, fit_native_objective, held_diagnostics
from spea.provenance import ROOT, sha256, utcnow, write_json


class Budget:
    def __init__(self, config):
        self.config = config
        self.start = time.process_time()
        self.reason = None

    def check(self):
        if time.process_time() - self.start > 60 * self.config["cpu_minutes"]:
            raise RuntimeError("Declared CPU budget reached")
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > self.config["peak_gib"] * 1024 ** 3:
            raise RuntimeError("Declared memory budget reached")

    def after_iteration(self, info):
        try:
            self.check()
            return True
        except RuntimeError as error:
            self.reason = str(error)
            return False


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/climate_statistics_comparison_v1.json")
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    if Path(args.run_id).name != args.run_id:
        raise ValueError("Run ID must be one component")
    config_path = ROOT / args.config
    cfg = json.loads(config_path.read_text())
    if cfg["families"] != list(FAMILIES) or cfg["family_dimensions"] != [len(v) for v in FAMILIES.values()]:
        raise ValueError("Frozen feature families changed")
    if catboost.__version__ != "1.2.10":
        raise ValueError("Require audited CatBoost1.2.10")
    out = ROOT / "runs" / args.run_id
    out.mkdir(exist_ok=False)
    for folder in ["models", "borders", "predictions"]:
        (out / folder).mkdir()
    start = time.monotonic()
    budget = Budget(cfg["budget"])
    hashes, input_stats = {}, {}
    metrics, fits, cellscores, diagnostics, unsupported = [], [], [], [], []
    record = dict(status="preparing", created_utc=utcnow(), command_argv=sys.argv,
                  config=cfg, input_hashes=hashes, empirical_fits=0, full_data_fits=0,
                  occurrence_changes=0, background_changes=0, scientific_gate_changes=0,
                  training_contribution_changes="declared training-only thinning arms",
                  historical_projections=0, genetic_reanalyses=0)

    def track(path, expected=None):
        rel = str(path.relative_to(ROOT))
        if rel not in hashes:
            value = sha256(path)
            if expected is not None and expected != value:
                raise ValueError("Changed input " + rel)
            hashes[rel] = value
            stat = path.stat()
            input_stats[rel] = (stat.st_size, stat.st_mtime_ns)
        elif expected is not None and hashes[rel] != expected:
            raise ValueError("Inconsistent source receipt " + rel)
        return path

    source_records = {}
    def asset(run, name):
        directory = ROOT / "runs" / run
        if run not in source_records:
            source_records[run] = json.loads(track(directory / "run.json").read_text())
            if source_records[run]["status"] != "completed":
                raise ValueError("Source run incomplete")
        return track(directory / name, source_records[run]["output_hashes"][name])

    def flush():
        for name, rows in [("metrics.csv", metrics), ("fit_audit.csv", fits),
                           ("cell_scores.csv.gz", cellscores), ("feature_diagnostics.json", diagnostics),
                           ("unsupported_slots.csv", unsupported)]:
            if rows:
                if name.endswith(".json"):
                    write_json(out / name, rows)
                else:
                    pd.DataFrame(rows).to_csv(out / name, index=False)
        record.update(cpu_seconds=time.process_time() - budget.start,
                      elapsed_seconds=time.monotonic() - start,
                      peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write_json(out / "run.json", record)

    try:
        track(config_path)
        write_json(out / "effective_config.json", cfg)
        write_json(out / "feature_registry.json", dict(source="CHELSA2.1 1981–2010",
            union_columns=UNION, families=FAMILIES, definitions=cfg["definitions"],
            conventions=cfg["numerics"], interpretation=cfg["predictor_semantics"]))
        bio = json.loads(asset(cfg["preparation_run"], "feature_registry.json").read_text())
        monthly = json.loads(asset(cfg["monthly_run"], "registry.json").read_text())
        reader = StatisticsReader(ROOT, bio, monthly)
        for path in reader.paths.values():
            rel = str(path.relative_to(ROOT))
            expected = source_records[cfg["monthly_run"]]["input_hashes"].get(rel)
            if expected is None:
                expected = source_records[cfg["preparation_run"]]["input_hashes"].get(rel)
            if expected is None:
                raise ValueError("No previous primary-array receipt " + rel)
            track(path, expected)
        genus_cells = pd.read_csv(asset(cfg["preparation_run"], "cells.csv.gz"))
        ns_cells = pd.read_csv(asset(cfg["north_south_run"], "species_cells.csv"))
        sampling_out = ROOT / "runs" / cfg["sampling_preparation_run"]
        sampling_record = json.loads(track(sampling_out / "run.json").read_text())
        if (sampling_record["status"] != "completed_diagnostic_no_fits"
                or sampling_record["input_hashes"][args.config] != hashes[args.config]
                or sampling_record["empirical_fits"] != 0):
            raise ValueError("Sampling design was not frozen under this config before fits")

        def sampling_asset(name):
            return track(sampling_out / name, sampling_record["output_hashes"][name])

        sampling_audit = pd.read_csv(sampling_asset("sampling_audit.csv"))
        sampling_asset("thinning_attempts.csv")
        sampling_asset("presence_feature_support.csv")
        splits = []
        for species in cfg["species"]:
            is_ns = species in ["hammondii_north", "spea_south"]
            run = cfg["north_south_run"] if is_ns else cfg["preparation_run"]
            for fold in range(3 if is_ns else 5):
                split_path = asset(run, f"splits/{species}_fold{fold}.npz")
                with np.load(split_path) as saved:
                    split = {k: saved[k] for k in ["train_native", "test_native", "train_cell_rows", "test_cell_rows",
                        "q_native", "q_weights_m2", "q_inclusion_probability", "target_native", "target_area_m2"]}
                if (not np.isin(split["train_native"], split["q_native"]).all()
                        or not np.isin(split["test_native"], split["target_native"]).all()
                        or np.isin(split["q_native"], split["target_native"]).any()
                        or np.any(np.diff(split["target_native"]) <= 0)):
                    raise ValueError("Inherited roles/integration support inconsistent")
                masks_path = sampling_asset(f"masks/{species}_fold{fold}.npz")
                with np.load(masks_path) as saved:
                    if (not np.array_equal(saved["eligible_native"], split["train_native"])
                            or not np.array_equal(saved["eligible_cell_rows"], split["train_cell_rows"])):
                        raise ValueError("Thinning eligibility changed")
                    masks = {arm: saved[arm+"_mask"].copy() for arm in cfg["sampling_arms"]}
                    supported = {arm: bool(saved[arm+"_supported"]) for arm in cfg["sampling_arms"]}
                for arm, mask in masks.items():
                    if mask.dtype != bool or mask.shape != split["train_native"].shape:
                        raise ValueError("Invalid frozen thinning mask")
                    row = sampling_audit[(sampling_audit.species == species) &
                                         (sampling_audit.fold == fold) & (sampling_audit.sampling == arm)]
                    if (len(row) != 1 or int(row.iloc[0].retained_train_cells) != mask.sum()
                            or bool(row.iloc[0].supported) != supported[arm]):
                        raise ValueError("Frozen audit and mask disagree")
                splits.append((species, fold, split_path, split, ns_cells if is_ns else genus_cells,
                               masks, supported))
        if len(splits) != 21:
            raise ValueError("Require all21 inherited folds")
        code = [Path(__file__), ROOT / "src/spea/climate_statistics.py", ROOT / "src/spea/seasonal_model.py",
                ROOT / "src/spea/native_objective_comparison.py", ROOT / "src/spea/objective_measures.py",
                ROOT / "src/spea/genus_model.py", ROOT / "src/spea/native_benchmark.py",
                ROOT / "tests/test_climate_statistics.py", ROOT / "catalog/review_gates.json",
                ROOT / "src/spea/spatial_thinning.py", ROOT / "tests/test_spatial_thinning.py",
                ROOT / "scripts/prepare_climate_statistics_sampling.py",
                ROOT / "uv.lock", ROOT / "pyproject.toml"]
        for path in code:
            track(path)
        with zipfile.ZipFile(out / "source_snapshot.zip", "w", zipfile.ZIP_DEFLATED) as archive:
            for path in code + [config_path]:
                archive.write(path, str(path.relative_to(ROOT)))
        expected_fits = sampling_record["declared_supported_fits"]
        record.update(status="fitting", design_frozen_utc=utcnow(), inherited_splits=21,
                      declared_model_slots=315, supported_model_count=expected_fits,
                      unsupported_sampling_slots=sampling_record["unsupported_slots"],
                      sampling_preparation_run=cfg["sampling_preparation_run"])
        flush()
        print(json.dumps(dict(stage="frozen", folds=21, fits=expected_fits,
            unsupported_sampling_slots=sampling_record["unsupported_slots"], raw_array_hashes_checked=43)), flush=True)
        for species, fold, split_path, split, source_cells, masks, supported in splits:
            budget.check()
            prefix = f"{species}_fold{fold}"
            seed = cfg["model_seed"] + 100 * cfg["species"].index(species) + fold
            original_encoding = encode_objective("bce", split["train_native"], split["q_native"],
                                                 split["q_weights_m2"], cfg["model"]["leaf_damping_per_presence"])
            original_training, train_audit = reader.read(original_encoding.native)
            target, target_audit = reader.read(split["target_native"])
            positions = np.searchsorted(split["target_native"], split["test_native"])
            held_rows = inherited_cell_rows(source_cells, split["test_cell_rows"])
            held_blocks = held_rows.block50.to_numpy()
            if not np.array_equal(held_rows.native_flat, split["test_native"]):
                raise ValueError("Held biological cell rows changed")
            bank = out / "borders" / (prefix + "_union.tsv")
            save_common_borders(original_training, bank, seed=seed)
            diagnostics.append(dict(species=species, fold=fold, split=str(split_path.relative_to(ROOT)),
                training=train_audit, held_landscape=target_audit,
                training_rows_include_both_sources=True, target_cells=len(target)))
            audit_idx = np.unique(np.linspace(0, len(target)-1, min(256, len(target)), dtype=int))
            borders_by_family = {}
            for family, names in FAMILIES.items():
                borders_by_family[family] = out / "borders" / (prefix + "_" + family + ".tsv")
                subset_borders(bank, borders_by_family[family], UNION, names)
            for arm in cfg["sampling_arms"]:
                if not supported[arm]:
                    for family in cfg["families"]:
                        unsupported.append(dict(species=species, fold=fold, sampling=arm, family=family,
                            retained_train_cells=int(masks[arm].sum()),
                            status="unsupported_insufficient_training_support"))
                    continue
                presence = split["train_native"][masks[arm]]
                encoding = encode_objective("bce", presence, split["q_native"],
                    split["q_weights_m2"], cfg["model"]["leaf_damping_per_presence"])
                frame_rows = np.r_[np.flatnonzero(masks[arm]),
                    np.arange(len(split["train_native"]), len(original_encoding.native))]
                if not np.array_equal(original_encoding.native[frame_rows], encoding.native):
                    raise ValueError("Thinned features and objective do not align")
                training = original_training.iloc[frame_rows].reset_index(drop=True)
                for family in cfg["families"]:
                    budget.check()
                    if record["empirical_fits"] >= cfg["budget"]["maximum_fits"]:
                        raise RuntimeError("Declared fit-count budget reached")
                    names, borders = FAMILIES[family], borders_by_family[family]
                    stem = prefix + "_" + arm + "_" + family
                    record["empirical_fits"] += 1
                    tick = time.monotonic()
                    budget.reason = None
                    model, details = fit_native_objective(training[names], encoding, "bce", cfg["model"],
                        borders, seed, callbacks=[budget], threads=cfg["budget"]["threads"])
                    if budget.reason:
                        raise RuntimeError(budget.reason)
                    model.save_model(str(out / "models" / (stem + ".cbm")))
                    eta = np.empty(len(target))
                    for a in range(0, len(target), cfg["budget"]["prediction_chunk_size"]):
                        budget.check()
                        b = a + cfg["budget"]["prediction_chunk_size"]
                        eta[a:b] = model.predict(target.iloc[a:b][names], prediction_type="RawFormulaVal",
                                                 thread_count=cfg["budget"]["threads"])
                    if not np.isfinite(eta).all():
                        raise ValueError("Nonfinite prediction")
                    held_frame = target.iloc[positions][names]
                    novel = np.any((held_frame.to_numpy() < training[names].min().to_numpy()) |
                                   (held_frame.to_numpy() > training[names].max().to_numpy()), axis=1)
                    result, scores, captures = held_diagnostics(eta[positions], eta,
                        split["target_area_m2"], held_blocks, novel)
                    result.update(species=species, fold=fold, sampling=arm, family=family, features=len(names),
                                  train_cells=len(presence), held_cells=len(positions),
                                  held_blocks=len(np.unique(held_blocks)), target_cells=len(target), stem=stem)
                    metrics.append(result)
                    for native, block, value in zip(split["test_native"], held_blocks, scores):
                        cellscores.append(dict(species=species, fold=fold, sampling=arm, family=family,
                            native_flat=int(native), block50=int(block), score=float(value)))
                    np.savez_compressed(out / "predictions" / (stem + ".npz"), target_raw_eta=eta,
                        held_raw_eta=eta[positions], held_scores=scores,
                        audit_indices=audit_idx, audit_raw_eta=eta[audit_idx])
                    details.update(species=species, fold=fold, sampling=arm, family=family, feature_names=names,
                        candidate=cfg["model"], seed=seed, source_split=str(split_path.relative_to(ROOT)),
                        thinning_mask_file=f"runs/{cfg['sampling_preparation_run']}/masks/{species}_fold{fold}.npz",
                        retained_training_cells=len(presence), original_training_cells=len(split["train_native"]),
                        borders=str(borders.relative_to(out)), borders_sha256=sha256(borders),
                        target_log_normalizer=result["target_log_normalizer"], held_area_captures=captures,
                        background_cells_area_inclusion_weights_unchanged=True)
                    write_json(out / "models" / (stem + "_metadata.json"), details)
                    fits.append(dict(species=species, fold=fold, sampling=arm, family=family, stem=stem,
                        status="completed", seconds=time.monotonic()-tick,
                        training_metric_error=details["training_metric_error"]))
                del training
            del target, original_training, eta
            flush()
            if sum(p.stat().st_size for p in out.rglob("*") if p.is_file()) > cfg["budget"]["output_gib"] * 1024**3:
                raise RuntimeError("Declared output budget reached")
            print(json.dumps(dict(stage="fold_complete", species=species, fold=fold,
                fits=record["empirical_fits"], cpu_seconds=round(record["cpu_seconds"], 2))), flush=True)
        for rel, expected_stat in input_stats.items():
            stat = (ROOT / rel).stat()
            if (stat.st_size, stat.st_mtime_ns) != expected_stat:
                raise ValueError("Input modified during run " + rel)
        if len(metrics) != expected_fits or len(fits) != expected_fits:
            raise ValueError("Incomplete comparison")
        record.update(status="completed", finished_utc=utcnow())
        flush()
        # Finalize tables before hashing them; gzip rewrites can change headers.
        record["output_hashes"] = {str(p.relative_to(out)): sha256(p) for p in out.rglob("*")
                                  if p.is_file() and p.name != "run.json"}
        write_json(out / "run.json", record)
        print(json.dumps(dict(status="completed", fits=expected_fits, cpu_seconds=record["cpu_seconds"],
                             peak_resident_bytes=record["peak_resident_bytes"])), flush=True)
    except Exception as error:
        record.update(status="failed_partial_results_preserved", error=str(error), traceback=traceback.format_exc())
        flush()
        raise


if __name__ == "__main__":
    main()
