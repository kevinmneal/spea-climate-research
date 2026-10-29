#!/usr/bin/env python3
"""One soft weighting comparison; all original presences and cutoffs retained."""
import os
for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "MKL_NUM_THREADS"]:
    os.environ.setdefault(key, "1")
import argparse
from dataclasses import replace
import json
from pathlib import Path
import resource
import shutil
import sys
import time
import traceback
import zipfile

import catboost
import numpy as np
import pandas as pd

from run_climate_statistics_comparison import Budget
from spea.climate_statistics import FAMILIES, StatisticsReader
from spea.spatial_thinning import inherited_cell_rows
from spea.native_objective_comparison import encode_objective, fit_native_objective, held_diagnostics
from spea.provenance import ROOT, sha256, utcnow, write_json
from spea.objective_measures import binary_class_weights


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/soft_spatial_weighting_v1.json")
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    config_path = ROOT / args.config
    cfg = json.loads(config_path.read_text())
    parent = ROOT / "runs" / cfg["parent_run"]
    parent_record = json.loads((parent / "run.json").read_text())
    if (parent_record["status"] != "completed" or catboost.__version__ != "1.2.10"
            or cfg["families"] != ["bioclim", "monthly"]
            or cfg["sampling_arms"] != {"soft_geo_10km": 10} or cfg["model"] != parent_record["config"]["model"]):
        raise ValueError("Declared parent, families, weighting arm or model changed")
    out = ROOT / "runs" / args.run_id
    out.mkdir(exist_ok=False)
    for folder in ["models", "borders", "predictions"]:
        (out / folder).mkdir()
    budget = Budget(cfg["budget"])
    start = time.monotonic()
    inputs, stats, metrics, fits, cellscores = {}, {}, [], [], []
    record = dict(status="preparing", created_utc=utcnow(), config=cfg, command_argv=sys.argv,
        parent_run=cfg["parent_run"], input_hashes=inputs, empirical_fits=0, full_data_fits=0,
        occurrence_changes=0, background_changes=0, scientific_gate_changes=0, historical_projections=0)

    def track(path, expected=None):
        rel = str(path.relative_to(ROOT))
        if rel not in inputs:
            inputs[rel] = sha256(path)
            stat = path.stat(); stats[rel] = (stat.st_size, stat.st_mtime_ns)
        if expected is not None and inputs[rel] != expected:
            raise ValueError("Changed frozen asset " + rel)
        return path

    def inherited(run, name):
        path = ROOT / "runs" / run / name
        return track(path, parent_record["input_hashes"][str(path.relative_to(ROOT))])

    def flush():
        for name, values in [("metrics.csv", metrics), ("fit_audit.csv", fits), ("cell_scores.csv.gz", cellscores)]:
            if values:
                pd.DataFrame(values).to_csv(out / name, index=False)
        record.update(cpu_seconds=time.process_time()-budget.start, elapsed_seconds=time.monotonic()-start,
            peak_resident_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        write_json(out / "run.json", record)

    try:
        track(config_path); track(parent / "run.json")
        sampling_dir = ROOT / "runs" / cfg["sampling_preparation_run"]
        srec = json.loads(track(sampling_dir / "run.json").read_text())
        if (srec["status"] != "completed_diagnostic_no_fits" or srec["unsupported_slots"] != 0
                or srec["input_hashes"][args.config] != inputs[args.config]):
            raise ValueError("Follow-up sampling was not frozen or lacks declared support")
        bio = json.loads(inherited(cfg["preparation_run"], "feature_registry.json").read_text())
        monthly = json.loads(inherited(cfg["monthly_run"], "registry.json").read_text())
        reader = StatisticsReader(ROOT, bio, monthly)
        for path in reader.paths.values():
            track(path, parent_record["input_hashes"][str(path.relative_to(ROOT))])
        genus = pd.read_csv(inherited(cfg["preparation_run"], "cells.csv.gz"))
        ns = pd.read_csv(inherited(cfg["north_south_run"], "species_cells.csv"))
        code = [Path(__file__), ROOT / "scripts/run_climate_statistics_comparison.py",
                ROOT / "src/spea/climate_statistics.py", ROOT / "src/spea/spatial_thinning.py", ROOT / "src/spea/spatial_weights.py",
                ROOT / "src/spea/native_objective_comparison.py", ROOT / "src/spea/objective_measures.py",
                ROOT / "src/spea/native_benchmark.py", ROOT / "catalog/review_gates.json", ROOT / "uv.lock"]
        for path in code:
            expected = parent_record["input_hashes"].get(str(path.relative_to(ROOT)))
            track(path, expected)
        write_json(out / "effective_config.json", cfg)
        with zipfile.ZipFile(out / "source_snapshot.zip", "w", zipfile.ZIP_DEFLATED) as archive:
            for path in code + [config_path]:
                archive.write(path, str(path.relative_to(ROOT)))
        record.update(status="fitting", design_frozen_utc=utcnow(), declared_model_slots=42,
                      supported_model_count=42, unsupported_sampling_slots=0)
        flush()
        for species_index, species in enumerate(cfg["species"]):
            is_ns = species in ["hammondii_north", "spea_south"]
            run = cfg["north_south_run"] if is_ns else cfg["preparation_run"]
            for fold in range(3 if is_ns else 5):
                budget.check()
                prefix = f"{species}_fold{fold}"
                split_path = inherited(run, f"splits/{prefix}.npz")
                with np.load(split_path) as z:
                    split = {k: z[k] for k in ["train_native", "q_native", "q_weights_m2", "test_native",
                                              "target_native", "target_area_m2", "test_cell_rows"]}
                weights_path = track(sampling_dir / "weights" / (prefix+".npz"), srec["output_hashes"]["weights/"+prefix+".npz"])
                with np.load(weights_path) as z:
                    positive_weights = z["presence_weights"].copy()
                    if not np.array_equal(z["eligible_native"], split["train_native"]):
                        raise ValueError("Weighted source membership changed")
                n = len(split["train_native"])
                if (positive_weights.shape != (n,) or not np.isfinite(positive_weights).all()
                        or np.any(positive_weights <= 0) or not np.isclose(positive_weights.sum(), n)):
                    raise ValueError("Every original presence needs positive normalized weight")
                mask = np.ones(n, dtype=bool)
                encoding = encode_objective("bce", split["train_native"], split["q_native"],
                    split["q_weights_m2"], cfg["model"]["leaf_damping_per_presence"])
                weights = binary_class_weights(positive_weights, split["q_weights_m2"], loss_mass=n)
                empirical = np.zeros(len(split["q_native"]))
                empirical[np.searchsorted(split["q_native"], split["train_native"])] = weights[:n]
                encoding = replace(encoding, weight=weights, empirical_on_q=empirical,
                                   native_l2_leaf_reg=encoding.effective_leaf_damping/weights.mean())
                training, _ = reader.read(encoding.native)
                target, _ = reader.read(split["target_native"])
                positions = np.searchsorted(split["target_native"], split["test_native"])
                held_rows = inherited_cell_rows(ns if is_ns else genus, split["test_cell_rows"])
                if not np.array_equal(held_rows.native_flat, split["test_native"]):
                    raise ValueError("Held source rows changed")
                blocks = held_rows.block50.to_numpy()
                seed = cfg["model_seed"] + species_index*100 + fold
                audit_idx = np.unique(np.linspace(0, len(target)-1, min(256, len(target)), dtype=int))
                for family in cfg["families"]:
                    budget.check()
                    if record["empirical_fits"] >= 42:
                        raise RuntimeError("Follow-up fit limit reached")
                    border_name = prefix+"_"+family+".tsv"
                    parent_border = track(parent / "borders" / border_name,
                        parent_record["output_hashes"]["borders/"+border_name])
                    borders = out / "borders" / border_name
                    shutil.copy2(parent_border, borders)
                    names = FAMILIES[family]
                    stem = prefix+"_soft_geo_10km_"+family
                    record["empirical_fits"] += 1
                    budget.reason = None
                    model, details = fit_native_objective(training[names], encoding, "bce", cfg["model"],
                        borders, seed, callbacks=[budget], threads=2)
                    if budget.reason:
                        raise RuntimeError(budget.reason)
                    model.save_model(str(out / "models" / (stem+".cbm")))
                    eta = np.empty(len(target))
                    for a in range(0, len(target), cfg["budget"]["prediction_chunk_size"]):
                        budget.check(); b = a+cfg["budget"]["prediction_chunk_size"]
                        eta[a:b] = model.predict(target.iloc[a:b][names], prediction_type="RawFormulaVal", thread_count=2)
                    novel = np.any((target.iloc[positions][names].to_numpy() < training[names].min().to_numpy()) |
                                   (target.iloc[positions][names].to_numpy() > training[names].max().to_numpy()), axis=1)
                    result, scores, captures = held_diagnostics(eta[positions], eta, split["target_area_m2"], blocks, novel)
                    result.update(species=species, fold=fold, sampling="soft_geo_10km", family=family,
                        features=len(names), train_cells=int(mask.sum()), held_cells=len(positions),
                        held_blocks=len(np.unique(blocks)), target_cells=len(target), stem=stem)
                    metrics.append(result)
                    for native, block, score in zip(split["test_native"], blocks, scores):
                        cellscores.append(dict(species=species, fold=fold, sampling="soft_geo_10km", family=family,
                            native_flat=int(native), block50=int(block), score=float(score)))
                    np.savez_compressed(out / "predictions" / (stem+".npz"), target_raw_eta=eta,
                        held_raw_eta=eta[positions], held_scores=scores, audit_indices=audit_idx, audit_raw_eta=eta[audit_idx])
                    details.update(species=species, fold=fold, sampling="soft_geo_10km", family=family, feature_names=names,
                        candidate=cfg["model"], seed=seed, borders="borders/"+border_name,
                        borders_sha256=sha256(borders), inherited_borders_sha256=sha256(parent_border),
                        source_split=str(split_path.relative_to(ROOT)), target_log_normalizer=result["target_log_normalizer"],
                        background_cells_area_inclusion_weights_unchanged=True,
                        all_original_presence_rows_retained=True, positive_weights_training_only=True,
                        frozen_weight_vector_sha256=sha256(weights_path),
                        normalized_positive_weight_min=float(weights[:n].min()),
                        normalized_positive_weight_max=float(weights[:n].max()))
                    write_json(out / "models" / (stem+"_metadata.json"), details)
                    fits.append(dict(species=species, fold=fold, sampling="soft_geo_10km", family=family,
                        stem=stem, status="completed", training_metric_error=details["training_metric_error"]))
                del target, training, eta
                flush()
                if sum(p.stat().st_size for p in out.rglob("*") if p.is_file()) > cfg["budget"]["output_gib"] * 1024**3:
                    raise RuntimeError("Follow-up output budget reached")
                print(json.dumps(dict(species=species, fold=fold, fits=record["empirical_fits"])), flush=True)
        for rel, expected in stats.items():
            stat = (ROOT / rel).stat()
            if (stat.st_size, stat.st_mtime_ns) != expected:
                raise ValueError("Input modified during follow-up")
        if len(metrics) != 42:
            raise ValueError("Incomplete follow-up")
        record.update(status="completed", finished_utc=utcnow())
        flush()
        record["output_hashes"] = {str(p.relative_to(out)): sha256(p) for p in out.rglob("*")
                                   if p.is_file() and p.name != "run.json"}
        write_json(out / "run.json", record)
        print(json.dumps(dict(status="completed", fits=42, cpu_seconds=record["cpu_seconds"])), flush=True)
    except Exception as error:
        record.update(status="failed_partial_results_preserved", error=str(error), traceback=traceback.format_exc())
        flush()
        raise


if __name__ == "__main__":
    main()
