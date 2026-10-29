#!/usr/bin/env python3
"""Freeze training-only spatial thinning masks before the factorial fits."""
import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd

from spea.climate_statistics import StatisticsReader, UNION
from spea.spatial_thinning import geodesic_distances_km, thin_attempts, minimum_retained_distance, inherited_cell_rows
from spea.provenance import ROOT, sha256, utcnow, write_json


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", default="config/climate_statistics_comparison_v1.json")
    args = parser.parse_args()
    config_path = ROOT / args.config
    cfg = json.loads(config_path.read_text())
    out = ROOT / "runs" / cfg["sampling_preparation_run"]
    out.mkdir(exist_ok=False)
    (out / "masks").mkdir()
    inputs, records = {args.config: sha256(config_path)}, {}
    audit, attempts, support = [], [], []
    started = time.monotonic()

    def asset(run, name):
        directory = ROOT / "runs" / run
        if run not in records:
            records[run] = json.loads((directory / "run.json").read_text())
            inputs[str((directory / "run.json").relative_to(ROOT))] = sha256(directory / "run.json")
            if records[run]["status"] != "completed":
                raise ValueError("Incomplete source run")
        path = directory / name
        value = sha256(path)
        if value != records[run]["output_hashes"][name]:
            raise ValueError("Changed source asset " + str(path))
        inputs[str(path.relative_to(ROOT))] = value
        return path

    bio = json.loads(asset(cfg["preparation_run"], "feature_registry.json").read_text())
    monthly = json.loads(asset(cfg["monthly_run"], "registry.json").read_text())
    reader = StatisticsReader(ROOT, bio, monthly)
    # Primary arrays are checked once against existing receipts by the fitter.
    # Retain receipt links here; preparation neither changes nor admits cells.
    source_array_receipts = {}
    for path in reader.paths.values():
        rel = str(path.relative_to(ROOT))
        value = records[cfg["monthly_run"]]["input_hashes"].get(rel)
        if value is None:
            value = records[cfg["preparation_run"]]["input_hashes"].get(rel)
        if value is None:
            raise ValueError("Missing array receipt")
        source_array_receipts[rel] = value
    genus = pd.read_csv(asset(cfg["preparation_run"], "cells.csv.gz"))
    ns = pd.read_csv(asset(cfg["north_south_run"], "species_cells.csv"))
    for species_index, species in enumerate(cfg["species"]):
        is_ns = species in ["hammondii_north", "spea_south"]
        run = cfg["north_south_run"] if is_ns else cfg["preparation_run"]
        source_cells = ns if is_ns else genus
        for fold in range(3 if is_ns else 5):
            path = asset(run, f"splits/{species}_fold{fold}.npz")
            with np.load(path) as split:
                native = split["train_native"].copy()
                rows = split["train_cell_rows"].copy()
                test_count = len(split["test_native"])
            cells = inherited_cell_rows(source_cells, rows)
            if not np.array_equal(cells.native_flat, native):
                raise ValueError("Training source rows do not match native IDs")
            d = geodesic_distances_km(cells.longitude, cells.latitude)
            features, _ = reader.read(native)
            saved = dict(eligible_native=native, eligible_cell_rows=rows)
            seeds = [cfg["thinning"]["seed"] + species_index * 1000 + fold * 20 + i
                     for i in range(cfg["thinning"]["attempts"])]
            for arm, radius in cfg["sampling_arms"].items():
                if radius == 0:
                    masks, chosen = np.ones((1, len(native)), bool), 0
                    arm_seeds = [None]
                else:
                    masks, chosen = thin_attempts(d, radius, seeds)
                    arm_seeds = seeds
                mask = masks[chosen]
                retained = cells.iloc[np.flatnonzero(mask)]
                blocks = retained.block50.nunique()
                supported = (mask.sum() >= cfg["thinning"]["minimum_train_cells"]
                             and blocks >= cfg["thinning"]["minimum_train_blocks50"])
                minimum = minimum_retained_distance(d, mask)
                if radius and minimum is not None and minimum < radius:
                    raise ValueError("Selected thinning violates declared spacing")
                saved[arm + "_candidates"] = masks
                saved[arm + "_mask"] = mask
                saved[arm + "_chosen_attempt"] = np.array(chosen)
                saved[arm + "_supported"] = np.array(supported)
                for index, candidate in enumerate(masks):
                    intersection = np.sum(candidate & mask)
                    union = np.sum(candidate | mask)
                    attempts.append(dict(species=species, fold=fold, sampling=arm, attempt=index,
                        seed=arm_seeds[index], selected=index == chosen,
                        retained_cells=int(candidate.sum()),
                        retained_blocks50=int(cells.iloc[np.flatnonzero(candidate)].block50.nunique()),
                        minimum_spacing_km=minimum_retained_distance(d, candidate),
                        jaccard_with_selected=float(intersection / union)))
                country_counts, source_counts = {}, {}
                for field, counts in [("countries", country_counts), ("source_ids", source_counts)]:
                    for values in retained[field]:
                        for value in json.loads(values):
                            counts[value] = counts.get(value, 0) + 1
                audit.append(dict(species=species, fold=fold, sampling=arm, minimum_radius_km=radius,
                    eligible_train_cells=len(native), retained_train_cells=int(mask.sum()),
                    retained_fraction=float(mask.mean()), eligible_train_blocks50=int(cells.block50.nunique()),
                    retained_train_blocks50=int(blocks), held_cells=test_count, supported=bool(supported),
                    status="supported" if supported else "unsupported_insufficient_training_support",
                    chosen_attempt=chosen, chosen_seed=arm_seeds[chosen],
                    minimum_spacing_km=minimum,
                    contributor_records=int(retained.contributors.sum()),
                    unknown_radius_contributor_records=int(retained.unknown_radius_contributors.sum()),
                    fraction_cells_with_unknown_radius=float((retained.unknown_radius_contributors > 0).mean()),
                    median_guard_radius_unknown1km_m=float(retained.radius_unknown1km_m.median()),
                    maximum_guard_radius_unknown1km_m=float(retained.radius_unknown1km_m.max()),
                    country_cell_counts=json.dumps(country_counts, sort_keys=True),
                    source_cell_counts=json.dumps(source_counts, sort_keys=True),
                    source_split=str(path.relative_to(ROOT)), source_split_sha256=inputs[str(path.relative_to(ROOT))]))
                for feature in UNION:
                    original, kept = features[feature].to_numpy(), features.loc[mask, feature].to_numpy()
                    support.append(dict(species=species, fold=fold, sampling=arm, feature=feature,
                        eligible_minimum=float(original.min()), eligible_maximum=float(original.max()),
                        retained_minimum=float(kept.min()), retained_maximum=float(kept.max()),
                        retained_mean=float(kept.mean()), retained_sd=float(kept.std())))
            np.savez_compressed(out / "masks" / f"{species}_fold{fold}.npz", **saved)
            print(json.dumps(dict(species=species, fold=fold,
                retained=[int(saved[arm+"_mask"].sum()) for arm in cfg["sampling_arms"]])), flush=True)
    for name, values in [("sampling_audit.csv", audit), ("thinning_attempts.csv", attempts),
                         ("presence_feature_support.csv", support)]:
        pd.DataFrame(values).to_csv(out / name, index=False)
    for name in ["scripts/prepare_climate_statistics_sampling.py", "src/spea/spatial_thinning.py",
                 "src/spea/climate_statistics.py", "tests/test_spatial_thinning.py"]:
        inputs[name] = sha256(ROOT / name)
    write_json(out / "effective_config.json", cfg)
    write_json(out / "run.json", dict(status="completed_diagnostic_no_fits", finished_utc=utcnow(),
        config=cfg, input_hashes=inputs, primary_array_receipts=source_array_receipts,
        elapsed_seconds=time.monotonic()-started,
        sampling_slots=len(audit), supported_slots=int(sum(r["supported"] for r in audit)),
        unsupported_slots=int(sum(not r["supported"] for r in audit)),
        declared_supported_fits=int(sum(r["supported"] for r in audit)*len(cfg["families"])),
        empirical_fits=0, occurrence_changes=0, background_changes=0, held_target_changes=0,
        output_hashes={str(p.relative_to(out)): sha256(p) for p in out.rglob("*") if p.is_file()}))


if __name__ == "__main__":
    main()
