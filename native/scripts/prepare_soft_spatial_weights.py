"""Freeze one soft weighting policy on original training cells before fitting."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from spea.spatial_thinning import inherited_cell_rows, geodesic_distances_km
from spea.spatial_weights import proximity_weights, normalize_presence_weights, kish_effective_n
from spea.provenance import ROOT, sha256, utcnow, write_json


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", default="config/soft_spatial_weighting_v1.json")
    args = ap.parse_args()
    path = ROOT / args.config
    cfg = json.loads(path.read_text())
    parent = json.loads((ROOT / "runs" / cfg["parent_run"] / "run.json").read_text())
    out = ROOT / "runs" / cfg["sampling_preparation_run"]
    out.mkdir(exist_ok=False); (out / "weights").mkdir()
    (out / "effective_config.json").write_bytes(path.read_bytes())
    (out / "source_snapshot.py").write_bytes(Path(__file__).read_bytes())
    inputs = {args.config: sha256(path)}

    def asset(run, name):
        p = ROOT / "runs" / run / name
        value = sha256(p)
        assert value == parent["input_hashes"][str(p.relative_to(ROOT))]
        inputs[str(p.relative_to(ROOT))] = value
        return p

    genus = pd.read_csv(asset(cfg["preparation_run"], "cells.csv.gz"))
    ns = pd.read_csv(asset(cfg["north_south_run"], "species_cells.csv"))
    audits, provenance = [], []
    arm = next(iter(cfg["sampling_arms"]))
    for species in cfg["species"]:
        is_ns = species in ["hammondii_north", "spea_south"]
        run = cfg["north_south_run"] if is_ns else cfg["preparation_run"]
        for fold in range(3 if is_ns else 5):
            prefix = f"{species}_fold{fold}"
            with np.load(asset(run, "splits/" + prefix + ".npz")) as split:
                selected = inherited_cell_rows(ns if is_ns else genus, split["train_cell_rows"])
                assert np.array_equal(selected.native_flat, split["train_native"])
                p = split["train_native"].copy()
            d = geodesic_distances_km(selected.longitude, selected.latitude)
            raw = proximity_weights(d, cfg["weighting"]["support_km"], cfg["weighting"]["alpha"])
            w = normalize_presence_weights(raw)
            assert (w > 0).all() and np.isclose(w.sum(), len(p))
            np.savez_compressed(out / "weights" / (prefix + ".npz"), eligible_native=p,
                raw_weights=raw, presence_weights=w)
            audits.append(dict(species=species, fold=fold, sampling=arm, train_cells=len(p),
                retained_cells=len(p), train_blocks50=selected.block50.nunique(),
                raw_min=float(raw.min()), raw_max=float(raw.max()),
                normalized_min=float(w.min()), normalized_max=float(w.max()),
                kish_effective_n=kish_effective_n(w), kish_fraction=kish_effective_n(w)/len(p),
                top_weight_share=float(w.max()/w.sum()), raw_total=float(raw.sum()), positive_mass=float(w.sum())))
            for field in ["countries", "source_ids"]:
                group_count, group_mass = {}, {}
                for values, mass in zip(selected[field], w):
                    for value in json.loads(values):
                        group_count[value] = group_count.get(value, 0) + 1
                        group_mass[value] = group_mass.get(value, 0.) + mass
                for value in group_count:
                    provenance.append(dict(species=species, fold=fold, field=field, category=value,
                        overlapping_cells=group_count[value], overlapping_positive_mass=group_mass[value],
                        baseline_cell_fraction=group_count[value]/len(p), weighted_fraction=group_mass[value]/len(p)))
    pd.DataFrame(audits).to_csv(out / "weight_audit.csv", index=False)
    pd.DataFrame(provenance).to_csv(out / "source_country_weight_audit.csv", index=False)
    write_json(out / "run.json", dict(status="completed_diagnostic_no_fits", finished_utc=utcnow(),
        empirical_fits=0, supported_slots=len(audits), unsupported_slots=0, input_hashes=inputs,
        script_sha256=sha256(Path(__file__)), weight_module_sha256=sha256(ROOT / "src/spea/spatial_weights.py"),
        output_hashes={str(p.relative_to(out)):sha256(p) for p in out.rglob("*") if p.is_file()}))
    print(json.dumps(dict(status="completed_diagnostic_no_fits", frozen_weight_vectors=len(audits), discarded_presences=0)))


if __name__ == "__main__":
    main()
