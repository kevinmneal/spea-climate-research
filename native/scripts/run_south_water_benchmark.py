#!/usr/bin/env python3
"""Run the frozen public-point climate/water comparison; never alter prior runs."""
from __future__ import annotations

import os

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(variable, "1")

import argparse
import hashlib
import json
import resource
import time
import traceback
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pandas as pd
from affine import Affine
from catboost import CatBoostClassifier
from pyproj import CRS, Transformer
from scipy.optimize import minimize
from scipy.spatial import cKDTree
from scipy.special import expit, logsumexp
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save_json(path: Path, value) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def ridge_loss_gradient(theta, matrix, target, weight, n_presence, penalty):
    eta = theta[0] + matrix @ theta[1:]
    loss = np.dot(weight, np.logaddexp(0, eta) - target * eta) / n_presence
    loss += 0.5 * penalty * np.dot(theta[1:], theta[1:])
    residual = weight * (expit(eta) - target) / n_presence
    gradient = np.r_[residual.sum(), matrix.T @ residual + penalty * theta[1:]]
    return float(loss), gradient


def quadratic_basis(features, climate_count=19):
    """Additive climate curvature; water has a single log-distance response."""
    climate = features[:, :climate_count]
    water = features[:, climate_count:]
    return np.column_stack((climate, climate ** 2, water))


def relative_scores(presence_eta, background_eta):
    normalizer = float(logsumexp(background_eta) - np.log(len(background_eta)))
    return np.asarray(presence_eta) - normalizer, normalizer


def prepare_data(cfg, output):
    p = pd.read_csv(ROOT / cfg["presence_table"])
    q_all = pd.read_csv(ROOT / cfg["background_table"])
    assert len(p) == 409 and len(q_all) == 28931 and p.ID.is_unique
    xml = ET.parse(ROOT / cfg["source_metadata"])
    assert xml.findtext(".//geodetic/horizdn") == "North_American_Datum_1983"
    assert xml.findtext(".//planar/planci/plandu") == "meters"
    albers = xml.find(".//albers")
    parallels = [float(x.text) for x in albers.findall("stdparll")]
    proj = (
        f"+proj=aea +lat_1={parallels[0]} +lat_2={parallels[1]} "
        f"+lat_0={albers.findtext('latprjo')} +lon_0={albers.findtext('longcm')} "
        f"+x_0={albers.findtext('feast')} +y_0={albers.findtext('fnorth')} "
        "+datum=NAD83 +units=m +no_defs"
    )
    transform = Transformer.from_crs(CRS.from_proj4(proj), "EPSG:4326", always_xy=True)
    reverse = Transformer.from_crs("EPSG:4326", CRS.from_proj4(proj), always_xy=True)
    registry = json.loads((ROOT / cfg["climate_registry"]).read_text())
    assert registry["baseline"] == "1981-2010"
    affine = Affine(*registry["window_grid"]["transform"][:6])
    fields = registry["fields"]
    assert [f["key"] for f in fields] == cfg["climate_keys"]

    arrays = {}
    exclusions = []
    for name, table in (("presence", p), ("background", q_all)):
        xy = table[["X", "Y"]].to_numpy(float)
        lon, lat = transform.transform(xy[:, 0], xy[:, 1])
        rx, ry = reverse.transform(lon, lat)
        assert np.max(np.abs(np.column_stack((rx, ry)) - xy)) < 0.01
        col, row = (~affine) * (lon, lat)
        row, col = np.floor(row).astype(int), np.floor(col).astype(int)
        height, width = registry["window_grid"]["shape"]
        assert ((row >= 0) & (row < height) & (col >= 0) & (col < width)).all()
        climate = np.column_stack([
            np.load(ROOT / field["physical_path"], mmap_mode="r")[row, col]
            for field in fields
        ]).astype(float)
        water = table[cfg["water_keys"]].to_numpy(float)
        assert np.isfinite(water).all() and (water >= 0).all()
        valid = np.isfinite(climate).all(axis=1)
        for index in np.flatnonzero(~valid):
            exclusions.append({"table": name, "source_row": int(index),
                               "reason": "incomplete_native_BIOCLIM_at_public_point"})
        if name == "presence":
            assert valid.all()
        arrays[name + "_features"] = np.column_stack((climate[valid], np.log1p(water[valid])))
        arrays[name + "_xy"] = xy[valid]
        arrays[name + "_source_row"] = np.flatnonzero(valid)
        arrays[name + "_native_row"] = row[valid]
        arrays[name + "_native_col"] = col[valid]
        arrays[name + "_block"] = table.block.to_numpy(int)[valid]
    arrays["presence_year"] = p.year.to_numpy(int)
    np.savez_compressed(output / "public_point_inputs.npz", **arrays)
    save_json(output / "excluded_public_rows.json", exclusions)
    q = q_all.iloc[arrays["background_source_row"]].reset_index(drop=True)
    assert len(q) == 28880 and len(exclusions) == 51
    save_json(output / "preflight.json", {
        "presence_rows": len(p), "complete_background_points": len(q),
        "excluded_background_points": len(exclusions), "projection_from_source_XML": proj,
        "climate_baseline": registry["baseline"], "climate_keys": cfg["climate_keys"],
        "water_transform": "log1p(distance_km)", "source_units_preserved": True,
        "private_survey_mapping_inferred": False,
    })
    return arrays, q_all


def build_splits(data, q_all, cfg):
    px, qx = data["presence_xy"], data["background_xy"]
    pb, qb = data["presence_block"], data["background_block"]
    early = data["presence_year"] < cfg["temporal_boundary_year"]
    splits = []
    guard = cfg["spatial_guard_m"]
    for block in cfg["held_source_blocks"]:
        hp, hq = early & (pb == block), qb == block
        held_xy = np.vstack((px[hp], q_all.loc[q_all.block == block, ["X", "Y"]].to_numpy(float)))
        tree = cKDTree(held_xy)
        p_distance, q_distance = tree.query(px)[0], tree.query(qx)[0]
        tp = early & (pb != block) & (p_distance > guard)
        tq = (qb != block) & (q_distance > guard)
        splits.append((f"spatial_{block}", "spatial", block, tp, tq, hp, hq))
    late = ~early
    tree = cKDTree(px[late])
    tp = early & (tree.query(px)[0] > guard)
    tq = tree.query(qx)[0] > guard
    splits.append(("temporal_guarded", "temporal", None, tp, tq, late, np.ones(len(qx), bool)))
    return splits


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    parser.add_argument("--run-id", required=True)
    args = parser.parse_args()
    cfg = json.loads((ROOT / args.config).read_text())
    output = ROOT / "runs" / args.run_id
    output.mkdir(exist_ok=False)
    for folder in ("splits", "models", "held", "training"):
        (output / folder).mkdir()
    save_json(output / "effective_config.json", cfg)
    state = {"status": "preflight", "model_fits": 0, "cpu_seconds": 0,
             "maximum_fits": cfg["maximum_fits"], "prior_model_changes": 0,
             "occurrence_admissions": 0, "interpretation": "relative_recording_density_only"}
    cpu = time.process_time()
    metrics, audits = [], []

    def flush():
        state["cpu_seconds"] = time.process_time() - cpu
        state["peak_resident_bytes"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        save_json(output / "run.json", state)
        save_json(output / "fit_audit.json", audits)
        pd.DataFrame(metrics).to_csv(output / "metrics.csv", index=False)

    try:
        for path, expected in cfg["input_hashes"].items():
            assert digest(ROOT / path) == expected, path
        data, q_all = prepare_data(cfg, output)
        features_p, features_q = data["presence_features"], data["background_features"]
        splits = build_splits(data, q_all, cfg)
        assert len(splits) == 6 and sum(int(s[5].sum()) for s in splits[:5]) == 324
        state["status"] = "fitting_frozen_comparison"
        for split_index, (split, kind, block, tp, tq, hp, hq) in enumerate(splits):
            pi, qi, hi, ai = [np.flatnonzero(x) for x in (tp, tq, hp, hq)]
            assert min(map(len, (pi, qi, hi, ai))) > 0
            assert not np.intersect1d(pi, hi).size
            if kind == "spatial":
                assert not np.intersect1d(qi, ai).size
            np.savez_compressed(output / "splits" / f"{split}.npz",
                                train_presence=pi, train_background=qi,
                                held_presence=hi, held_background=ai)
            n_presence = len(pi)
            target = np.r_[np.ones(len(pi)), np.zeros(len(qi))]
            weights = np.r_[np.ones(len(pi)), np.full(len(qi), n_presence / len(qi))]
            assert np.isclose(weights[target == 1].sum(), weights[target == 0].sum())
            # Only eligible training background supplies scales and centering.
            mean = features_q[qi].mean(axis=0)
            scale = features_q[qi].std(axis=0)
            scale[scale < 1e-8] = 1
            z_p, z_q = (features_p - mean) / scale, (features_q - mean) / scale
            for with_water in (False, True):
                width = 21 if with_water else 19
                keys = cfg["climate_keys"] + (cfg["water_keys"] if with_water else [])
                xp, xq = z_p[:, :width], z_q[:, :width]
                train = np.vstack((xp[pi], xq[qi]))
                for family in cfg["families"]:
                    stem = f"{split}_{family}_{'water' if with_water else 'climate'}"
                    assert state["model_fits"] < cfg["maximum_fits"]
                    state["model_fits"] += 1
                    flush()
                    if family == "quadratic_ridge":
                        bp, bq = quadratic_basis(xp), quadratic_basis(xq)
                        basis_mean, basis_scale = bq[qi].mean(axis=0), bq[qi].std(axis=0)
                        basis_scale[basis_scale < 1e-8] = 1
                        bp, bq = (bp - basis_mean) / basis_scale, (bq - basis_mean) / basis_scale
                        matrix = np.vstack((bp[pi], bq[qi]))
                        penalty = cfg["quadratic_ridge"]["penalty_per_presence"]
                        result = minimize(
                            ridge_loss_gradient, np.zeros(matrix.shape[1] + 1),
                            args=(matrix, target, weights, n_presence, penalty), jac=True,
                            method="L-BFGS-B", options={"maxiter": 2000, "gtol": 1e-9, "ftol": 1e-13},
                        )
                        loss, grad = ridge_loss_gradient(result.x, matrix, target, weights, n_presence, penalty)
                        assert result.success and np.max(np.abs(grad)) < 1e-5, result.message
                        theta = result.x
                        pe, qe = theta[0] + bp[hi] @ theta[1:], theta[0] + bq[ai] @ theta[1:]
                        np.savez_compressed(output / "models" / f"{stem}.npz", theta=theta,
                                            mean=mean[:width], scale=scale[:width],
                                            basis_mean=basis_mean, basis_scale=basis_scale)
                        train_eta = theta[0] + matrix @ theta[1:]
                        audit = {"optimizer_success": bool(result.success), "optimizer_iterations": int(result.nit),
                                 "gradient_max_abs": float(np.max(np.abs(grad))), "penalized_loss": loss}
                    else:
                        setting = cfg["catboost"]
                        model = CatBoostClassifier(
                            loss_function="Logloss", iterations=setting["iterations"], depth=setting["depth"],
                            learning_rate=setting["learning_rate"],
                            l2_leaf_reg=setting["leaf_damping_per_presence"] * n_presence,
                            border_count=128, feature_border_type="GreedyLogSum", bootstrap_type="No",
                            random_strength=0, boosting_type="Plain", leaf_estimation_method="Newton",
                            leaf_estimation_iterations=1, boost_from_average=False,
                            random_seed=cfg["seed"] + split_index, thread_count=cfg["threads"],
                            allow_writing_files=False, verbose=False, use_best_model=False,
                        )
                        model.fit(pd.DataFrame(train, columns=keys), target, sample_weight=weights)
                        assert model.tree_count_ == setting["iterations"]
                        model.save_model(str(output / "models" / f"{stem}.cbm"))
                        model.save_model(str(output / "models" / f"{stem}.json"), format="json")
                        np.savez_compressed(output / "models" / f"{stem}_scaler.npz",
                                            mean=mean[:width], scale=scale[:width])
                        pe = model.predict(pd.DataFrame(xp[hi], columns=keys), prediction_type="RawFormulaVal")
                        qe = model.predict(pd.DataFrame(xq[ai], columns=keys), prediction_type="RawFormulaVal")
                        train_eta = model.predict(pd.DataFrame(train, columns=keys), prediction_type="RawFormulaVal")
                        audit = {"tree_count": model.tree_count_, "parameters": model.get_all_params()}
                    assert np.isfinite(pe).all() and np.isfinite(qe).all()
                    scores, norm = relative_scores(pe, qe)
                    training_loss = float(np.dot(weights, np.logaddexp(0, train_eta) - target * train_eta) / weights.sum())
                    held_blocks = data["presence_block"][hi]
                    block_mean = float(np.mean([scores[held_blocks == b].mean() for b in np.unique(held_blocks)]))
                    auc_target = np.r_[np.ones(len(hi)), np.zeros(len(ai))]
                    auc_weight = np.r_[np.full(len(hi), 1 / len(hi)), np.full(len(ai), 1 / len(ai))]
                    auc = float(roc_auc_score(auc_target, np.r_[pe, qe], sample_weight=auc_weight))
                    metrics.append({"split": split, "kind": kind, "source_block": block,
                                    "family": family, "water": with_water, "train_presence": len(pi),
                                    "train_background": len(qi), "held_presence": len(hi), "held_background": len(ai),
                                    "mean_record_score": float(scores.mean()), "mean_source_block_score": block_mean,
                                    "background_log_normalizer": norm, "presence_background_auc": auc})
                    audits.append({"model": stem, "family": family, "water": with_water,
                                   "training_bce": training_loss, "presence_mass": n_presence,
                                   "background_mass": float(weights[target == 0].sum()),
                                   "scales_from_training_background_only": True, **audit})
                    np.savez_compressed(output / "held" / f"{stem}.npz", presence_eta=pe,
                                        background_eta=qe, scores=scores, normalizer=norm,
                                        presence_index=hi, background_index=ai)
                    np.savez_compressed(output / "training" / f"{stem}.npz", raw_eta=train_eta)
                    flush()
            print(f"Completed {split}: {len(pi)} training / {len(hi)} held records", flush=True)
        assert state["model_fits"] == cfg["maximum_fits"] == 24
        state["status"] = "completed_exploratory_public_point_benchmark"
        flush()
        manifest = {str(p.relative_to(output)): digest(p) for p in output.rglob("*") if p.is_file()}
        save_json(output / "manifest.json", manifest)
    except Exception:
        state["status"] = "failed_preserved"
        state["error"] = traceback.format_exc()
        flush()
        raise


if __name__ == "__main__":
    main()
