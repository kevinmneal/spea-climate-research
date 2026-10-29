#!/usr/bin/env python3
"""Match saved four-family fits and summarize within-fold representation contrasts.

Reads metrics and small fit metadata only; never loads models or predicts.
"""
import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/private/tmp/spea-bioclim-reference-mpl")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from spea.climate_statistics import BIOCLIM
from spea.provenance import ROOT, sha256, utcnow, write_json
from spea.seasonal_model import ANNUAL, VARIABILITY, MONTHLY

NAMES = {"annual": ANNUAL, "variability": VARIABILITY,
         "bioclim": BIOCLIM, "monthly": MONTHLY}
LABELS = ["Intermontana", "Bombifrons", "Multiplicata",
          "Northern hammondii", "Southern species"]
COUNTS = ["train_cells", "held_cells", "held_blocks", "target_cells"]
KEYS = ["species", "fold", "sampling", "family"]


def load_source(name):
    path = ROOT / "runs" / name
    record = json.loads((path / "run.json").read_text())
    if record["status"] != "completed":
        raise ValueError(f"Incomplete source: {name}")
    digest = sha256(path / "metrics.csv")
    if digest != record["output_hashes"]["metrics.csv"]:
        raise ValueError(f"Changed saved metrics: {name}")
    return record, pd.read_csv(path / "metrics.csv"), digest


def summarize(frame, groups, fields):
    rows = []
    for key, part in frame.groupby(groups, sort=False):
        row = dict(zip(groups, key))
        row["folds"] = len(part)
        for field in fields:
            row[field] = float(part[field].mean())
            row[field + "_min"] = float(part[field].min())
            row[field + "_max"] = float(part[field].max())
            row[field + "_positive_folds"] = int((part[field] > 0).sum())
        rows.append(row)
    return pd.DataFrame(rows)


def figure(contrasts, cfg, output):
    """All folds and both arms for the two interpretively distinct contrasts."""
    fields = cfg["metrics"]
    titles = ["Occupied-block log score (nats)", "Recorded-cell log score (nats)",
              "Source ROC AUC", "Source AP (artificial 1:1 mass)"]
    pairs = ["monthly_minus_variability", "monthly_minus_bioclim"]
    colors = [plt.colormaps["viridis"](.22), plt.colormaps["viridis"](.68)]
    fig, axes = plt.subplots(4, 2, figsize=(12.4, 13.3), sharey=True)
    for i, (field, title) in enumerate(zip(fields, titles)):
        values = contrasts[contrasts.contrast.isin(pairs)][field].to_numpy()
        span = max(abs(values).max(), .01) * 1.09
        for j, pair in enumerate(pairs):
            ax = axes[i, j]
            for y, species in enumerate(cfg["species"]):
                p = contrasts[(contrasts.species == species) & (contrasts.contrast == pair)]
                equal = p[p.sampling == "equal_cells"].sort_values("fold")[field].to_numpy()
                soft = p[p.sampling == "soft_geo_10km"].sort_values("fold")[field].to_numpy()
                jitter = np.linspace(-.042, .042, len(equal))
                for k in range(len(equal)):
                    ax.plot([equal[k], soft[k]], [y - .15 + jitter[k], y + .15 + jitter[k]],
                            color=".8", lw=.6, zorder=1)
                for v, offset, color in zip([equal, soft], [-.15, .15], colors):
                    ax.scatter(v, y + offset + jitter, s=23, color=color, alpha=.6,
                               edgecolor="none", zorder=2)
                    ax.scatter(v.mean(), y + offset, s=60, marker="D", color=color,
                               edgecolor="white", linewidth=.6, zorder=3)
            ax.axvline(0, color=".5", ls="--", lw=.8)
            ax.set_xlim(-span, span)
            ax.set_ylim(4.5, -.5)
            ax.set_yticks(range(5), LABELS, fontsize=9)
            ax.set_xlabel(title, fontsize=10)
            ax.spines[["top", "right", "left"]].set_visible(False)
            ax.tick_params(axis="y", length=0)
            ax.grid(axis="x", color=".92", lw=.6)
            if i == 0:
                ax.set_title(["Monthly − annual plus variability",
                              "Monthly − full BIOCLIM"][j], fontsize=12, pad=11)
    handles = [plt.Line2D([], [], color=c, marker="D", ls="none", label=t)
               for c, t in zip(colors, ["Equal training-cell weights", "Soft training weights, 10 km"])]
    fig.suptitle("Monthly profiles against simpler controls and full BIOCLIM",
                 x=.025, y=.982, ha="left", fontsize=16)
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.025, .954),
               ncol=2, frameon=False, fontsize=10)
    fig.text(.025, .021,
             "Positive: monthly scores higher. Circles: every fold; diamonds: equal-fold means; lines: the same fold across arms.\n"
             "Same saved shallow procedure and geography. BIOCLIM is a non-nested reference with additional thermal information.\n"
             "Primary metric is occupied-block log score; secondary diagnostics can disagree. These are development comparisons.",
             fontsize=9.3, color=".3", va="bottom")
    fig.subplots_adjust(left=.18, right=.975, top=.884, bottom=.105, hspace=.38, wspace=.16)
    fig.savefig(output / "seasonal_bioclim_reference.png", dpi=170, facecolor="white")
    with plt.rc_context({"svg.fonttype": "none"}):
        fig.savefig(output / "seasonal_bioclim_reference.svg", facecolor="white")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    config_path = ROOT / args.config
    cfg = json.loads(config_path.read_text())
    source_names = list(dict.fromkeys([cfg["matched_run"], *cfg["bioclim_sources"].values()]))
    sources = {name: load_source(name) for name in source_names}
    matched = sources[cfg["matched_run"]][1].copy()
    matched["metrics_source_run"] = cfg["matched_run"]
    parts = [matched]
    for arm, name in cfg["bioclim_sources"].items():
        part = sources[name][1]
        part = part[(part.sampling == arm) & (part.family == "bioclim")].copy()
        part["metrics_source_run"] = name
        part["model_source_run"] = name
        part["reused"] = True
        parts.append(part)
    metrics = pd.concat(parts, ignore_index=True)
    expected = {(s, k, a, f) for s, n in zip(cfg["species"], cfg["fold_counts"])
                for k in range(n) for a in cfg["sampling_arms"] for f in cfg["families"]}
    if len(metrics) != len(expected) or set(metrics[KEYS].itertuples(index=False, name=None)) != expected:
        raise ValueError("Incomplete or duplicate family/arm/fold slots")
    if not np.isfinite(metrics[cfg["metrics"]].to_numpy()).all():
        raise ValueError("Nonfinite comparison metric")
    for _, part in metrics.groupby(["species", "fold"]):
        if any(part[field].nunique() != 1 for field in COUNTS):
            raise ValueError("Unmatched source or target counts")

    metadata, hashes, alignment = {}, {}, []
    for row in metrics.itertuples(index=False):
        path = ROOT / "runs" / row.model_source_run / "models" / f"{row.stem}_metadata.json"
        record = sources[row.model_source_run][0]
        rel = str(path.relative_to(ROOT / "runs" / row.model_source_run))
        digest = sha256(path)
        if digest != record["output_hashes"][rel]:
            raise ValueError(f"Changed metadata: {rel}")
        hashes[str(path.relative_to(ROOT))] = digest
        m = json.loads(path.read_text())
        metadata[row.species, row.fold, row.sampling, row.family] = m
        if m["feature_names"] != NAMES[row.family] or row.features != len(NAMES[row.family]):
            raise ValueError("Feature order or dimension changed")
        seed = cfg["model_seed"] + 100 * cfg["species"].index(row.species) + row.fold
        if m["seed"] != seed or m["candidate"] != cfg["model"] or m["objective"] != "bce":
            raise ValueError("Unmatched fixed fitting procedure")
        if not m["no_held_data_in_fit"] or m["statistical_early_stopping"]:
            raise ValueError("Unexpected held-data use")
        if m["presence_mass"] != row.train_cells or not np.isclose(
                m["effective_leaf_damping"], .03 * row.train_cells, rtol=1e-12):
            raise ValueError("Unmatched presence mass or effective damping")
        m["saved_split_digest"] = record["input_hashes"][m["source_split"]]
        if m.get("source_split_sha256", m["saved_split_digest"]) != m["saved_split_digest"]:
            raise ValueError("Split receipt mismatch")

    for (species, fold, arm), part in metrics.groupby(["species", "fold", "sampling"], sort=False):
        ms = [metadata[species, fold, arm, f] for f in cfg["families"]]
        if len({m["source_split"] for m in ms}) != 1 or len({m["saved_split_digest"] for m in ms}) != 1:
            raise ValueError("Different saved source memberships")
        weights = {m.get("frozen_weight_vector_sha256") for m in ms}
        if len(weights) != 1 or (arm == "soft_geo_10km" and None in weights):
            raise ValueError("Different training weighting vectors")
        row = dict(species=species, fold=int(fold), sampling=arm,
                   source_split=ms[0]["source_split"], saved_split_sha256=ms[0]["saved_split_digest"],
                   seed=ms[0]["seed"], fixed_procedure_matched=True,
                   frozen_weight_vector_sha256=next(iter(weights)))
        row.update({field: int(part[field].iloc[0]) for field in COUNTS})
        alignment.append(row)

    # The inherited common-column cutoffs were verified in the source runs.
    # Their saved receipts must also remain identical across these weighting arms.
    for species, count in zip(cfg["species"], cfg["fold_counts"]):
        for fold in range(count):
            for family in cfg["families"]:
                e = metadata[species, fold, "equal_cells", family]
                s = metadata[species, fold, "soft_geo_10km", family]
                if e["borders_sha256"] != s["borders_sha256"]:
                    raise ValueError("Different quantization receipts across arms")

    fold_rows = []
    for (species, fold, arm), part in metrics.groupby(["species", "fold", "sampling"], sort=False):
        values = part.set_index("family")
        for contrast, (a, b) in cfg["contrasts"].items():
            row = dict(species=species, fold=int(fold), sampling=arm, contrast=contrast)
            row.update({field: float(values.loc[a, field] - values.loc[b, field])
                        for field in cfg["metrics"]})
            fold_rows.append(row)
    contrasts = pd.DataFrame(fold_rows)
    summary = summarize(contrasts, ["species", "sampling", "contrast"], cfg["metrics"])
    family_summary = metrics.groupby(["species", "sampling", "family"], sort=False)[cfg["metrics"]].mean().reset_index()

    registry_path = ROOT / cfg["environment_registry"]
    registry = json.loads(registry_path.read_text())
    fields = {field["key"]: field for field in registry["fields"]}
    lineage = []
    for family in cfg["families"]:
        source_keys = ([f"bio{i:02}" for i in range(1, 20)] if family == "bioclim"
                       else [f"{v}_{i:02}" for v in ["tas", "pr"] for i in range(1, 13)])
        lineage.append(dict(family=family, dimensions=len(NAMES[family]),
            feature_names=json.dumps(NAMES[family]),
            input_physical_paths=json.dumps([fields[k]["physical_path"] for k in source_keys]),
            derivation=("Existing provider BIOCLIM rasters; log1p BIO12/13/14/16/17/18/19"
                        if family == "bioclim" else "Deterministic prefix of monthly_features in src/spea/seasonal_model.py"),
            nesting=("Non-nested reference; includes diurnal/minimum/maximum thermal information"
                     if family == "bioclim" else "annual (2) subset variability (4) subset monthly (28)")))
    if any(k.startswith(("tasmin_", "tasmax_")) for k in fields):
        raise ValueError("Update lineage: monthly registry now contains min/max temperature fields")

    output = ROOT / cfg["output"]
    run = ROOT / "runs" / cfg["run"]
    if output.exists() or run.exists():
        raise FileExistsError("Use new output/run paths; preserve earlier packets")
    output.mkdir(parents=True)
    run.mkdir(parents=True)
    (run / "summary_script.py").write_text(Path(__file__).read_text())
    write_json(run / "config.json", cfg)
    write_json(run / "source_metadata_hashes.json", hashes)
    metrics.to_csv(output / "fold_metrics.csv", index=False)
    family_summary.to_csv(output / "family_summary.csv", index=False)
    contrasts.to_csv(output / "contrast_folds.csv", index=False)
    summary.to_csv(output / "contrast_summary.csv", index=False)
    pd.DataFrame(alignment).to_csv(output / "source_alignment.csv", index=False)
    pd.DataFrame(lineage).to_csv(output / "family_lineage.csv", index=False)
    table = ["| Unit | Weight arm | Monthly − variability | Monthly − full BIOCLIM |",
             "|---|---|---:|---:|"]
    for species, label in zip(cfg["species"], LABELS):
        for arm in cfg["sampling_arms"]:
            p = summary[(summary.species == species) & (summary.sampling == arm)].set_index("contrast")
            entries = []
            for contrast in ["monthly_minus_variability", "monthly_minus_bioclim"]:
                r = p.loc[contrast]
                entries.append(f"{r.mean_block_log_score:+.4f} ({int(r.mean_block_log_score_positive_folds)}/{int(r.folds)})")
            table.append(f"| {label} | {'Equal' if arm == 'equal_cells' else 'Soft 10 km'} | {entries[0]} | {entries[1]} |")
    (output / "main_table.md").write_text("\n".join(table) + "\n")
    figure(contrasts, cfg, output)
    write_json(output / "effective_config.json", cfg)
    provenance = dict(status="aggregated_pending_figure_inspection", created_utc=utcnow(),
        config_sha256=sha256(config_path), summary_script_sha256=sha256(Path(__file__)),
        source_runs={name: {"run_sha256": sha256(ROOT / "runs" / name / "run.json"),
                           "metrics_sha256": source[2]} for name, source in sources.items()},
        environment_registry_sha256=sha256(registry_path),
        slots=len(metrics), source_alignment_groups=len(alignment), contrasts=len(contrasts),
        new_fits=0, new_predictions=0, historical_projections=0,
        source_and_target_counts_match=True, saved_split_receipts_match=True,
        frozen_weight_vectors_match=True, fixed_procedure_matches=True,
        preserved_all_folds=True, evaluation_normalization="Each saved fit's own held-area normalizer",
        full_bioclim_nested_in_monthly=False, independent_confirmation=False,
        original_source_hashes_recomputed=False, raw_models_or_arrays_loaded=False)
    write_json(output / "provenance.json", provenance)
    write_json(run / "run.json", provenance)
    print(json.dumps({"status": "aggregated", "slots": len(metrics),
                      "alignment_groups": len(alignment), "fold_contrasts": len(contrasts),
                      "output": cfg["output"], "new_fits": 0}), flush=True)


if __name__ == "__main__":
    main()
