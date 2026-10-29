"""Portable main-figure rebuild from packaged aggregates (requires matplotlib/numpy).

Usage: python build_main_figures.py --output-dir /tmp/spea-figure-replay
S1 cannot be rebuilt here: native supports remain in the original local archive.
"""
import argparse
import csv
import hashlib
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", tempfile.mkdtemp(prefix="spea-manuscript-mpl-"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator

ROOT = Path(__file__).resolve().parents[1]
AUDIT = "source_data/"
INPUTS = set()
NAVY = matplotlib.colormaps["viridis"](0.25)
TEAL = matplotlib.colormaps["viridis"](0.55)
INK = "#263238"
parser = argparse.ArgumentParser()
parser.add_argument("--figure", choices=["figure1_controlled_responses", "figure2_seasonal_representations", "figure3_mapped_water"])
parser.add_argument("--output-dir", type=Path, required=True, help="Temporary or new figure output directory; package figures are never overwritten")
args = parser.parse_args()
requested_figure = args.figure
OUT = args.output_dir.resolve()
assert OUT not in [ROOT / "figures", ROOT / "source_data", ROOT / "replay"], "Choose a separate output directory"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.titlesize": 10, "axes.labelsize": 9,
    "xtick.labelsize": 9, "ytick.labelsize": 9, "legend.fontsize": 9,
    "axes.spines.top": False, "axes.spines.right": False,
    "axes.spines.left": False, "axes.edgecolor": "#7b858a",
    "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.size": 0, "xtick.major.size": 3,
    "svg.fonttype": "none", "savefig.facecolor": "white",
    "pdf.fonttype": 42,
})


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read(path):
    INPUTS.add(path)
    with (ROOT / path).open(newline="") as f:
        return list(csv.DictReader(f))


def style(ax, xlim):
    ax.set_xlim(*xlim)
    ax.axvline(0, color="#667278", lw=0.8, zorder=1)
    ax.grid(axis="x", color="#e6e9ea", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    ax.tick_params(axis="y", pad=5)


def ticks(ax, step):
    low, high = ax.get_xlim()
    ax.set_xticks(np.arange(np.ceil(low / step), np.floor(high / step) + 1) * step)


def save(fig, stem):
    if requested_figure is not None and requested_figure != stem:
        plt.close(fig)
        return
    assert abs(fig.get_figwidth() - 6.5) < 1e-10
    fig.savefig(OUT / (stem + ".png"), dpi=300)
    fig.savefig(OUT / (stem + ".svg"))
    # Inspect text extents at the final physical figure dimensions.
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    outside = []
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text():
            continue
        assert text.get_fontsize() >= 9, (stem, text.get_text(), text.get_fontsize())
        box = text.get_window_extent(renderer)
        if box.width and box.height and (box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.width + 1 or box.y1 > fig.bbox.height + 1):
            outside.append(text.get_text())
    assert not outside, (stem, "text outside canvas", outside)
    dimensions[stem] = {"width_inches": fig.get_figwidth(), "height_inches": fig.get_figheight(), "png_dpi": 300, "minimum_text_points": 9}
    plt.close(fig)


dimensions = {}
audit_manifest_before = sha(ROOT / AUDIT / "source_manifest.json")
manifest_files_before = {}
for name, identity in json.loads((ROOT / AUDIT / "source_manifest.json").read_text())["packaged_files"].items():
    actual = sha(ROOT / AUDIT / name)
    assert actual == identity["sha256"], name
    manifest_files_before[name] = actual

# Figure 1: convex reference only; nonlinear numerical history is omitted.
targets = read(AUDIT + "support_primary_structured_targets.csv")
structures = read(AUDIT + "table1_controlled_responses.csv")
assert len(targets) == 6 and len(structures) == 4
fig, axes = plt.subplots(1, 2, figsize=(6.5, 4.05), gridspec_kw={"width_ratios": [1.12, 1]})
fig.subplots_adjust(left=0.245, right=0.98, top=0.89, bottom=0.24, wspace=0.77)
ax = axes[0]
names = ["S", "L", "C", "CL"]
offsets = {("North", "0"): -0.21, ("South", "0"): -0.13, ("North", "1"): -0.05, ("South", "1"): 0.05, ("North", "2"): 0.13, ("South", "2"): 0.21}
markers = {"0": "o", "1": "s", "2": "^"}
for i, name in enumerate(names):
    for r in targets:
        ax.scatter(float(r[name]), i + offsets[(r["group"], r["fold"])], s=25,
                   color=NAVY if r["group"] == "North" else TEAL,
                   marker=markers[r["fold"]], linewidths=0.45,
                   edgecolors="white", zorder=3)
    saved_mean = float(next(r for r in structures if r["structure"] == name)["score_nats_relative_recording_density"])
    ax.scatter(saved_mean, i, marker="D", s=36, color=INK, edgecolors="white", linewidths=0.6, zorder=4)
ax.set_yticks(range(4), ["S: shared", "L: shared + unit", "C: shared + gradual", "CL: gradual + unit"])
ax.set_ylim(3.5, -0.5)
style(ax, (-0.48, 0.62))
ticks(ax, 0.25)
ax.set_xlabel("Log density score (nats)")
ax.set_title("A  Absolute scores", loc="left", pad=10)

ax = axes[1]
ordered = sorted(targets, key=lambda r: (int(r["fold"]), r["group"] != "North"))
for i, r in enumerate(ordered):
    value = float(r["CL_minus_C"])
    color = NAVY if r["group"] == "North" else TEAL
    ax.plot([0, value], [i, i], color=color, lw=1.2, alpha=0.55, zorder=2)
    ax.scatter(value, i, s=31, color=color, marker=markers[r["fold"]], edgecolors="white", linewidths=0.45, zorder=3)
primary_mean = float(next(r for r in structures if r["structure"] == "CL")["mean_gain_nats"])
ax.axhline(5.65, color="#ccd2d4", lw=0.7)
ax.scatter(primary_mean, 6.3, marker="D", s=41, color=INK, edgecolors="white", linewidths=0.6, zorder=4)
labels = [r["group"] + " " + str(int(r["fold"]) + 1) for r in ordered] + ["Mean"]
ax.set_yticks(list(range(6)) + [6.3], labels)
ax.set_ylim(6.9, -0.6)
style(ax, (-0.11, 0.14))
ticks(ax, 0.10)
ax.set_xlabel("CL − C (nats)")
ax.set_title("B  Added unit response", loc="left", pad=10)
handles = [Line2D([], [], marker="o", ls="", color=NAVY, label="North"),
           Line2D([], [], marker="o", ls="", color=TEAL, label="South"),
           Line2D([], [], marker="D", ls="", color=INK, label="Mean"),
           *[Line2D([], [], marker=markers[str(i)], ls="", color=INK, label=f"Region {i + 1}") for i in range(3)]]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.54, 0.018), ncol=3, frameon=False, columnspacing=1.9, handletextpad=0.55)
save(fig, "figure1_controlled_responses")

# Figure 2: every fold and saved equal-fold mean; both evaluation targets.
folds = read(AUDIT + "seasonal_main_figure_fold_contrasts_long.csv")
summaries = read(AUDIT + "support_seasonal_all_metric_contrasts.csv")
assert len(folds) == 168
units = ["intermontana", "bombifrons", "multiplicata", "hammondii_north", "spea_south"]
labels = [r"$S.\ intermontana$", r"$S.\ bombifrons$", r"$S.\ multiplicata$", "North", "South"]
arm_colors = {"equal_cells": NAVY, "soft_geo_10km": TEAL}
arm_offsets = {"equal_cells": -0.12, "soft_geo_10km": 0.12}
mean_offsets = {"equal_cells": -0.35, "soft_geo_10km": 0.35}
contrasts = ["monthly_minus_variability", "monthly_minus_bioclim"]
fig, axes = plt.subplots(2, 2, figsize=(6.5, 6.8), sharey=True)
fig.subplots_adjust(left=0.245, right=0.98, top=0.855, bottom=0.092, hspace=0.33, wspace=0.13)
count = 0
for row, evaluation in enumerate(["occupied_block", "occupied_cell"]):
    for col, contrast in enumerate(contrasts):
        ax = axes[row, col]
        for i, unit in enumerate(units):
            if i % 2 == 0:
                ax.axhspan(i - 0.49, i + 0.49, color="#f6f7f7", zorder=0)
            common = [r for r in folds if r["unit"] == unit and r["contrast"] == contrast and r["evaluation"] == evaluation]
            n = len(common) // 2
            jitter = {str(j): float(v) for j, v in enumerate(np.linspace(-0.045, 0.045, n))}
            for fold in jitter:
                pair = {r["training_weights"]: r for r in common if r["fold"] == fold}
                ax.plot([float(pair[a]["gain_nats"]) for a in arm_colors],
                        [i + arm_offsets[a] + jitter[fold] for a in arm_colors],
                        color="#bfc6c8", lw=0.55, alpha=0.65, zorder=1)
            for arm, color in arm_colors.items():
                rr = [r for r in common if r["training_weights"] == arm]
                ax.scatter([float(r["gain_nats"]) for r in rr],
                           [i + arm_offsets[arm] + jitter[r["fold"]] for r in rr],
                           s=20, color=color, marker="o", linewidths=0.35,
                           edgecolors="white", zorder=3)
                count += len(rr)
                saved = next(r for r in summaries if r["species"] == unit and r["sampling"] == arm and r["contrast"] == contrast)
                metric = "mean_block_log_score" if row == 0 else "mean_cell_log_score"
                ax.scatter(float(saved[metric]), i + mean_offsets[arm], color=color, marker="D", s=39,
                           edgecolors=INK, linewidths=0.65, zorder=4)
        ax.set_yticks(range(5), labels)
        ax.set_ylim(4.6, -0.6)
        style(ax, (-0.43, 0.52) if row == 0 else (-0.84, 0.87))
        ticks(ax, 0.2 if row == 0 else 0.4)
        if col == 1:
            ax.tick_params(axis="y", labelleft=False)
        if row == 1:
            ax.set_xlabel("Score increment (nats)")
        ax.set_title(f"{'ABCD'[row * 2 + col]}  " + ("Primary: block" if row == 0 else "Secondary: cell"), loc="left", pad=8)
assert count == 168
fig.text(0.417, 0.982, "Monthly − variability", ha="center", va="top", fontsize=10)
fig.text(0.807, 0.982, "Monthly − full BIOCLIM", ha="center", va="top", fontsize=10)
handles = [Line2D([], [], marker="o", color=NAVY, ls="", label="Equal weights"),
           Line2D([], [], marker="o", color=TEAL, ls="", label="Soft 10 km"),
           Line2D([], [], marker="o", color=INK, ls="", label="Fold"),
           Line2D([], [], marker="D", color=INK, ls="", label="Mean")]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.575, 0.948), ncol=4,
           frameon=False, columnspacing=1.15, handletextpad=0.4)
save(fig, "figure2_seasonal_representations")

# Figure 3: distinct observation units and scales, without comparing effect sizes.
water = read("source_data/occurrence_water_paired_gains.csv")
survey = read("source_data/survey_water_paired_gains.csv")
water_means = read(AUDIT + "table3a_presence_water_density.csv")
survey_mean = read(AUDIT + "table3b_survey_water_history.csv")
fig, axes = plt.subplots(1, 2, figsize=(6.5, 4.2))
fig.subplots_adjust(left=0.14, right=0.98, top=0.88, bottom=0.255, wspace=0.55)
ax = axes[0]
families = [("quadratic_ridge", NAVY, -0.14), ("catboost", TEAL, 0.14)]
for family, color, offset in families:
    rr = sorted([r for r in water if r["kind"] == "spatial" and r["family"] == family], key=lambda r: r["split"])
    assert len(rr) == 5
    for r in rr:
        i = int(r["split"].split("_")[-1]) - 1
        value = float(r["water_gain"])
        ax.scatter(value, i + offset, s=30, marker="o", facecolors="white" if i == 2 else color,
                   edgecolors=color, linewidths=1 if i == 2 else 0.4, zorder=3)
    saved = next(r for r in water_means if r["family"] == family)
    ax.scatter(float(saved["water_density_gain_nats"]), 5.25 + offset, marker="D", s=43,
               color=color, edgecolors=INK, linewidths=0.6, zorder=4)
ax.axhline(4.65, color="#ccd2d4", lw=0.7)
ax.set_yticks([0, 1, 2, 3, 4, 5.25], ["Block 1", "Block 2", "Block 3\n(n = 3)", "Block 4", "Block 5", "Mean"])
ax.set_ylim(5.85, -0.6)
style(ax, (-1.57, 0.78))
ticks(ax, 0.5)
ax.set_xlabel("Relative-density gain\n(nats)")
ax.set_title("A  Recorded occurrence", loc="left", pad=10)
ax = axes[1]
for r in survey:
    i = int(r["held_block"]) - 1
    value = float(r["joint_gain"])
    ax.plot([0, value], [i, i], color=INK, lw=1, alpha=0.5, zorder=2)
    ax.scatter(value, i, s=30, marker="o", facecolors="white" if i == 4 else INK,
               edgecolors=INK, linewidths=1 if i == 4 else 0.4, zorder=3)
ax.axhline(4.65, color="#ccd2d4", lw=0.7)
ax.scatter(float(survey_mean[0]["water_joint_gain_nats_per_cell_history"]), 5.25,
           marker="D", s=43, color=INK, edgecolors="white", linewidths=0.6, zorder=4)
ax.set_yticks([0, 1, 2, 3, 4, 5.25], ["Block 1", "Block 2", "Block 3", "Block 4", "Block 5\n(n = 5)", "Mean"])
ax.set_ylim(5.85, -0.6)
style(ax, (-0.085, 0.004))
ticks(ax, 0.02)
ax.set_xlabel("Joint-history gain\n(nats per cell history)")
ax.set_title("B  Survey histories", loc="left", pad=10)
handles = [Line2D([], [], marker="o", ls="", color=NAVY, label="Quadratic"),
           Line2D([], [], marker="o", ls="", color=TEAL, label="CatBoost"),
           Line2D([], [], marker="D", ls="", color=INK, label="Equal-block mean")]
fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.55, 0.015), ncol=3,
           frameon=False, columnspacing=1.6, handletextpad=0.45)
save(fig, "figure3_mapped_water")

assert sha(ROOT / AUDIT / "source_manifest.json") == audit_manifest_before
for name, expected in manifest_files_before.items():
    assert sha(ROOT / AUDIT / name) == expected, name

rendered_dimensions = dimensions.copy()
if requested_figure is not None and (OUT / "provenance.json").exists():
    prior = json.loads((OUT / "provenance.json").read_text())
    dimensions = {**prior["dimensions"], **dimensions}
provenance = {
    "status": "rendered_pending_visual_inspection",
    "scope": "Three main manuscript figures using saved aggregate scores and saved means only; no fitting, predictions, new uncertainty estimates, or native-data loads.",
    "script_sha256": sha(Path(__file__)),
    "input_sha256": {p: sha(ROOT / p) for p in sorted(INPUTS)},
    "output_sha256": {stem + ext: sha(OUT / (stem + ext)) for stem in dimensions for ext in [".png", ".svg"]},
    "dimensions": dimensions,
    "rendered_figures_this_execution": list(rendered_dimensions),
    "software": {"matplotlib": matplotlib.__version__, "numpy": np.__version__},
    "packaged_source_preservation": {"manifest_sha256_before_and_after": audit_manifest_before, "sealed_files_checked_before_and_after": len(manifest_files_before), "all_unchanged": True},
    "definitions": {
        "figure1": "Primary 591-cell convex reference only. S: shared response; L: shared plus North/South response; C: shared plus gradual geographic response; CL: gradual plus North/South response. Every structure includes the geographic main effect. All six group-region held scores and saved equally weighted means; paired CL-C increments have their own x scale. Regions 1-3 correspond to saved folds 0-2.",
        "figure2": "All 168 fold/contrast/evaluation observations from the audited 84 fold contrasts: five units, equal and soft 10-km training weights, monthly minus variability and monthly minus non-nested full BIOCLIM, primary occupied-block and secondary occupied-cell evaluation. Diamonds are SAVED equal-fold means. Fold pairing lines are visual links, not intervals. Same x scale within each evaluation row; different scales between rows.",
        "figure3": "Left: water-minus-climate relative recording-density scores for fixed quadratic and shallow CatBoost controls in five source blocks; block 3 has only three records. Right: water-minus-terrestrial complete surveyed-history log-probability per cell, five source blocks; block 5 has five cells. Open symbols mark the sparse blocks. Diamonds are saved equal-block means. The two panels use different score definitions and x scales and must not be pooled.",
    },
    "plotted_counts": {"figure1_absolute_scores": 24, "figure1_absolute_means": 4, "figure1_paired_gains": 6, "figure1_gain_mean": 1, "figure2_fold_points": count, "figure2_saved_mean_points": 40, "figure3_occurrence_block_points": 10, "figure3_survey_block_points": 5, "figure3_saved_mean_points": 3},
    "confidence_intervals": False, "new_fits": 0, "new_predictions": 0,
    "automated_layout_check": "All visible text at least 9 points, no text outside 6.5-inch canvas. Figures saved without tight-cropping so physical dimensions remain exact.",
}
(OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
print(json.dumps({"status": provenance["status"], "dimensions": dimensions, "output_directory": str(OUT)}, indent=2))
