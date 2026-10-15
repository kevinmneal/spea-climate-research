#!/usr/bin/env python3
"""Plot published aggregate scores; no model fitting or native arrays."""
import argparse
import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
UNITS = [("intermontana", "S. intermontana"), ("bombifrons", "S. bombifrons"),
         ("multiplicata", "S. multiplicata"), ("hammondii_north", "Northern hammondii"),
         ("spea_south", "Southern hammondii")]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "replay_figure.png")
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError("Choose a new output path")
    with (ROOT / "data/matched/fold_metrics.csv").open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    slots = {(r["species"], r["sampling"], r["family"], r["fold"]): float(r["mean_block_log_score"]) for r in rows}
    colors = [plt.colormaps["viridis"](.22), plt.colormaps["viridis"](.69)]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.35), sharex=True, sharey=True)
    values = []
    for ax, control, title in zip(axes, ["variability", "bioclim"],
                                  ["Monthly − (annual + variability)", "Monthly − full BIOCLIM"]):
        for y, (species, label) in enumerate(UNITS):
            n = 5 if y < 3 else 3
            for arm, offset, color in zip(["equal_cells", "soft_geo_10km"], [-.15, .15], colors):
                v = [slots[(species, arm, "monthly", str(k))] - slots[(species, arm, control, str(k))] for k in range(n)]
                jitter = [(.08 * k / (n - 1)) - .04 for k in range(n)]
                ax.scatter(v, [y + offset + j for j in jitter], color=color, alpha=.6, s=24, edgecolors="none")
                ax.scatter(math.fsum(v)/n, y + offset, color=color, marker="D", s=68, edgecolors="white", linewidths=.7, zorder=3)
                values.extend(v)
        ax.axvline(0, color=".5", ls="--", lw=.8)
        ax.set_title(title, fontsize=11.5, loc="left", pad=12)
        ax.set_xlabel("Change in occupied-block log score (nats)", fontsize=10)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
        ax.grid(axis="x", color=".92", lw=.6)
        ax.set_axisbelow(True)
    span = max(abs(v) for v in values) * 1.1
    axes[0].set_xlim(-span, span)
    axes[0].set_ylim(4.5, -.6)
    axes[0].set_yticks(range(5), [label for _, label in UNITS], fontsize=10)
    handles = [plt.Line2D([], [], marker="D", color=c, ls="none", label=label)
               for c, label in zip(colors, ["Equal training-cell weights", "Soft spatial weights (10 km)"])]
    fig.suptitle("What do monthly climate profiles add?", x=.025, y=.975, ha="left", fontsize=17)
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(.018, .921), ncol=2, frameon=False, fontsize=10)
    fig.text(.025, .022,
             "Positive values favour monthly profiles. Circles: held folds; diamonds: equal-fold means.\n"
             "Same fixed shallow learner and held geography. BIOCLIM includes additional thermal information. Developmental comparisons.",
             color=".3", fontsize=9, va="bottom")
    fig.subplots_adjust(left=.19, right=.975, top=.76, bottom=.20, wspace=.14)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=180, facecolor="white")
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()
