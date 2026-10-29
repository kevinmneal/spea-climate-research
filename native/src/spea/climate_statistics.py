"""Fixed summaries of twelve climate-normal months, not weather variability."""
from pathlib import Path

import numpy as np
import pandas as pd

from .genus_model import LOG_KEYS
from .seasonal_model import MONTHLY, VARIABILITY, monthly_features

SHAPE = [
    "temperature_median_minus_mean_c", "temperature_iqr_c", "temperature_skew",
    "precipitation_median_over_mean", "precipitation_iqr_over_mean", "precipitation_skew",
]
SEASONAL = [
    "precipitation_normalized_entropy", "wettest_circular_3month_fraction",
    "driest_circular_3month_fraction", "precipitation_harmonic_sin",
    "precipitation_harmonic_cos", "temperature_harmonic_sin_c",
    "temperature_harmonic_cos_c", "precipitation_temperature_offset_c",
]
STATISTICS = SHAPE + SEASONAL
BIOCLIM = ["log1p_" + k if k in LOG_KEYS else k for k in [f"bio{i:02}" for i in range(1, 20)]]
UNION = BIOCLIM + MONTHLY + STATISTICS
FAMILIES = {
    "bioclim": BIOCLIM,
    "bioclim_statistics": BIOCLIM + STATISTICS,
    "compact_summary": VARIABILITY + STATISTICS,
    "monthly": MONTHLY,
    "monthly_statistics": MONTHLY + STATISTICS,
}


def _population_skew(values):
    centered = values - values.mean(axis=1, keepdims=True)
    scale = np.max(np.abs(centered), axis=1)
    constant = scale == 0
    z = np.zeros_like(centered)
    np.divide(centered, scale[:, None], out=z, where=~constant[:, None])
    m2 = np.mean(z * z, axis=1)
    result = np.zeros(len(values))
    np.divide(np.mean(z ** 3, axis=1), m2 ** 1.5, out=result, where=~constant)
    return result, constant


def statistical_features(precipitation, temperature):
    """Return fourteen physical/normalized summaries and constant-row counts.

    Equal calendar-month weights, population moment skewness, linear empirical
    quantiles. Exact constant rows have skew zero; no small-variance threshold
    removes legitimate variation. January has angle zero. Circular windows
    include December/January. Zero precipitation terms contribute zero entropy.
    """
    p, t = np.asarray(precipitation, dtype=float), np.asarray(temperature, dtype=float)
    if p.ndim != 2 or p.shape[1] != 12 or p.shape != t.shape:
        raise ValueError("Aligned rows of twelve monthly precipitation/temperature values required")
    if not np.isfinite(p).all() or not np.isfinite(t).all() or np.any(p < 0):
        raise ValueError("Finite temperatures and nonnegative precipitation required")
    total = p.sum(axis=1)
    if np.any(total <= 0):
        raise ValueError("Positive annual precipitation required; never drop cells")
    mean_p, mean_t = total / 12, t.mean(axis=1)
    q = p / total[:, None]
    anomaly = t - mean_t[:, None]
    tq = np.quantile(t, [.25, .5, .75], axis=1, method="linear")
    pq = np.quantile(p, [.25, .5, .75], axis=1, method="linear")
    ts, tc = _population_skew(t)
    ps, pc = _population_skew(p)
    logq = np.zeros_like(q)
    np.log(q, out=logq, where=q > 0)
    entropy = -(q * logq).sum(axis=1) / np.log(12)
    quarters = q + np.roll(q, -1, axis=1) + np.roll(q, -2, axis=1)
    theta = 2 * np.pi * np.arange(12) / 12
    s, c = np.sin(theta), np.cos(theta)
    frame = np.column_stack([
        tq[1] - mean_t, tq[2] - tq[0], ts,
        pq[1] / mean_p, (pq[2] - pq[0]) / mean_p, ps,
        entropy, quarters.max(axis=1), quarters.min(axis=1),
        q @ s, q @ c, anomaly @ s / 6, anomaly @ c / 6,
        (q * anomaly).sum(axis=1),
    ])
    if not np.isfinite(frame).all():
        raise ValueError("Nonfinite statistical summary")
    return pd.DataFrame(frame.astype("float32"), columns=STATISTICS), {
        "rows": len(p), "exact_constant_temperature_rows": int(tc.sum()),
        "exact_constant_precipitation_rows": int(pc.sum()),
    }


def subset_borders(source, destination, union_names, chosen_names):
    """Remap crossed feature subsets, preserving shared training-only cutoffs."""
    source, destination = Path(source), Path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    if len(set(union_names)) != len(union_names) or len(set(chosen_names)) != len(chosen_names):
        raise ValueError("Unique feature names required")
    positions = {name: i for i, name in enumerate(chosen_names)}
    if not set(positions) <= set(union_names):
        raise ValueError("Unknown chosen feature")
    lines = []
    for line in source.read_text().splitlines():
        fields = line.split("\t")
        i = int(fields[0])
        if not 0 <= i < len(union_names):
            raise ValueError("Border index outside union")
        name = union_names[i]
        if name in positions:
            fields[0] = str(positions[name])
            lines.append((positions[name], "\t".join(fields)))
    lines.sort(key=lambda item: item[0])
    destination.write_text("\n".join(line for _, line in lines) + "\n")


class StatisticsReader:
    """Read existing aligned physical arrays without retaining global mmaps."""
    def __init__(self, root, bio_registry, monthly_registry):
        self.root = Path(root)
        if bio_registry["window_grid"] != monthly_registry["window_grid"]:
            raise ValueError("Different climate grids")
        self.shape = tuple(bio_registry["window_grid"]["shape"])
        by_key = {f["key"]: f for f in bio_registry["fields"] + monthly_registry["fields"]}
        keys = [f"bio{i:02}" for i in range(1, 20)]
        keys += [f"{v}_{m:02}" for v in ["pr", "tas"] for m in range(1, 13)]
        self.paths = {key: self.root / by_key[key]["physical_path"] for key in keys}

    def _read(self, key, ids):
        source = np.load(self.paths[key], mmap_mode="r")
        try:
            if source.shape != self.shape:
                raise ValueError("Climate source shape changed")
            result = np.asarray(source.ravel()[ids], dtype="float64")
        finally:
            source._mmap.close()
        return result

    def read(self, native):
        ids = np.asarray(native, dtype="int64")
        if ids.ndim != 1 or np.any(ids < 0) or np.any(ids >= np.prod(self.shape)):
            raise ValueError("Invalid native indices")
        bio = np.column_stack([self._read(f"bio{i:02}", ids) for i in range(1, 20)])
        if not np.isfinite(bio).all():
            raise ValueError("Missing BIOCLIM source values")
        for i, key in enumerate([f"bio{i:02}" for i in range(1, 20)]):
            if key in LOG_KEYS:
                if np.any(bio[:, i] < 0):
                    raise ValueError("Negative BIOCLIM precipitation")
                bio[:, i] = np.log1p(bio[:, i])
        p = np.column_stack([self._read(f"pr_{m:02}", ids) for m in range(1, 13)])
        t = np.column_stack([self._read(f"tas_{m:02}", ids) for m in range(1, 13)])
        stats, audit = statistical_features(p, t)
        result = pd.concat([pd.DataFrame(bio.astype("float32"), columns=BIOCLIM),
                            monthly_features(p, t), stats], axis=1)
        if list(result.columns) != UNION:
            raise ValueError("Feature order changed")
        return result, audit
