"""Soft geographic declustering; weights are not measured sampling effort."""
import numpy as np


def proximity_weights(distances_km, support_km, alpha=1.):
    """enmSdmX weightByDist policy, using explicit physical distances.

    A point with no neighbors inside support gets one raw unit; nearer
    neighbors gradually reduce its weight. No point gets zero weight.
    This describes geographic overrepresentation, not biological independence.
    """
    d = np.asarray(distances_km, float)
    if (d.ndim != 2 or d.shape[0] != d.shape[1] or not len(d)
            or not np.isfinite(d).all() or np.any(d < 0)
            or not np.allclose(d, d.T, rtol=0, atol=1e-8)
            or np.any(np.diag(d) != 0) or not np.isfinite(support_km)
            or support_km <= 0 or not np.isfinite(alpha) or alpha <= 0):
        raise ValueError("Symmetric finite distances and positive support/alpha required")
    influence = np.maximum(1. - d / support_km, 0.) ** alpha
    np.fill_diagonal(influence, 0.)
    return 1. / (1. + influence.sum(axis=1))


def normalize_presence_weights(raw):
    raw = np.asarray(raw, float)
    if raw.ndim != 1 or not len(raw) or not np.isfinite(raw).all() or np.any(raw <= 0):
        raise ValueError("Every retained presence requires positive finite weight")
    return raw * (len(raw) / raw.sum())


def kish_effective_n(weights):
    """Weight concentration diagnostic, not the number of independent sites."""
    w = normalize_presence_weights(weights)
    return float(w.sum() ** 2 / np.dot(w, w))
