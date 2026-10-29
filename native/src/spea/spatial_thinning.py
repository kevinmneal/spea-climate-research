"""Training-only minimum-distance declustering of existing climate cells.

Cell-center separation is a sampling scenario, not coordinate precision or
proof of independent observations. This module never changes background area.
"""
import numpy as np
from pyproj import Geod


def inherited_cell_rows(cells, row_ids):
    """Resolve preserved source-row IDs in a filtered N/S or genus table."""
    if "source_cell_row" in cells.columns:
        if not cells.source_cell_row.is_unique:
            raise ValueError("Ambiguous source cell row namespace")
        return cells.set_index("source_cell_row", drop=False).loc[row_ids]
    return cells.iloc[row_ids]


def geodesic_distances_km(longitude, latitude):
    """Pairwise WGS84 distances, including high latitudes and the date line."""
    lon, lat = np.asarray(longitude, float), np.asarray(latitude, float)
    if (lon.ndim != 1 or lon.shape != lat.shape or not len(lon)
            or not np.isfinite(lon).all() or not np.isfinite(lat).all()
            or np.any(np.abs(lon) > 180) or np.any(np.abs(lat) > 90)):
        raise ValueError("Finite aligned WGS84 coordinates required")
    geod = Geod(ellps="WGS84")
    n = len(lon)
    result = np.empty((n, n), dtype=float)
    for a in range(0, n, 128):
        b = min(a + 128, n)
        x1, x2 = np.broadcast_arrays(lon[a:b, None], lon[None, :])
        y1, y2 = np.broadcast_arrays(lat[a:b, None], lat[None, :])
        result[a:b] = geod.inv(x1, y1, x2, y2)[2] / 1000
    np.fill_diagonal(result, 0.)
    return result


def thin_attempts(distances_km, minimum_km, seeds):
    """Remove a maximum-neighbor-degree member, with seeded random ties.

    Return every candidate mask and the first largest of the declared attempts.
    The selected subset is not guaranteed to have maximum possible cardinality.
    Equality at the declared radius is allowed; strictly closer pairs conflict.
    """
    d = np.asarray(distances_km, float)
    seeds = list(seeds)
    if (d.ndim != 2 or d.shape[0] != d.shape[1] or not len(d)
            or not np.isfinite(d).all() or np.any(d < 0)
            or not np.allclose(d, d.T, rtol=0, atol=1e-8)
            or np.any(np.diag(d) != 0) or not np.isfinite(minimum_km)
            or minimum_km <= 0 or not seeds or len(set(seeds)) != len(seeds)):
        raise ValueError("Symmetric nonnegative distances, positive radius and unique seeds required")
    neighbors = d < minimum_km
    np.fill_diagonal(neighbors, False)
    masks = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        active = np.ones(len(d), dtype=bool)
        degree = neighbors.sum(axis=1).astype(int)
        while degree.max() > 0:
            remove = int(rng.choice(np.flatnonzero(degree == degree.max())))
            active[remove] = False
            degree -= neighbors[remove] & active
            degree[~active] = -1
        masks.append(active)
    masks = np.asarray(masks)
    selected = int(np.argmax(masks.sum(axis=1)))
    return masks, selected


def minimum_retained_distance(distances_km, mask):
    ids = np.flatnonzero(mask)
    if len(ids) < 2:
        return None
    d = distances_km[np.ix_(ids, ids)]
    return float(d[np.triu_indices(len(ids), k=1)].min())
