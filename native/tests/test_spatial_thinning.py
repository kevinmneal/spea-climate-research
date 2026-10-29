"""Distance, graph, held-role and source-mass checks for declustering."""
import numpy as np
import pandas as pd
import pytest
from pyproj import Geod

from spea.spatial_thinning import geodesic_distances_km, thin_attempts, minimum_retained_distance, inherited_cell_rows
from spea.native_objective_comparison import encode_objective


def test_filtered_species_table_resolves_source_ids_not_row_positions():
    cells = pd.DataFrame({"source_cell_row": [602, 715, 930], "native_flat": [110, 210, 310]})
    np.testing.assert_array_equal(inherited_cell_rows(cells, [930, 602]).native_flat, [310, 110])
    genus = cells.drop(columns="source_cell_row")
    np.testing.assert_array_equal(inherited_cell_rows(genus, [2, 0]).native_flat, [310, 110])


def test_geodesic_distance_across_dateline_at_high_latitude():
    geod = Geod(ellps="WGS84")
    lon2, lat2, _ = geod.fwd(179.99, 75., 90., 10000.)
    d = geodesic_distances_km([179.99, lon2], [75., lat2])
    assert lon2 < 0
    np.testing.assert_allclose(d, [[0, 10], [10, 0]], atol=1e-8)


def test_max_degree_removal_preserves_separated_star_leaves():
    # Euclidean distances suffice for this algorithmic graph test: central
    # point has four conflicts, while all outer leaves are mutually eligible.
    xy = np.array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1]], float)
    d = np.linalg.norm(xy[:, None] - xy[None], axis=2)
    masks, selected = thin_attempts(d, 1.1, range(10))
    np.testing.assert_array_equal(masks, np.tile([False, True, True, True, True], (10, 1)))
    assert selected == 0
    assert minimum_retained_distance(d, masks[selected]) >= 1.1


def test_seeds_reproduce_valid_masks_without_claiming_optimality():
    rng = np.random.default_rng(5)
    xy = rng.normal(size=(60, 2))
    d = np.linalg.norm(xy[:, None] - xy[None], axis=2)
    a, chosen = thin_attempts(d, .4, range(10))
    b, again = thin_attempts(d, .4, range(10))
    np.testing.assert_array_equal(a, b)
    assert chosen == again == np.argmax(a.sum(axis=1))
    for mask in a:
        assert minimum_retained_distance(d, mask) >= .4


def test_radius_equality_and_isolated_points_are_retained():
    d = np.abs(np.array([0., 5., 10.])[:, None] - np.array([0., 5., 10.])[None, :])
    masks, _ = thin_attempts(d, 5., [10, 11])
    assert masks.all()


def test_thinning_keeps_dropped_occupied_cells_on_area_quadrature():
    # Positive membership changes; the available-area measure and coordinates
    # do not. The loss normalization and effective damping scale with new n.
    q, area = np.arange(10), np.arange(1, 11, dtype=float)
    full = encode_objective("bce", [1, 3, 7], q, area, .03)
    thin = encode_objective("bce", [1, 7], q, area, .03)
    np.testing.assert_array_equal(thin.native[2:], q)
    assert 3 in thin.native[2:]
    np.testing.assert_allclose(thin.weight[2:] / thin.weight[2:].sum(), area / area.sum())
    assert thin.weight[:2].sum() == pytest.approx(2)
    assert thin.weight[2:].sum() == pytest.approx(2)
    assert thin.effective_leaf_damping == pytest.approx(.06)
    assert full.effective_leaf_damping == pytest.approx(.09)


@pytest.mark.parametrize("radius,seeds", [(0, [1]), (5, []), (5, [1, 1])])
def test_invalid_thinning_requests_raise(radius, seeds):
    with pytest.raises(ValueError):
        thin_attempts(np.array([[0., 1.], [1., 0.]]), radius, seeds)
