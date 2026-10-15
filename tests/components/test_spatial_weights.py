import numpy as np
import pytest

from spea_snapshot.spatial_weights import proximity_weights, normalize_presence_weights, kish_effective_n


def test_coincident_records_share_one_raw_unit():
    assert np.allclose(proximity_weights(np.zeros((3, 3)), 10), 1/3)


def test_neighbor_influence_tapers_and_stops():
    for distance, expected in [(0, .5), (5, 2/3), (10, 1), (11, 1)]:
        assert np.allclose(proximity_weights([[0, distance], [distance, 0]], 10), expected)


def test_all_presences_survive_and_mass_is_preserved():
    d = np.array([[0, 1, 30], [1, 0, 29], [30, 29, 0]])
    raw = proximity_weights(d, 10)
    w = normalize_presence_weights(raw)
    assert np.all(w > 0) and np.isclose(w.sum(), 3)
    assert np.isclose(kish_effective_n(w), kish_effective_n(raw))
    assert 2 < kish_effective_n(w) < 3


def test_normalization_preserves_relative_weights():
    w = normalize_presence_weights([.1, .2, 1])
    assert np.allclose(w/w[0], [1, 2, 10])
    assert kish_effective_n(np.ones(8)) == 8


def test_invalid_distance_and_nonpositive_mass_rejected():
    with pytest.raises(ValueError):
        proximity_weights([[0, 2], [3, 0]], 10)
    with pytest.raises(ValueError):
        normalize_presence_weights([0, 1])
