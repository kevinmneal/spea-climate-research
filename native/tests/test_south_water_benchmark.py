"""Checks for scientific invariants in the public-point water comparison."""
import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

path = Path(__file__).resolve().parents[1] / "scripts/run_south_water_benchmark.py"
spec = importlib.util.spec_from_file_location("south_water_benchmark", path)
benchmark = importlib.util.module_from_spec(spec)
spec.loader.exec_module(benchmark)


def test_ridge_gradient_against_finite_differences_and_unpenalized_intercept():
    rng = np.random.default_rng(41)
    matrix = rng.normal(size=(12, 4))
    target = np.r_[np.ones(3), np.zeros(9)]
    weight = np.r_[np.ones(3), np.full(9, 1 / 3)]
    theta = rng.normal(size=5)
    args = (matrix, target, weight, 3, 0.03)
    _, gradient = benchmark.ridge_loss_gradient(theta, *args)
    step = 1e-6
    numerical = np.empty(5)
    for j in range(5):
        delta = np.zeros(5)
        delta[j] = step
        numerical[j] = (
            benchmark.ridge_loss_gradient(theta + delta, *args)[0]
            - benchmark.ridge_loss_gradient(theta - delta, *args)[0]
        ) / (2 * step)
    np.testing.assert_allclose(gradient, numerical, rtol=1e-6, atol=1e-8)
    _, without_penalty = benchmark.ridge_loss_gradient(theta, matrix, target, weight, 3, 0)
    assert gradient[0] == without_penalty[0]


def test_held_density_normalizes_and_is_invariant_to_intercept():
    presence = np.array([-2.0, 0.0, 1.0])
    available = np.array([-1.0, -0.5, 0.5, 2.0])
    scores, normalizer = benchmark.relative_scores(presence, available)
    assert np.isclose(np.exp(available - normalizer).mean(), 1.0)
    shifted, _ = benchmark.relative_scores(presence + 700, available + 700)
    np.testing.assert_allclose(scores, shifted, atol=1e-12)
    np.testing.assert_allclose(scores, np.log(np.exp(presence) / np.exp(available).mean()))


def test_guard_excludes_nearby_training_and_keeps_late_records_out():
    data = {
        "presence_xy": np.array([[-10000., 0], [-4900., 0], [15000., 0], [25000., 0], [-10000., 1]]),
        "background_xy": np.array([[-10000., 0], [0., 0], [15000., 0], [25000., 0], [-9000., 0]]),
        "presence_block": np.array([1, 2, 2, 2, 1]),
        "background_block": np.array([1, 1, 2, 2, 1]),
        "presence_year": np.array([2000, 2000, 2000, 2019, 2019]),
    }
    q = pd.DataFrame(data["background_xy"], columns=["X", "Y"])
    q["block"] = data["background_block"]
    cfg = {"temporal_boundary_year": 2017, "spatial_guard_m": 5000, "held_source_blocks": [1, 2]}
    splits = benchmark.build_splits(data, q, cfg)
    np.testing.assert_array_equal(np.flatnonzero(splits[0][3]), [2])
    np.testing.assert_array_equal(np.flatnonzero(splits[0][5]), [0])
    np.testing.assert_array_equal(np.flatnonzero(splits[-1][3]), [1, 2])
    np.testing.assert_array_equal(np.flatnonzero(splits[-1][5]), [3, 4])
    assert all(not split[3][data["presence_year"] >= 2017].any() for split in splits)


def test_water_has_no_quadratic_tail():
    climate = np.full((2, 19), 2.)
    water = np.array([[3., 4.], [5., 6.]])
    design = benchmark.quadratic_basis(np.column_stack((climate, water)))
    assert design.shape == (2, 40)
    np.testing.assert_array_equal(design[:, -2:], water)
