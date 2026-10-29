"""Checks of latent-state marginalization and training-only feature transforms."""
from itertools import product

import numpy as np
import pytest
from scipy.special import expit

from spea.survey_water import HabitatTransform, history_terms, objective, observation_design


def example():
    x = np.array([[1., -.8], [1., .4], [1., 1.1]])
    w = np.array([[[1., .2], [1., -.5], [1., 1.2]],
                  [[1., 0.], [1., .8], [1., -.6]],
                  [[1., .4], [1., -.2], [1., .7]]])
    y = np.array([[1., 0., 0.], [0., 0., 0.], [0., 1., 1.]])
    m = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]], dtype=bool)
    return np.array([-.3, .6, -.9, .7]), x, w, y, m


def test_direct_latent_state_mixture_and_score_decomposition():
    theta, x, w, y, m = example()
    t = history_terms(theta, x, w, y, m)
    psi = expit(x @ theta[:2]); p = expit(w @ theta[2:])
    likelihood = []
    for i in range(len(x)):
        pattern = np.prod([p[i, j] if y[i, j] else 1-p[i, j] for j in np.flatnonzero(m[i])])
        absent = 1-psi[i] if not y[i, m[i]].any() else 0.
        likelihood.append(absent + psi[i] * pattern)
    np.testing.assert_allclose(np.exp(t['log_history']), likelihood, rtol=1e-13)
    np.testing.assert_allclose(t['log_history'], t['log_any'] + t['log_pattern_given_any'], atol=1e-13)


def test_all_binary_histories_sum_to_one():
    theta, x, w, _, m = example()
    for i in range(len(x)):
        values = []
        for pattern in product([0., 1.], repeat=int(m[i].sum())):
            y = np.zeros((1, 3)); y[0, m[i]] = pattern
            values.append(np.exp(history_terms(theta, x[i:i+1], w[i:i+1], y, m[i:i+1])['log_history'][0]))
        assert sum(values) == pytest.approx(1., abs=1e-13)


@pytest.mark.parametrize('precision', [0., 1.])
def test_analytic_gradient_matches_finite_difference(precision):
    theta, x, w, y, m = example()
    _, gradient = objective(theta, x, w, y, m, precision)
    numeric = []
    for j in range(len(theta)):
        d = np.zeros(len(theta)); d[j] = 1e-6
        numeric.append((objective(theta+d, x, w, y, m, precision)[0] -
                        objective(theta-d, x, w, y, m, precision)[0]) / 2e-6)
    np.testing.assert_allclose(gradient, numeric, atol=2e-9)


def test_untrapped_source_zeros_contribute_no_likelihood():
    theta, x, w, y, m = example()
    changed = y.copy(); changed[~m] = np.nan
    a = objective(theta, x, w, y, m); b = objective(theta, x, w, changed, m)
    np.testing.assert_allclose(a[0], b[0], atol=0.)
    np.testing.assert_allclose(a[1], b[1], atol=0.)
    m[0] = False
    with pytest.raises(ValueError, match='at least one surveyed year'):
        history_terms(theta, x, w, y, m)


def test_habitat_terms_use_only_training_cells():
    values = np.array([[.2, 2.], [.5, 3.], [.9, 8.], [100., 1000.]])
    transform = HabitatTransform.fit(values[:3], 1)
    train = transform.apply(values[:3])
    np.testing.assert_allclose(train[:, 1:].mean(axis=0), 0., atol=1e-14)
    np.testing.assert_allclose(train[:, 1:].std(axis=0), 1., atol=1e-14)
    before = transform.apply(values[:3])
    values[-1] *= 1e6
    np.testing.assert_array_equal(transform.apply(values[:3]), before)


def test_observation_scaling_excludes_held_cells_and_untrapped_years():
    values = np.array([[[1., 3.], [3., 7.]], [[np.nan, np.nan], [5., 11.]], [[999., 999.], [999., 999.]]])
    observed = np.array([[1, 1], [0, 1], [1, 1]], bool)
    w, mean, scale = observation_design(values, observed, np.array([0, 1]))
    np.testing.assert_allclose(mean, [3., 7.])
    np.testing.assert_allclose(w[:2][observed[:2]][:, 1:].mean(axis=0), 0., atol=1e-14)
    np.testing.assert_array_equal(w[1, 0], [1., 0., 0.])
    assert np.isfinite(scale).all()
