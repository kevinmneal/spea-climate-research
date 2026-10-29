"""Convex finite-basis and exact Bernoulli weighting reference calculations.

These helpers are mathematical controls for synthetic design checks. They do
not establish a biological response, a calibrated null test, or a new ecological
method. Availability and empirical mass are normalized separately within each
observed group; groups contribute equally to the conditional objective.
"""
from __future__ import annotations

from functools import lru_cache

import numpy as np
from scipy.special import roots_legendre


def _vector(value, size, name):
    result = np.asarray(value, dtype=float)
    if result.shape != (size,) or not np.isfinite(result).all():
        raise ValueError(f'{name} must be an aligned finite vector')
    return result


def _groups(value, size):
    result = np.asarray(value)
    if result.shape != (size,) or not np.isin(result, [-1, 1]).all():
        raise ValueError('Groups must be an aligned vector of -1 and/or +1')
    return result


def _normalized_positive(values):
    """Normalize nonnegative mass without an overflowing unscaled sum."""
    maximum = values.max()
    if maximum <= 0:
        raise ValueError('Every group needs positive empirical mass')
    scaled = values / maximum
    return scaled / scaled.sum()


def _log_normalized_area(values):
    # Keep arbitrarily small positive areas in log space: a finite score can
    # counteract an extreme area ratio even when that ratio itself underflows.
    log_area = np.log(values)
    shifted = log_area - log_area.max()
    return shifted - np.log(np.exp(shifted).sum())


def _prepare(design, groups, area, empirical_mass, ridge):
    x = np.asarray(design, dtype=float)
    if x.ndim != 2 or min(x.shape) < 1 or not np.isfinite(x).all():
        raise ValueError('Design must be a finite nonempty matrix')
    g = _groups(groups, len(x))
    a = _vector(area, len(x), 'Area')
    mass = _vector(empirical_mass, len(x), 'Empirical mass')
    if np.any(a <= 0) or np.any(mass < 0):
        raise ValueError('Area must be positive and empirical mass nonnegative')
    if not np.isfinite(ridge) or ridge <= 0:
        raise ValueError('A finite positive coefficient ridge is required')
    strata = []
    for group in np.unique(g):
        use = g == group
        strata.append((x[use], _log_normalized_area(a[use]), _normalized_positive(mass[use])))
    return x.shape[1], strata


def _statistics(beta, strata, ridge):
    value = .5 * ridge * np.dot(beta, beta)
    gradient = ridge * beta.copy()
    hessian = ridge * np.eye(len(beta))
    group_weight = 1. / len(strata)
    for x, log_area, mass in strata:
        eta = x @ beta
        # A group-specific shift cancels from the conditional log likelihood.
        # Subtracting it here avoids loss of precision for redundant intercepts.
        centered = eta - eta.max()
        weighted = log_area + centered
        shift = weighted.max()
        exponential = np.exp(weighted - shift)
        distribution = exponential / exponential.sum()
        log_normalizer = shift + np.log(exponential.sum())
        value += group_weight * (log_normalizer - np.dot(mass, centered))
        mean = distribution @ x
        gradient += group_weight * (mean - mass @ x)
        centered_x = x - mean
        hessian += group_weight * (centered_x.T @ (distribution[:, None] * centered_x))
    return float(value), gradient, hessian


def value_gradient_hessian(beta, design, groups, area, empirical_mass, ridge):
    """Return penalized conditional loss, exact gradient and full Hessian.

    For each group the loss is ``log E_area[exp(X beta)] - E_mass[X beta]``.
    The objective averages these group losses and adds ``ridge * ||beta||²/2``.
    Its Hessian is the group-averaged model-weighted covariance of the design
    plus the ridge identity, including all cross-coefficient terms.
    """
    n_features, strata = _prepare(design, groups, area, empirical_mass, ridge)
    coefficient = _vector(beta, n_features, 'Coefficients')
    return _statistics(coefficient, strata, ridge)


def fit_conditional_basis(design, groups, area, empirical_mass,
                          ridge=1e-5, tol=1e-8):
    """Fit the strictly convex finite-basis reference with damped exact Newton.

    Rank-redundant columns, intercepts and group-constant columns are allowed;
    the positive ridge supplies a unique minimum. No feature scaling is done
    internally, because that would change the declared coefficient penalty.
    An unconverged fit raises rather than silently returning usable parameters.
    """
    n_features, strata = _prepare(design, groups, area, empirical_mass, ridge)
    if not np.isfinite(tol) or not 0 < tol <= 1e-6:
        raise ValueError('Tolerance must be positive and no larger than 1e-6')
    beta = np.zeros(n_features)
    value, gradient, hessian = _statistics(beta, strata, ridge)
    for iteration in range(201):
        gradient_inf = float(np.max(np.abs(gradient)))
        if gradient_inf <= tol:
            return dict(coefficients=beta, objective=value,
                        gradient_inf=gradient_inf, iterations=iteration,
                        success=True)
        if iteration == 200:
            break
        direction = -np.linalg.solve(hessian, gradient)
        directional_derivative = float(gradient @ direction)
        if not np.isfinite(direction).all() or directional_derivative >= 0:
            raise RuntimeError('Finite-basis Newton direction failed')
        step = 1.
        for _ in range(60):
            candidate = beta + step * direction
            candidate_value, candidate_gradient, candidate_hessian = _statistics(candidate, strata, ridge)
            # Very near a minimum, the decrease can be below a floating-point
            # objective ULP. Accept only a lower gradient in that tiny window.
            rounding_window = 8 * np.finfo(float).eps * max(1., abs(value))
            sufficient_decrease = candidate_value <= value + 1e-4 * step * directional_derivative
            rounding_progress = (abs(candidate_value - value) <= rounding_window
                                 and np.max(np.abs(candidate_gradient)) < gradient_inf)
            if np.isfinite(candidate_value) and (sufficient_decrease or rounding_progress):
                beta = candidate
                value, gradient, hessian = candidate_value, candidate_gradient, candidate_hessian
                break
            step *= .5
        else:
            raise RuntimeError(f'Finite-basis line search did not converge; gradient={gradient_inf:.3g}')
    raise RuntimeError(f'Finite-basis fit did not converge; gradient={gradient_inf:.3g}')


def relative_log_density(eta, groups, area):
    """Group-normalized log density relative to uniform density per unit area.

    Within each group, ``sum(area * exp(result)) / sum(area) == 1``.
    Group-specific additive score shifts and changes of area units cancel.
    """
    score = np.asarray(eta, dtype=float)
    if score.ndim != 1 or not len(score) or not np.isfinite(score).all():
        raise ValueError('Scores must be a nonempty finite vector')
    g = _groups(groups, len(score))
    a = _vector(area, len(score), 'Area')
    if np.any(a <= 0):
        raise ValueError('Area must be positive')
    result = np.empty_like(score)
    for group in np.unique(g):
        use = g == group
        centered = score[use] - score[use].max()
        weighted = _log_normalized_area(a[use]) + centered
        shift = weighted.max()
        result[use] = centered - shift - np.log(np.exp(weighted - shift).sum())
    return result


def _probabilities(probability):
    p = np.asarray(probability, dtype=float)
    if (p.ndim != 1 or not len(p) or not np.isfinite(p).all()
            or np.any(p < 0) or np.any(p > 1)):
        raise ValueError('Probabilities must be a nonempty finite vector in [0, 1]')
    if not np.any(p > 0):
        raise ValueError('Conditional mass is undefined when all probabilities are zero')
    return p


@lru_cache(maxsize=32)
def _unit_quadrature(order):
    nodes, weights = roots_legendre(order)
    return (nodes + 1.) / 2., weights / 2.


def _inverse_count_integrals(probability):
    """For each i return integral prod_{j != i}(1-p_j+p_j*t) dt exactly.

    A polynomial of degree n-1 needs ceil(n/2) Gauss–Legendre nodes. Integrands
    are evaluated in bounded chunks using log products. Nodes are interior,
    so p=1 has no zero denominator. Zero-probability items may be included.
    """
    p = np.asarray(probability, dtype=float)
    nodes, weights = _unit_quadrature((len(p) + 1) // 2)
    result = np.zeros(len(p))
    # Bound temporary matrices to approximately 4 MiB each for large vectors.
    chunk = min(64, max(1, 524288 // len(p)))
    for start in range(0, len(nodes), chunk):
        t = nodes[start:start + chunk]
        log_factors = np.log1p(-p[None, :] * (1. - t[:, None]))
        excluded_products = np.exp(log_factors.sum(axis=1)[:, None] - log_factors)
        result += weights[start:start + chunk] @ excluded_products
    return result


def expected_balanced_cell_mass(probability, blocks):
    """Exact expected block-balanced cell weights conditional on any presence.

    Cells are independent Bernoulli events. For a realized nonempty field,
    each occupied block gets equal mass, split equally among its occupied
    cells. This function returns E[I_i / (K * N_block(i)) | K > 0].

    Within-block and between-block reciprocal counts are integrated separately
    by polynomial-exact Gauss–Legendre quadrature. This is an expectation of
    normalized random weights, generally different from normalized expected
    counts. Zero-probability cells are omitted algebraically, never pruned at a
    nonzero threshold. The cost is quadratic in the largest nonzero cell block
    and in the number of potentially occupied blocks; use a bounded grid.
    """
    p = _probabilities(probability)
    b = np.asarray(blocks)
    if b.shape != p.shape:
        raise ValueError('Blocks must be aligned with probabilities')
    # Reject missing labels; np.unique otherwise silently treats these as a
    # block despite NaN not comparing equal to itself in subsequent masks.
    if (b.dtype.kind in 'fc' and not np.isfinite(b).all()) or (b.dtype.kind == 'O' and any(v is None or v != v for v in b)):
        raise ValueError('Blocks cannot contain missing values')
    _, inverse = np.unique(b, return_inverse=True)
    active_blocks = np.unique(inverse[p > 0])
    q = np.empty(len(active_blocks))
    within = np.zeros(len(p))
    members = []
    for j, block in enumerate(active_blocks):
        idx = np.flatnonzero((inverse == block) & (p > 0))
        members.append(idx)
        with np.errstate(divide='ignore'):
            q[j] = -np.expm1(np.log1p(-p[idx]).sum())
        within[idx] = p[idx] * _inverse_count_integrals(p[idx])
    between = _inverse_count_integrals(q)
    with np.errstate(divide='ignore'):
        conditional_probability = -np.expm1(np.log1p(-p).sum())
    result = np.zeros(len(p))
    for j, idx in enumerate(members):
        # Divide first so extremely rare fields do not underflow in a product.
        result[idx] = (within[idx] / conditional_probability) * between[j]
    if not np.isfinite(result).all() or not np.isclose(result.sum(), 1., rtol=1e-9, atol=1e-12):
        raise RuntimeError('Bernoulli expected mass failed its normalization check')
    return result


def expected_cell_mass(probability):
    """E[I_i / total occupied cells | any occupied cell], without block weights."""
    p = _probabilities(probability)
    return expected_balanced_cell_mass(p, np.zeros(len(p), dtype=int))
