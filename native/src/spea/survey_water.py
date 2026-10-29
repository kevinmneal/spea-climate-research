"""Static survey-history mixture for a bounded habitat comparison.

The latent state is a modeling assumption over 1995--2016, not independently
verified occupancy. Annual encounters combine availability and detection.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import expit


@dataclass
class HabitatTransform:
    """Training-cell scales; quadratic terms also receive unit training variance."""

    mean: np.ndarray
    scale: np.ndarray
    square_mean: np.ndarray
    square_scale: np.ndarray
    terrestrial_count: int

    @classmethod
    def fit(cls, values, terrestrial_count):
        values = np.asarray(values, float)
        mean = values.mean(axis=0)
        scale = values.std(axis=0)
        if not np.isfinite(values).all() or np.any(scale <= 0):
            raise ValueError("Finite, varying training habitat required")
        squares = ((values[:, :terrestrial_count] - mean[:terrestrial_count]) /
                   scale[:terrestrial_count]) ** 2
        square_mean = squares.mean(axis=0)
        square_scale = squares.std(axis=0)
        if np.any(square_scale <= 0):
            raise ValueError("Quadratic training terms must vary")
        return cls(mean, scale, square_mean, square_scale, terrestrial_count)

    def apply(self, values):
        z = (np.asarray(values, float) - self.mean) / self.scale
        squares = (z[:, :self.terrestrial_count] ** 2 - self.square_mean) / self.square_scale
        return np.column_stack((np.ones(len(z)), z, squares))

    def as_dict(self):
        return {k: (v.tolist() if isinstance(v, np.ndarray) else v)
                for k, v in vars(self).items()}


def observation_design(values, observed, training_cells):
    """Scale effort/rain using only surveyed years in training cells."""
    values = np.asarray(values, float)
    train = values[training_cells][observed[training_cells]]
    mean = train.mean(axis=0)
    scale = train.std(axis=0)
    if not np.isfinite(train).all() or np.any(scale <= 0):
        raise ValueError("Finite, varying training observation predictors required")
    z = np.where(observed[:, :, None], (values - mean) / scale, 0.)
    return np.concatenate((np.ones((*observed.shape, 1)), z), axis=2), mean, scale


def history_terms(theta, habitat, observation, y, observed):
    """Marginalize the two latent states, retaining missing years as unobserved."""
    if np.any(observed.sum(axis=1) == 0):
        raise ValueError("Each scored history needs at least one surveyed year")
    if not np.isin(y[observed], [0, 1]).all():
        raise ValueError("Surveyed outcomes must be binary")
    dimension = habitat.shape[1]
    eta = habitat @ theta[:dimension]
    tau = observation @ theta[dimension:]
    log_psi = -np.logaddexp(0., -eta)
    log_not_psi = -np.logaddexp(0., eta)
    log_p = -np.logaddexp(0., -tau)
    log_not_p = -np.logaddexp(0., tau)
    yy = np.where(observed, y, 0.)
    log_pattern = np.where(observed, yy * log_p + (1. - yy) * log_not_p, 0.).sum(axis=1)
    log_zero = np.where(observed, log_not_p, 0.).sum(axis=1)
    positive = yy.sum(axis=1) > 0
    log_history = np.where(positive, log_psi + log_pattern,
                           np.logaddexp(log_not_psi, log_psi + log_zero))
    posterior_available = np.ones(len(eta))
    zero = ~positive
    posterior_available[zero] = np.exp(log_psi[zero] + log_zero[zero] - log_history[zero])
    # Probability of at least one encounter under the actual supplied exposure.
    log_q = log_psi + np.log(-np.expm1(log_zero))
    log_not_q = np.logaddexp(log_not_psi, log_psi + log_zero)
    log_any = np.where(positive, log_q, log_not_q)
    return dict(log_history=log_history, log_any=log_any,
                log_pattern_given_any=np.where(positive, log_history - log_any, 0.),
                psi=expit(eta), p=expit(tau), q=np.exp(log_q),
                posterior_available=posterior_available, positive=positive)


def objective(theta, habitat, observation, y, observed, ridge_precision=1.):
    """Sum of cell-history negative log likelihoods and fixed slope penalty."""
    terms = history_terms(theta, habitat, observation, y, observed)
    n_habitat = habitat.shape[1]
    r = terms['posterior_available']
    habitat_gradient = habitat.T @ (r - terms['psi'])
    residual = np.where(observed, y - terms['p'], 0.)
    observation_gradient = np.einsum('i,ij,ijk->k', r, residual, observation)
    penalty = np.ones(len(theta)) * ridge_precision
    penalty[[0, n_habitat]] = 0.  # The two intercepts are not penalized.
    value = -terms['log_history'].sum() + .5 * np.dot(penalty, theta ** 2)
    gradient = -np.r_[habitat_gradient, observation_gradient] + penalty * theta
    return float(value), gradient


def numerical_hessian(theta, function):
    """Diagnostic curvature only; no independent-year confidence intervals."""
    columns = []
    for j in range(len(theta)):
        delta = np.zeros(len(theta)); delta[j] = 1e-5 * max(1., abs(theta[j]))
        columns.append((function(theta + delta)[1] - function(theta - delta)[1]) / (2 * delta[j]))
    h = np.column_stack(columns)
    return (h + h.T) / 2
