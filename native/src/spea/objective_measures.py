"""Explicit weight encodings for an optional presence/background ablation.

These helpers do not change existing fits, select a landscape, or assign
occupancy labels. Callers remain responsible for the empirical and reference
measures, inclusion probabilities, strata, and library penalty conventions.
"""
from __future__ import annotations

import numpy as np


def _mass(values, name):
    mass = np.array(values, dtype=np.float64, copy=True)
    if (mass.ndim != 1 or not len(mass) or not np.isfinite(mass).all()
            or np.any(mass < 0) or not np.isfinite(mass.sum()) or mass.sum() <= 0):
        raise ValueError(f'{name} must be finite nonnegative masses with a positive total')
    return mass


def binary_class_weights(presence_mass, integration_mass, *, tau=1., loss_mass=1.):
    """Weights for concatenated presence rows then reference-landscape rows.

    Each input is normalized separately. Total positive weight is loss_mass;
    total background weight is tau * loss_mass, independent of row counts.
    tau expresses the source-class loss ratio, never biological prevalence.
    The same cell can legitimately occur in both classes. This helper uses a
    global mixture; separately conditioned regional objectives are different.
    """
    p = _mass(presence_mass, 'presence_mass')
    q = _mass(integration_mass, 'integration_mass')
    if not np.isfinite(tau) or tau <= 0 or not np.isfinite(loss_mass) or loss_mass <= 0:
        raise ValueError('tau and loss_mass must be finite and positive')
    with np.errstate(over='raise', invalid='raise'):
        weights = np.r_[loss_mass * (p / p.sum()),
                        (loss_mass * tau) * (q / q.sum())]
    if not np.isfinite(weights).all() or not np.isfinite(weights.sum()):
        raise ValueError('Binary loss weights overflow')
    return weights


def poisson_quadrature_encoding(empirical_mass, integration_mass):
    """Return target=e/w and weight=w on an aligned integration support.

    No implicit normalization or epsilon area is introduced. Empirical mass
    can be noninteger and is an optimization measure, not animal counts.
    Every node must have positive integration mass. Zero-area presence-only
    rows from independent quadrature need a different encoding. Callers must
    construct event-inclusive quadrature without double-counting area.
    """
    e = _mass(empirical_mass, 'empirical_mass')
    w = _mass(integration_mass, 'integration_mass')
    if e.shape != w.shape or np.any(w <= 0):
        raise ValueError('Aligned masses and positive integration mass at every node required')
    with np.errstate(over='ignore', divide='ignore', invalid='ignore'):
        target = e / w
    if not np.isfinite(target).all():
        raise ValueError('Poisson optimization target is not representable')
    return target, w
