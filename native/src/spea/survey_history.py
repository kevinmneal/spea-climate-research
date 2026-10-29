"""Source-preserving annual USGS histories and exact within-cell conditional scores.

These histories concern encounters during surveyed years. Conditioning on each
cell's total detections eliminates its intercept, but supplies no occupancy
probability and cannot identify time-varying availability versus detectability.
"""
from __future__ import annotations
from dataclasses import dataclass
from itertools import combinations
from math import comb, lgamma

import numpy as np
import pandas as pd
from scipy.special import logsumexp


def annual_usgs_history(survey: pd.DataFrame, effort: pd.DataFrame,
                        years=range(1995, 2017)) -> pd.DataFrame:
    """Join explicit source cell IDs and retain unobserved years as missing outcomes."""
    if not survey.cellID.is_unique or not effort.cellID.is_unique:
        raise ValueError('Survey cell identifiers must be unique in each source table')
    if survey.cellID.isna().any() or effort.cellID.isna().any():
        raise ValueError('Missing survey cell identifiers')
    if set(survey.cellID) != set(effort.cellID):
        raise ValueError('Survey and effort identifiers must agree exactly')
    aligned = effort.set_index('cellID').loc[survey.cellID]
    rows = []
    for year in years:
        suffix = f'{year % 100:02}'
        y = survey[f'SPHA.{suffix}'].to_numpy(float)
        days = aligned[f'eff.{suffix}'].to_numpy(float)
        rain = aligned[f'wy{suffix}'].to_numpy(float)
        if not np.isin(y, [0, 1]).all():
            raise ValueError('Annual source outcomes must be explicit binary values')
        if np.isinf(days).any() or np.any(days < 0):
            raise ValueError('Invalid trapping effort')
        observed = np.isfinite(days) & (days > 0)
        if np.any((y == 1) & ~observed):
            raise ValueError('Detection without positive documented effort')
        if not np.isfinite(rain[observed]).all() or np.any(rain[observed] < 0):
            raise ValueError('Observed years need valid reported precipitation')
        # This source's FGDC dictionary explicitly defines eff.* NA as no trapping.
        # Do not generalize that meaning to an unrelated source's missing effort.
        state = np.where(observed, np.where(y == 1, 'surveyed_detection', 'surveyed_nondetection'),
                         np.where(np.isnan(days), 'not_trapped_per_source_metadata', 'zero_effort'))
        rows.append(pd.DataFrame(dict(
            source_id='S35', source_cell_id=survey.cellID.to_numpy(), year=year,
            source_survey_row=np.arange(len(survey)) + 2,
            published_block=survey.block.to_numpy(), published_random_fold=survey.fold.to_numpy(),
            source_pitfall_arrays=survey.n.to_numpy(), source_n_det_unresolved=survey['n.det'].to_numpy(),
            source_any_method_detection=survey.SPHA.to_numpy(),
            source_annual_outcome=y.astype('int8'), trapping_days=days,
            water_year_precipitation_mm=rain, observed=observed,
            encounter_outcome=np.where(observed, y, np.nan), observation_state=state,
        )))
    result = pd.concat(rows, ignore_index=True).sort_values(['source_cell_id', 'year']).reset_index(drop=True)
    if result.duplicated(['source_cell_id', 'year']).any():
        raise ValueError('Duplicate annual observation key')
    return result


@dataclass
class ConditionalHistory:
    """Exact sums over annual binary patterns with the observed total held fixed."""
    source_cell_id: int
    n: int
    k: int
    pattern_sums: np.ndarray
    observed_sum: np.ndarray

    @classmethod
    def from_arrays(cls, source_cell_id, x, y, maximum_patterns=1_000_000):
        x = np.asarray(x, dtype=float); y = np.asarray(y)
        if x.ndim != 2 or y.shape != (len(x),) or not np.isfinite(x).all() or not np.isin(y, [0, 1]).all():
            raise ValueError('Finite predictors and aligned binary outcomes required')
        n = len(y); k = int(y.sum())
        if not 0 < k < n:
            raise ValueError('Conditional histories require both encounter states')
        count = comb(n, k)
        if count > maximum_patterns:
            raise ValueError('Exact enumeration exceeds declared pattern cap')
        patterns = np.asarray(list(combinations(range(n), k)), dtype='int64')
        sums = x[patterns].sum(axis=1)
        return cls(int(source_cell_id), n, k, sums, x[y == 1].sum(axis=0))

    def nll_gradient_hessian(self, beta):
        beta = np.asarray(beta, dtype=float)
        eta = self.pattern_sums @ beta
        normalizer = logsumexp(eta)
        weights = np.exp(eta-normalizer)
        mean = weights @ self.pattern_sums
        centered = self.pattern_sums-mean
        hessian = centered.T @ (weights[:, None]*centered)
        return (float(normalizer-self.observed_sum @ beta), mean-self.observed_sum, hessian)

    @property
    def uniform_nll(self):
        return lgamma(self.n+1)-lgamma(self.k+1)-lgamma(self.n-self.k+1)


def objective(histories, beta):
    dimension = len(beta)
    value = 0.; gradient = np.zeros(dimension); hessian = np.zeros((dimension, dimension))
    for history in histories:
        v, g, h = history.nll_gradient_hessian(beta)
        value += v; gradient += g; hessian += h
    return value, gradient, hessian
