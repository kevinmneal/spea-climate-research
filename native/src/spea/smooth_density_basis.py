"""Training-availability basis for a bounded conditional-density comparison.

This supplies quadratic climate responses whose coefficients may vary smoothly
with signed geographic distance and/or differ between independently assigned
groups. It does not turn the geographic group contrast into a discontinuity
test. All learned transformations use only the rows passed to ``fit``.
"""
from __future__ import annotations

import numpy as np
import sklearn
from sklearn.preprocessing import SplineTransformer


def weighted_quantiles(values, weights, probabilities):
    """Linearly interpolate the empirical weighted CDF at distinct values.

    Duplicate values first receive their summed mass. CDF coordinates are the
    cumulative mass at each distinct value (not midpoint ranks). Quantiles
    below the first CDF coordinate equal the minimum; 0 and 1 return the range
    endpoints. This definition is deterministic under row permutation and
    splitting a weighted row into identical weighted pieces.
    """
    x, w, q = (np.asarray(v, dtype=float) for v in (values, weights, probabilities))
    if (x.ndim != 1 or not len(x) or w.shape != x.shape or q.ndim != 1
            or not all(np.isfinite(v).all() for v in (x, w, q))
            or np.any(w <= 0) or np.any((q < 0) | (q > 1))):
        raise ValueError('Quantiles need finite values, positive weights and probabilities in [0, 1]')
    unique, inverse = np.unique(x, return_inverse=True)
    mass = np.bincount(inverse, weights=w / w.max(), minlength=len(unique))
    cdf = np.cumsum(mass)
    cdf /= cdf[-1]
    return np.interp(q, cdf, unique, left=unique[0], right=unique[-1])


def _inputs(climate, distance_m, groups):
    x = np.asarray(climate, dtype=float)
    s = np.asarray(distance_m, dtype=float)
    g = np.asarray(groups)
    if (x.ndim != 2 or x.shape[1] != 4 or not len(x)
            or not np.isfinite(x).all()):
        raise ValueError('Climate must be a nonempty finite matrix with four columns')
    if s.shape != (len(x),) or not np.isfinite(s).all():
        raise ValueError('Distance must be an aligned finite vector in meters')
    if g.shape != (len(x),) or not np.isin(g, [-1, 1]).all():
        raise ValueError('Groups must be an aligned vector of -1 and/or +1')
    return x, s, g


def _equal_group_weights(weights, groups):
    w = np.asarray(weights, dtype=float)
    if (w.shape != groups.shape or not np.isfinite(w).all()
            or np.any(w <= 0) or set(np.unique(groups)) != {-1, 1}):
        raise ValueError('Fit requires both groups and aligned finite positive availability weights')
    result = np.empty_like(w)
    for group in (-1, 1):
        use = groups == group
        scaled = w[use] / w[use].max()
        result[use] = .5 * scaled / scaled.sum()
    return result


class SmoothDensityBasis:
    """Fit availability-weighted transformations; reuse them on held-out rows.

    Supply four already physically/log-transformed climate axes. ``weights``
    must be training quadrature availability, never occurrence weights. Each
    group supplies half of the transformation mass regardless of its area.
    ``transform`` returns the full CL design. Slice with ``columns[model]``
    for S (shared), L (group interaction), C (geographic interaction), or CL.

    Climate basis h = [z, z²], where z uses training weighted mean and SD.
    Geography is either one linear distance column, or five cubic B-spline
    columns with four training weighted-quantile knots and linear extrapolation.
    Geographic columns are centered before multiplication by h. The full raw
    order is [h, geography, geography*h, group*h], with geography-major
    interaction order. Every column is then centered within each training group
    and scaled by the equal-group weighted RMS. Group centering changes only
    conditional normalizer constants; scaling defines the coefficient penalty.
    """

    distance_scale_m = 800_000.
    schema_version = 1

    def _set_columns(self):
        ng = 5 if self.geography == 'spline' else 1
        base_end, smooth_end = 8 + ng, 8 + ng + 8 * ng
        self.columns = {
            'S': list(range(base_end)),
            'L': list(range(base_end)) + list(range(smooth_end, smooth_end + 8)),
            'C': list(range(smooth_end)),
            'CL': list(range(smooth_end + 8)),
        }
        h_names = [f'z{i}' for i in range(4)] + [f'z{i}^2' for i in range(4)]
        geo_names = [f'geo{i}' for i in range(ng)]
        self.column_names = (h_names + geo_names
                             + [f'{geo}*{h}' for geo in geo_names for h in h_names]
                             + [f'group*{h}' for h in h_names])

    def _make_spline(self):
        self._spline = None
        if self.geography == 'spline':
            if (self.knots_scaled_.shape != (4,)
                    or not np.isfinite(self.knots_scaled_).all()
                    or not np.all(np.diff(self.knots_scaled_) > 0)):
                raise ValueError('Spline requires four distinct increasing training knots')
            self._spline = SplineTransformer(
                n_knots=4, degree=3, knots=self.knots_scaled_[:, None],
                extrapolation='linear', include_bias=False,
            ).fit(self.knots_scaled_[:, None])

    def _geography(self, distance_m):
        scaled = distance_m / self.distance_scale_m
        return (self._spline.transform(scaled[:, None])
                if self.geography == 'spline' else scaled[:, None])

    def _raw(self, climate, distance_m, groups):
        z = (climate - self.climate_mean_) / self.climate_scale_
        h = np.column_stack([z, z * z])
        geo = self._geography(distance_m) - self.geography_mean_
        return np.column_stack([
            h, geo, (geo[:, :, None] * h[:, None, :]).reshape(len(h), -1),
            groups[:, None] * h,
        ])

    def fit(self, climate, distance_m, groups, weights, geography='spline'):
        """Learn from training availability rows only and return ``self``."""
        if geography not in ('spline', 'linear'):
            raise ValueError('Geography must be spline or linear')
        x, s, g = _inputs(climate, distance_m, groups)
        w = _equal_group_weights(weights, g)
        self.geography = geography
        self.sklearn_version_ = sklearn.__version__
        self.climate_mean_ = w @ x
        climate_sd = np.sqrt(w @ np.square(x - self.climate_mean_))
        self.constant_climate_ = climate_sd <= 1e-12
        self.climate_scale_ = np.where(self.constant_climate_, 1., climate_sd)
        self.knots_scaled_ = (weighted_quantiles(s / self.distance_scale_m, w, [0, 1/3, 2/3, 1])
                              if geography == 'spline' else np.empty(0))
        self._make_spline()
        self.geography_mean_ = w @ self._geography(s)
        self._set_columns()
        raw = self._raw(x, s, g)
        self.group_mean_ = np.array([
            np.average(raw[g == group], axis=0, weights=w[g == group])
            for group in (-1, 1)
        ])
        centered = raw - self.group_mean_[(g == 1).astype(int)]
        rms = np.sqrt(w @ np.square(centered))
        self.constant_columns_ = rms <= 1e-12
        self.column_scale_ = np.where(self.constant_columns_, 1., rms)
        return self

    def transform(self, climate, distance_m, groups):
        """Return all CL columns without learning anything from these rows."""
        if not hasattr(self, 'column_scale_'):
            raise ValueError('Fit the basis before transforming data')
        x, s, g = _inputs(climate, distance_m, groups)
        result = (self._raw(x, s, g) - self.group_mean_[(g == 1).astype(int)]) / self.column_scale_
        if not np.isfinite(result).all():
            raise ValueError('Basis overflowed on supplied prediction values')
        return result

    def metadata(self):
        """Return JSON-ready fitted parameters; no training rows or pickle."""
        if not hasattr(self, 'column_scale_'):
            raise ValueError('Fit the basis before exporting metadata')
        return {
            'schema_version': self.schema_version, 'geography': self.geography,
            'sklearn_version': self.sklearn_version_,
            'distance_scale_m': self.distance_scale_m,
            'knot_quantile_method': 'distinct_value_weighted_empirical_cdf_linear_interpolation',
            'spline_degree': 3, 'spline_include_bias': False, 'spline_extrapolation': 'linear',
            'group_order': [-1, 1], 'climate_mean': self.climate_mean_.tolist(),
            'climate_scale': self.climate_scale_.tolist(),
            'constant_climate': self.constant_climate_.tolist(),
            'knots_scaled': self.knots_scaled_.tolist(),
            'geography_mean': self.geography_mean_.tolist(),
            'group_mean': self.group_mean_.tolist(),
            'column_scale': self.column_scale_.tolist(),
            'constant_columns': self.constant_columns_.tolist(),
            'columns': {k: v.copy() for k, v in self.columns.items()},
            'column_names': self.column_names.copy(),
        }

    @classmethod
    def from_metadata(cls, metadata):
        """Restore fixed transformations from their JSON representation."""
        if (metadata.get('schema_version') != cls.schema_version
                or metadata.get('geography') not in ('linear', 'spline')
                or metadata.get('distance_scale_m') != cls.distance_scale_m
                or metadata.get('group_order') != [-1, 1]
                or metadata.get('spline_degree') != 3
                or metadata.get('spline_include_bias') is not False
                or metadata.get('spline_extrapolation') != 'linear'):
            raise ValueError('Unsupported basis metadata parameters')
        obj = cls()
        obj.geography = metadata['geography']
        obj.sklearn_version_ = metadata['sklearn_version']
        obj._set_columns()
        size, ng = len(obj.columns['CL']), 5 if obj.geography == 'spline' else 1
        shapes = {
            'climate_mean': (4,), 'climate_scale': (4,), 'constant_climate': (4,),
            'knots_scaled': (4,) if obj.geography == 'spline' else (0,),
            'geography_mean': (ng,), 'group_mean': (2, size),
            'column_scale': (size,), 'constant_columns': (size,),
        }
        for name, shape in shapes.items():
            value = np.asarray(metadata[name], dtype=bool if name.startswith('constant_') else float)
            if value.shape != shape or not np.isfinite(value).all():
                raise ValueError(f'Invalid saved {name}')
            setattr(obj, name + '_', value)
        if np.any(obj.climate_scale_ <= 0) or np.any(obj.column_scale_ <= 0):
            raise ValueError('Saved scales must be positive')
        if metadata.get('columns') != obj.columns or metadata.get('column_names') != obj.column_names:
            raise ValueError('Saved column definitions do not match basis schema')
        obj._make_spline()
        return obj
