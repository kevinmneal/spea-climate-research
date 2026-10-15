"""Nested monthly-climate representations for the seasonal information test."""
from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd

ANNUAL = ['annual_mean_temperature_c', 'log1p_annual_precipitation_mm']
VARIABILITY = ANNUAL + ['monthly_temperature_range_c', 'monthly_precipitation_cv']
MONTHLY = VARIABILITY + [f'temperature_anomaly_{m:02}' for m in range(1, 13)] + [
    f'precipitation_fraction_{m:02}' for m in range(1, 13)]
FAMILIES = {'annual': ANNUAL, 'variability': VARIABILITY, 'monthly': MONTHLY}


def monthly_features(precipitation, temperature):
    """Physical normals -> exact nested summaries, returned as float32.

    All 12 months have equal weight in temperature summaries and population
    precipitation CV. Precipitation is total precipitation, including snow.
    Monthly fractions and temperature anomalies retain calendar timing while
    separating it algebraically from annual amount/warmth. There is no fitting
    or scaling against held or training distributions in these definitions.
    """
    p, t = np.asarray(precipitation, dtype=float), np.asarray(temperature, dtype=float)
    if p.ndim != 2 or p.shape[1] != 12 or p.shape != t.shape:
        raise ValueError('Aligned rows of 12 monthly precipitation and temperature values required')
    if not np.isfinite(p).all() or not np.isfinite(t).all() or np.any(p < 0):
        raise ValueError('Finite temperatures and nonnegative precipitation required')
    annual = p.sum(axis=1)
    if np.any(annual <= 0):
        raise ValueError('Positive annual precipitation required; do not drop cells')
    mean_t = t.mean(axis=1)
    frame = np.column_stack([mean_t, np.log1p(annual), np.ptp(t, axis=1),
                             np.std(p, axis=1, ddof=0)/(annual/12),
                             t-mean_t[:, None], p/annual[:, None]])
    return pd.DataFrame(frame.astype('float32'), columns=MONTHLY)


class MonthlyReader:
    def __init__(self, root, registry):
        self.root = Path(root)
        self.shape = tuple(registry['window_grid']['shape'])
        fields = {f['key']: f for f in registry['fields']}
        self.paths = {v: [self.root/fields[f'{v}_{m:02}']['physical_path'] for m in range(1,13)]
                      for v in ['pr','tas']}

    def read(self, native):
        native = np.asarray(native, dtype=np.int64)
        if native.ndim != 1 or np.any(native < 0) or np.any(native >= np.prod(self.shape)):
            raise ValueError('Valid flat native-grid cell IDs required')
        values = {}
        for variable in ['pr','tas']:
            frame = np.empty((len(native),12), dtype='float32')
            for month,path in enumerate(self.paths[variable]):
                source = np.load(path, mmap_mode='r')
                if source.shape != self.shape:
                    raise ValueError('Monthly source grid changed')
                frame[:, month] = source.ravel()[native]
                source._mmap.close()
            values[variable] = frame
        return monthly_features(values['pr'], values['tas'])


def prefix_borders(source, destination, feature_count):
    """Reuse identical training-only cutoffs for shared nested columns."""
    source, destination = Path(source), Path(destination)
    if destination.exists():
        raise FileExistsError(destination)
    lines = [line for line in source.read_text().splitlines()
             if line.strip() and int(line.split('\t')[0]) < feature_count]
    destination.write_text('\n'.join(lines)+'\n')
