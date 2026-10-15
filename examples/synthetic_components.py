"""Synthetic climate and distance examples; no records, fits or ecological inference."""
import json
import numpy as np
from spea_snapshot.seasonal_model import monthly_features
from spea_snapshot.climate_statistics import statistical_features
from spea_snapshot.spatial_weights import proximity_weights, normalize_presence_weights
from spea_snapshot.spatial_thinning import thin_attempts

temperature = np.array([[8, 9, 12, 16, 20, 24, 27, 26, 22, 17, 12, 9]] * 2, dtype=float)
winter_rain = np.array([90, 75, 60, 25, 10, 2, 1, 2, 8, 30, 60, 85], dtype=float)
precipitation = np.stack([winter_rain, np.roll(winter_rain, 6)])
monthly = monthly_features(precipitation, temperature)
statistics, audit = statistical_features(precipitation, temperature)
distances = np.abs(np.array([0., 1., 10.])[:, None] - np.array([0., 1., 10.])[None, :])
masks, chosen = thin_attempts(distances, 2., range(10))
print(json.dumps(dict(
    example="SYNTHETIC: two monthly profiles and three points on an abstract kilometre axis",
    monthly_columns=len(monthly.columns), added_summary_columns=len(statistics.columns),
    equal_annual_controls=bool(np.array_equal(monthly.iloc[0, :2], monthly.iloc[1, :2])),
    precipitation_temperature_offset_c=statistics["precipitation_temperature_offset_c"].tolist(),
    normalized_soft_weights=normalize_presence_weights(proximity_weights(distances, 10.)).tolist(),
    selected_thinning_mask=masks[chosen].tolist(), new_fits=0), indent=2))
