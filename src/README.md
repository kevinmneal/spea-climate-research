# Selected scientific components

These are small numerical components from the working research code, distributed to make the implemented methods inspectable and reusable. They are not the complete fitting pipeline or a packaged general-purpose species-distribution-modeling library.

| Module | Implemented work |
|---|---|
| `spea_snapshot/seasonal_model.py` | Nested 2/4/28-column representations; prefix-preserving training cutoffs; aligned-array reader |
| `spea_snapshot/climate_statistics.py` | Fourteen additional distributional and circular seasonal summaries; exact feature-family definitions; shared-cutoff remapping |
| `spea_snapshot/spatial_thinning.py` | WGS84 distances; repeated maximum-conflict-degree removal; first largest retained set, without claiming global optimality |
| `spea_snapshot/spatial_weights.py` | Soft distance weights, fixed-total normalization and a weight-concentration diagnostic |

The original source files and their SHA-256 values are recorded in [provenance.json](../provenance.json). Three modules are exact copies. `climate_statistics.py` uses the identical `LOG_KEYS` constant inline instead of importing unrelated model/raster code. The module namespace changes to `spea_snapshot`; the formulas do not change.

The scientific tests retain the original analytic checks, including circular month windows, skew against SciPy, cutoff identity, distance behavior at high latitude/across the date line, graph separation, and weight normalization. One full-model objective integration test is excluded with its absent dependency. Passing these tests does not validate ecological assumptions or reproduce fitted predictions.

See [optional commands](../REPRODUCIBILITY.md#optional-code-and-figure). Array-reader classes need user-supplied aligned arrays and registries; this repository ships neither. The [synthetic example](../examples/synthetic_components.py) needs no biological or geographic inputs and fits no model.
