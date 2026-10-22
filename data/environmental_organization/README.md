# Environmental-organization aggregates

Saved results from the completed 20 October 2026 *S. intermontana* screen. Read the [scientific note](../../docs/environmental_organization.md) and [frozen design adaptation](../../metadata/environmental_organization_design.json). These files contain no occurrence coordinates, native-cell indices or block identifiers.

| File | Rows or contents | Meaning |
|---|---|---|
| [mixture_fold_scores.csv](mixture_fold_scores.csv) | 4,485 rows | Five fold contributions for each of 897 observed/simulated representative datasets and role layouts |
| [organization_statistics.csv](organization_statistics.csv) | 897 rows | Sum of `gain_sum` divided by sum of `held_blocks`; includes six observed and 891 simulated/layout combinations |
| [envelope_summary.csv](envelope_summary.csv) | Nine rows | Three generators × three selections; 5th/50th/95th percentiles, counts at or above the observation, and upper-envelope flag |
| [association_fold_scores.csv](association_fold_scores.csv) | Ten rows | Equal-block mean conditional score for each of five folds and two association models; held cells and blocks retained |
| [association_summary.csv](association_summary.csv) | Two rows | Fold means weighted by held-block count, covering the same 134 blocks and 358 cells for each model |
| [feature_variance.csv](feature_variance.csv) | Six rows | Saved fraction of available-climate variance within blocks, plus the percentage representation |
| [cohort_summary.json](cohort_summary.json) | Aggregate counts | Original, in-mask and precise cohorts; outside-mask and precision exclusions; native availability counts |

`selection=0` is primary; `1` and `2` are prespecified choices of a representative recorded cell in each block. `generator=observed` has `replicate=-1`. Simulated `replicate` values run from 0 to 98 within `uniform`, `shared` and `smooth`; each simulated dataset is evaluated under all three role layouts. There are 297 simulated datasets, not 891 independent biological samples. The two observed cohorts are `primary` and `broader_centers`. The broader cohort has no simulated controls.

`gain_mean` is the two-component minus one-component held mean log density, and `gain_sum` is that difference summed over held representatives. `train_blocks` and `held_blocks` retain the original representative eligibility for the fold. `smallest_component_mass` is the smaller fitted mixture weight times its training count; `small_component` flags values below five. Iteration counts are recorded diagnostics, not a new convergence proof in the public replay.

Association `mean_log_density` is relative to physical-area availability within each block. The public exporter averages the native saved block scores by fold/model and omits their block IDs. Reconstructing the all-held-block score requires weighting the five exported means by their `held_blocks`, rather than averaging folds equally. That differs from the seasonal manuscript's equal-fold summary target.

Simulation percentiles use linear interpolation at `(n - 1) * q`, matching the saved numerical analysis. Exceedance counts include equality. These envelopes condition on fitted parameters and do not yield calibrated p-values. Feature-variance percentages are checked against their saved fractions; native climate variances are not recomputed here.

```sh
python3 scripts/replay_environmental_organization.py --output environmental_replay_output
```

Run that command from the repository root. Outputs include observed and conditional-association Markdown tables, reconstructed statistic/envelope/association CSVs and a receipt with source hashes and explicit scope. The main [table replay](../../REPRODUCIBILITY.md) runs this module automatically. The existing seasonal data and the new screen are related analyses and must not be pooled as independent replicates.
