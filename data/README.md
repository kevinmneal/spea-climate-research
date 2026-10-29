# Aggregate results

The integrated manuscript's complete three-strand evidence is in [manuscript/source_data](../manuscript/evidence_index.md). This directory preserves the earlier seasonal, sensitivity, decomposition and separate environmental-organization modules, with their existing replay commands.

All supplied rows summarize recorded cells or fitted models. They contain no individual occurrence coordinates or reconstructable native-cell lists. Seasonal CSV source bytes are preserved, including historical column names and run identifiers. The new environmental-organization export preserves existing aggregate rows and documents its additional fold aggregation. Public document wording does not rename historical keys.

## Seasonal units and metrics

| Key | Meaning |
|---|---|
| `intermontana`, `bombifrons`, `multiplicata` | Nominal *Spea* fitting units; five held folds each |
| `hammondii_north`, `spea_south` | Northern and southern *S. hammondii* operational units; three held folds each |
| `mean_block_log_score` | Primary log relative recording density, equally averaged across occupied 50-km blocks within each held fold |
| `mean_cell_log_score` | Secondary log relative recording density, equally averaged across held recorded cells within a fold |
| `area_weighted_source_auc` | Secondary ROC AUC for recorded-source versus available-area classification, with block-weighted records and area-weighted background |
| `source_average_precision_prior_half` | Secondary source average precision at an artificial 1:1 recorded/area source mass |

Log scores are in nats, with zero representing uniform density per physical area within the held target. Larger is better for this prediction task. They are not occupancy probabilities. Means across folds give each fold equal weight; different units have different geographical tasks, so absolute means are not a common species-specialization scale. See [Methods](../docs/methods.md#score-and-diagnostic-weights).

## Original comparison

[`original/`](original) contains five files consumed by the original portable replay: `primary_folds.csv`, `primary_summary.csv`, `primary_contrasts.csv`, `procedure_contrasts.csv`, and `support_by_fold.csv`.

The 126 slots cross three families (annual 2 columns, variability 4, monthly 28), five units, and two procedures. `fixed` is the original 300-tree depth-four baseline. `selected` uses eligible inner-fold selection, with the declared baseline fallback in two southern outer folds. These are distinct from the later fixed shallow comparison. Some slots are exact model aliases, not independent fits.

## Matched fixed shallow comparison

[`matched/`](matched) contains `fold_metrics.csv`, `family_summary.csv`, `contrast_folds.csv`, and `contrast_summary.csv`. Every one of 168 family/arm/fold slots is supplied, including adverse and mixed results. The main replay reconstructs all four tables' core scores from `fold_metrics.csv`.

| Dimension | Values |
|---|---|
| Family | `annual` (2), `variability` (4), `monthly` (28), `bioclim` (19) |
| Training arm | `equal_cells`; `soft_geo_10km` (all cells retained with distance weights) |
| Learner | Fixed shallow CatBoost, 150 trees / depth 2 / learning rate 0.03 |
| Primary contrast | Monthly minus variability; monthly minus full BIOCLIM is the non-nested reference |

`features` is the number of predictor columns. `train_cells`, `held_cells`, `held_blocks`, and `target_cells` describe support. `stem`, `model_source_run`, `metrics_source_run` and `reused` preserve source lineage identifiers; they are not paths to public model assets. Other fields are saved diagnostics not rederived by this aggregate replay. In particular, recomputing `target_log_normalizer`, probability-mass checks, or training-marginal support requires excluded native arrays.

All contrast names specify left minus right. Summary `_min` and `_max` are fold extrema; `_positive_folds` counts strictly positive differences. None is a confidence interval or a significance test. Full BIOCLIM uses partly different thermal information, so its reference contrast is not a nested ablation of monthly inputs.

## Additional summaries and sampling sensitivities

[`sensitivities/`](sensitivities) retains the supplied feature, hard-thinning and soft-weight results, with `south_procedure.csv` separating a learner-setting change from a data or sampling gain. These are fixed-setting sensitivities. The main replay does not reconstruct their original models or every summary. No optimal thinning radius or inferred sampling-effort surface is selected by this snapshot.

## Saved-score decomposition

[`decomposition/`](decomposition) contains aggregate fold, stratum and mean components for intermontana's matched monthly-minus-variability contrast. The two partitions describe cooler/warmer precipitation association and inside/outside the monthly training-quadrature marginal ranges. Allocation and within-stratum contributions sum to total score change under the original evaluation weights.

These tables support the manuscript's accounting discussion; no individual cell or held-block coordinates are included. The main replay does not regenerate the partition from native predictions. A marginal-range flag does not establish joint environmental support, and an accounting identity is not a causal explanation of seasonal biology.

## Environmental organization

[`environmental_organization/`](environmental_organization/README.md) contains the completed screen's 4,485 mixture fold rows, 897 organization statistics, nine control envelopes, ten association fold aggregates, two association summaries, six saved variance fractions and cohort counts. The main replay reconstructs the statistic/envelope/association arithmetic. Its gains are aggregated across held representatives or blocks, rather than giving folds equal weight. Read the module's definitions before comparing these scores with the seasonal tables.

All three precise-cohort observed selections, all three broader-center sensitivities and every simulated/layout result are retained, including the sensitivity that exceeds both fitted upper envelopes. No association block identifiers or native climate-cell identities are exported. The six variance fractions cannot be recomputed from these aggregates alone.
