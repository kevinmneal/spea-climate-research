# Reproducibility of this snapshot

Version 0.2.0 exposes the arithmetic behind the central seasonal comparisons and the separate environmental-organization screen, plus selected reusable methods. It is not a complete occurrence-to-model release. The numerical results were produced in the larger working archive; this snapshot does not imply that those inputs can be reconstructed exactly from today's provider downloads.

## Included levels

| Level | Available here? | What the check establishes |
|---|---|---|
| Original table arithmetic | Yes | Reconstructs 30 family means and 126 paired fold contrasts from 126 selected/fixed family/fold slots; verifies 21 fold-support rows and 3,914 held cells |
| Matched weighting and BIOCLIM arithmetic | Yes | Reconstructs 40 family means, 168 paired fold contrasts and 40 contrast summaries from 168 family/arm/fold slots; checks all four reported metrics |
| Environmental-organization arithmetic | Yes | Reconstructs 897 statistics from 4,485 fold rows, nine simulation envelopes and two conditional association aggregates from ten fold rows; checks cohort/diagnostic consistency and six saved variance-percentage conversions |
| Reusable predictor and sampling components | Yes, optional dependencies | Tests numerical transformations and geometric sampling rules on analytic or synthetic inputs |
| Public overview figure | Yes, optional Matplotlib | Replots the primary matched comparison from supplied fold scores |
| Other manuscript figures and sensitivity/decomposition tables | Saved outputs | Supplied for inspection; the main replay does not regenerate these figures or their underlying native predictions |
| Occurrence preparation, native predictions and model fitting | No | Frozen record-level inputs, accepted cell assignments, native arrays, models and other assets are absent |
| Historical routes or genetic subdivisions | No | The agenda is exploratory context; this replay tests no historical or genetic claim |
| Seasonal analogue pilot | Working-archive summary only | Brief context in the agenda; pair matching and its native inputs are not replayed here |

The original and matched sets reuse models and observations, and the environmental screen reuses recorded cells and developmental geography. Their slot counts must not be added as independent model experiments or biological replicates. The screen has 297 simulated datasets evaluated under three role layouts, not 891 independent biological samples. Neither a replay pass nor an input checksum validates the scientific assumptions or supplies fresh evaluation geography.

## Main replay

From the repository root, with Python 3.13:

```sh
python3 scripts/replay.py --output replay_output
python3 -m unittest discover -s tests -p 'test_replay*.py'
```

There is no dependency installation for these commands. The script refuses an existing output directory; use a new directory for another run. It checks saved arithmetic at absolute tolerance 1e-10, and rejects missing or duplicate folds and nonfinite or inconsistent scores. The tests also deliberately remove, duplicate and alter input rows to check that those errors are caught.

Outputs include:

- `original/table2.md` and `original/table3.md`: original selected/fallback and deeper fixed-baseline tables.
- `matched_primary_table.md`: current matched shallow comparison against variability and full BIOCLIM.
- `matched_family_means.csv`, `matched_fold_contrasts.csv`, `matched_contrast_summary.csv`: all four metrics, with every paired fold, means, extrema and positive-fold counts.
- `environmental_organization/observed_table.md` and `association_table.md`: all six observed representative outcomes and the separate conditional association comparison.
- `environmental_organization/statistics.csv`, `envelopes.csv` and `association_summary.csv`: reconstructed screen aggregates, including the isolated sensitivity discrepancy.
- `receipt.json`, `original/replay_receipt.json` and `environmental_organization/receipt.json`: runtime, input hashes and explicit scope of each replay.

The portable original script is an exact copy from the research archive. The later aggregate replay is new release code that reconstructs the already saved comparisons. It performs no fitting, predictor selection or new scientific experiment. Legacy original outputs use “Southern species”; this is the same `spea_south` operational unit called southern hammondii in the public text.

## Environmental-organization replay boundary

The [separate module](scripts/replay_environmental_organization.py) can also run alone:

```sh
python3 scripts/replay_environmental_organization.py --output environmental_replay_output
```

The organization statistic weights fold gains by their held representative counts. The association aggregate weights exported fold means by their held-block counts. Both differ from the seasonal manuscript's equal-fold averages. Linear simulation percentiles and empirical counts are recalculated from all saved statistics; those counts are not calibrated p-values.

The [design](metadata/environmental_organization_design.json) and [data definitions](data/environmental_organization/README.md) specify the fixed procedure. The public replay checks saved arithmetic and internal consistency. It does not recreate representative choices, native normalizers, model coefficients, GMM fits, simulated native cells or climate-variance estimates. The six variance fractions are supplied values whose percentage conversion is checked. The working archive's 1,405-check native verification is separately recorded in provenance and is not reproduced by this public CI.

## Optional code and figure

In a separate Python 3.13 environment:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-components.txt
PYTHONPATH=src .venv/bin/python -m pytest tests/components -q
PYTHONPATH=src .venv/bin/python examples/synthetic_components.py
.venv/bin/python scripts/plot_results.py --output replay_figure.png
```

The example is explicitly synthetic and is not a proxy for model evaluation. The selected components retain the original numerical functions; the only numerical-module adaptation replaces an imported constant with the identical literal set to avoid unrelated dependencies. Import changes and the omitted full-model integration test are recorded in [provenance.json](provenance.json).

The requirements file pins the installed versions used for selected-component verification. [The historical software snapshot](metadata/software_snapshot.json) separately records saved lock versions for the original fits; those locks do not establish a recovered per-fit installed-environment inventory. The environmental design record includes the installed versions captured by its own analysis receipt; its public aggregate replay needs no numerical-library dependencies.

## What full reproduction would require

An exact native replay would require the frozen source extracts and record adjudications; source uncertainty and event relationships; accepted native-cell assignments; aligned CHELSA 2.1 arrays and physical cell areas; complete available and held masks; guard decisions and source exclusions; sampled quadrature with inclusion probabilities; training-only cutoffs; fitted models and raw margins; and run-specific software/settings. Most of those assets are not distributed in this small companion.

[Methods](docs/methods.md), [mask definitions](metadata/mask_specification.json), [original method settings](metadata/method_specification.json), [references](metadata/references.bib) and fold summaries document the task. They cannot substitute for the absent frozen assets. Acquiring fresh data would constitute a new preparation requiring its own provenance and validation.

## Snapshot integrity and updates

`CHECKSUMS.sha256` lists the released files. `python3 scripts/check_snapshot.py` checks their bytes. `provenance.json` connects exported artifacts to the private source checkpoint without exposing its Git history. CI verifies the public replay and selected components, not the full research pipeline.

Release tags preserve dated snapshots. The main branch can evolve; cite a release or commit for comparisons. Scientific changes should update the manuscript, affected tables and scope statement together, with a correction or new analysis described in the change log.
