# Reproducing the reported results

Version 0.3.0 provides the integrated manuscript's saved numerical comparisons, portable arithmetic and main-figure builders, and the original native implementation with a source/access inventory. The immediately runnable layer checks saved results. Exact occurrence-to-fit reproduction still needs archived inputs that are not all public or redistributable.

## What is supplied

| Level | What runs from this repository | What it establishes |
| --- | --- | --- |
| Integrated manuscript arithmetic | Standard-library replay of 3,664 checks | Saved structured-response, pooling, matched seasonal, procedure-history, water and support identities |
| Three main figures | NumPy/Matplotlib builder from packaged tables | Displays the saved comparisons; adds no fit or uncertainty estimate |
| Editable manuscript | python-docx builder from Markdown and figures | Recreates document content and layout instructions; PDF rendering is separate |
| Older seasonal and environmental arithmetic | Existing standard-library replay and tests | Preserved v0.1.0/v0.2.0 table and screen calculations |
| Selected numerical components | Analytic/synthetic tests and example | Implemented transformations and sampling rules on controlled inputs |
| Native fitters and settings | Supplied for inspection; syntax/import/help checks performed | Original implementation and project dependencies; not a successful native refit |
| Figure S1 native geography or full model replay | Requires omitted native assets | Not reproduced by the portable checks or public CI |

The original and matched comparisons reuse observations, fits and geography. Their row counts must not be added as independent experiments. Passing a hash or arithmetic check does not validate the observations, model assumptions, population interpretation or novelty of the paper.

## Integrated manuscript replay

With Python 3.13, from the repository root:

```sh
python3 manuscript/replay/replay_aggregates.py --output replay_output/manuscript_receipt.json
```

No third-party package or network access is required. The default without `--output` writes under `manuscript/rebuilt/`, preserving the released receipt. The replay verifies 42 packaged source-file hashes, structured means/contrasts, all matched climate-family and weight-arm scores, pooling diagnostics, both water comparisons, original procedure history, saved alignment counts and 3,914 displayed held cells. Floating-point comparisons use absolute and relative tolerance 1e-12. [The evidence index](manuscript/evidence_index.md) maps claims to the actual tables and specifies each score's support.

Scores use natural logarithms, whose units are **nats**. The descriptive score name and averaging unit matter more than that label: relative recording-density scores and log probabilities per complete survey-cell history are different quantities and are never pooled into a single effect. Fold spread and positive-fold counts are descriptive, not confidence intervals or significance tests.

The earlier 2,510-check native-source audit is retained as a dated receipt. The public replay does not rerun it. Likewise, saved model/array hashes identify native assets; they do not make those assets available or verify their contents anew.

## Rebuild figures and the editable document

Use an isolated Python 3.13 environment for optional figure/component dependencies:

```sh
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements-components.txt
.venv/bin/python manuscript/replay/build_main_figures.py --output-dir replay_output/manuscript_figures
```

The builder reconstructs Figures 1–3 from packaged aggregate scores. It preserves the different water estimands, adverse results, all retained folds and North/South labels. It adds no fits, predictions, confidence intervals or new observations. Figure S1 is supplied as an inspected image; rebuilding it needs omitted native masks, support/split arrays and coordinate-bearing display inputs.

To rebuild the Word document:

```sh
.venv/bin/python -m pip install -r requirements-document.txt
.venv/bin/python manuscript/assembly/build_document.py
```

The output is `manuscript/rebuilt/Spea_Comparative_Manuscript.docx`. Transfer author edits from Word back into the Markdown source before rebuilding. Word, LibreOffice or another compatible renderer is required to make a PDF; typography and pagination can vary with fonts and renderer. The delivered PDF was rendered and visually checked. [Public release validation](manuscript/verification/public_release_validation.json) states the actual document and image checks.

## Preserve the earlier replays

The former seasonal and environmental-organization modules remain runnable:

```sh
python3 scripts/replay.py --output replay_output/legacy
python3 -m unittest discover -s tests -p 'test_replay*.py'
```

The script refuses an existing destination. It reconstructs the original selected/fallback and fixed-baseline seasonal tables, matched weighting/BIOCLIM tables, 897 environmental-organization statistics, nine simulation envelopes and two association summaries. Seasonal averages give folds equal weight; environmental statistics use their specified representative or block weighting. Its empirical simulation counts are not calibrated p-values, and it does not rerun the screen's native models or simulations.

The optional selected-component checks are unchanged:

```sh
PYTHONPATH=src .venv/bin/python -m pytest tests/components -q
PYTHONPATH=src .venv/bin/python examples/synthetic_components.py
.venv/bin/python scripts/plot_results.py --output replay_figure.png
```

The example uses synthetic inputs and is not ecological validation. [Component documentation](src/README.md) records the original source, small namespace/constant adaptations and the omitted full-model integration test.

## Inspect the native implementation

Start with [native/README.md](native/README.md), which gives strand-specific entry points, an inspection environment and conditional native commands. [code_selection.json](native/code_selection.json) records 51 original files, their roles/hashes, project import closure, historical identity comparisons and 455 omitted direct native/audit dependencies. It also identifies upstream preparation that remains outside this selection.

The original `native/pyproject.toml` and `native/uv.lock` are inert identity evidence used by historical producer guards. Do not install the broad original project merely to inspect these scripts; the native guide provides a smaller separate inspection environment. Native commands require restored parent inputs, source snapshots and audit artifacts, with original checks left intact. A successful `--help`, syntax check or synthetic test is not a successful fit. Generated native inputs/outputs can contain exact records and locations and are ignored by the public repository.

## Source access and exact restoration

[SOURCES_AND_ACCESS.md](SOURCES_AND_ACCESS.md) and its [machine-readable inventory](metadata/public_sources.json) record the GBIF API query pages and provider datasets, published and author-source additions, all 43 CHELSA objects, CEC/Natural Earth definitions, and Rose/USGS table/code identities. Follow those provider routes for public materials. They do not retroactively create an atomic GBIF snapshot or recover unpublished accepted-record decisions.

Exact restoration additionally requires the retained source bytes and adjudications, accepted-cell membership, physical arrays, available/held masks, guards and quadrature, training-only transformations, model artifacts and source-ID survey joins. Some inputs are request-dependent, subject to source rights and location sensitivity; this release does not promise their availability. Fresh downloads or a changed environment need their own provenance and validation. No survey coordinates are reconstructed from the released background grid.

## Release integrity

```sh
python3 scripts/check_snapshot.py
```

`CHECKSUMS.sha256` covers released files and the checker validates their bytes and local Markdown links. `provenance.json` records the source checkpoint, exact copies and documented export adaptations. Public CI runs the aggregate replays, selected-component tests, figure builders and native syntax checks. It does not download occurrences, refit native models or establish independent ecological replication.

Cite a versioned release or commit. [Version history](VERSION_HISTORY.md) links the earlier manuscript, while [CHANGELOG.md](CHANGELOG.md) describes each update. Existing release tags are preserved.
