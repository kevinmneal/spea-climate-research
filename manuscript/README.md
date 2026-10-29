# Recorded climatic differences and geographic prediction across Spea

Kevin M. Neal · Working manuscript · 29 October 2026 · Public companion v0.3.0

Read the [PDF](Spea_Comparative_Manuscript.pdf), comment in the [Word document](Spea_Comparative_Manuscript.docx), or inspect the [complete Markdown source](spea_comparative_manuscript.md). The document includes the main text and references, three main figures, supporting methods and tables, and the coarse evaluation-geography figure.

## Scientific scope

The paper connects three completed comparisons: North/South recorded climates versus controlled response gains; matched annual, variability, monthly and full-BIOCLIM representations under equal and soft spatial weights; and the separate southern occurrence and survey-history water comparisons. It retains the full-BIOCLIM primary advantage, secondary-score disagreements, source and penalty sensitivity, sparse water-block influence and the adverse survey result. These developmental comparisons do not establish intrinsic population responses, new boundaries or population history.

## Evidence and reproduction

The [evidence index](evidence_index.md) connects each claim and figure to its complete saved comparisons. `source_data/` contains aggregate scores, summaries, original source identities and provenance. `figures/` contains the four displayed figures as PNG and SVG. Numeric source files and figures retain their original bytes.

From the repository root, the standard-library replay checks all 3,664 packaged arithmetic and identity comparisons:

```sh
python3 manuscript/replay/replay_aggregates.py --output replay_output/manuscript_receipt.json
```

With the optional dependencies in `requirements-components.txt`, rebuild the three main figures:

```sh
python3 manuscript/replay/build_main_figures.py --output-dir replay_output/manuscript_figures
```

With `python-docx` installed, rebuild the editable document without overwriting the delivered copy:

```sh
python3 manuscript/assembly/build_document.py
```

The builder writes into `manuscript/rebuilt/`. Rendering to PDF requires a compatible Word/LibreOffice renderer. Different renderers or fonts can change pagination. [The release verification](verification/public_release_validation.json) records what was actually checked. Earlier native audit receipts are dated historical evidence, not checks newly executed by public CI.

The [native code](../native/README.md), [reproducibility guide](../REPRODUCIBILITY.md) and [source/access inventory](../SOURCES_AND_ACCESS.md) describe the original fitting implementation and missing archived inputs. Aggregate replay does not refit models, rebuild Figure S1's native geography, or independently validate observations or ecological inference. Exact records, native masks/arrays, model artifacts and third-party source papers are not bundled here.

## Manuscript history

This integrated paper supersedes the focused seasonal manuscript as the active manuscript. Its numerical findings remain incorporated with their qualifications. The older text is preserved in the [v0.2.0 snapshot](https://github.com/kevinmneal/spea-climate-research/blob/v0.2.0/docs/manuscript.md) and the [earlier releases](../VERSION_HISTORY.md). The intermontana case remains background synthesis; this update does not introduce a second overlapping paper.

The release updates availability/status wording and documents an earlier supplementary-methods clarification: precipitation log transformations and full-native-target integration apply to specified comparisons, not every analysis. No new scientific fit, prediction or occurrence admission was performed for this release. The manuscript has not been submitted to bioRxiv or presented as peer reviewed.
