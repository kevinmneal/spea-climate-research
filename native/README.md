# Native analysis code: scope and reproduction

This folder exposes the original fitting, transformation and saved-fit synthesis code behind the integrated *Spea* manuscript. It contains a small, project-import-complete selection: 37 Python files, ten frozen configurations, two small audit/source catalogs and the original whole-project environment identity baselines. Files retain their original paths and bytes under this folder. The public manuscript's portable aggregate replay remains the immediately executable route to its reported arithmetic and main figures.

**Native fitting does not run from a clean public clone.** The frozen record-level inputs, source adjudications, native grids and masks, spatial roles, parent-run manifests and other audit assets are not distributed here. This code lets readers inspect the implementation; an exact native rerun additionally requires the original authorized archive. Current provider downloads do not reconstruct that archive by themselves. No full clone-to-manuscript guarantee is made.

The accompanying [code inventory](code_selection.json) records every selected path, role, byte count and SHA-256; transitive project imports; the original producer run identities; and 455 omitted direct native input/audit dependencies found in saved receipts and pinned configurations. It distinguishes those direct dependencies from upstream raw-data ancestry, which was not exhaustively exported. Current selected code/configuration bytes matched all 101 applicable saved-source identity comparisons. Those checks establish identity, not biological validity.

## The manuscript's execution strands

| Strand | Actual entry points | What is inherited or qualified |
| --- | --- | --- |
| North/South controlled responses | `scripts/run_hammondii_working_comparison.py`, `config/hammondii_working_comparison.json` | Requires the saved full native reference geography, support refresh, source-membership audit and accepted-cell/BIOCLIM preparation. Fits S/L/C/CL and declared sensitivities. |
| Pooling and reciprocal transfer | `scripts/run_hammondii_reciprocal_transfer.py`, `config/hammondii_reciprocal_transfer.json` | Reuses the controlled comparison's saved splits and analysis grid, plus the registered climate/terrain arrays. `run_genus_environment_comparison.py` is included because it supplies imported helpers; its standalone broader experiment is not another manuscript result. |
| Matched seasonal representations and weights | `prepare_climate_statistics_sampling.py` -> `run_climate_statistics_comparison.py`; `prepare_soft_spatial_weights.py` -> `run_soft_spatial_weighting.py`; then `run_seasonal_weighting_interaction.py` and `summarize_seasonal_bioclim_reference.py` (all in `scripts/`) | The final 168 slots combine saved fits from several runs. The parent climate-statistics run also includes historical summary/thinning arms outside the main manuscript comparison. Its full execution is not a minimal new experiment. The matched run fits annual/variability arms and verifies/reuses monthly fits. The final synthesis reads scores and small model metadata, adding no fit or prediction. |
| Southern occurrence-water comparison | `scripts/run_south_water_benchmark.py`, `config/south_water_benchmark_v1.json` | Requires original public USGS presence/background points, metadata, fixed land-masked physical BIOCLIM arrays and pinned source-audit artifacts. It uses raw physical BIOCLIM inputs with training-background standardization, and log-transformed water distances. |
| Southern survey-history comparison | `scripts/run_south_survey_water.py`, `config/south_survey_water_v1.json` | Requires source-keyed survey/effort tables and pinned source-code/audit receipts. Survey coordinates are neither reconstructed nor needed by this diagnostic. It does not reproduce the released integrated USGS model. |

The original seasonal selected/fallback and deeper fixed-baseline procedures are retained in `run_genus_seasonal_comparison.py` and `run_hammondii_species_seasonal_comparison.py`, with their configurations. They explain the supporting procedure-history results; they do not replace the later matched fixed-shallow main comparison. The older pooled-hammondii slots in the genus procedure are not substituted for the paper's separate North/South slots.

## Inspection and conditional commands

From this folder (`native/` in the public companion), use Python 3.13 and preserve `scripts/`, `src/spea/` and `config/` together. Project imports use:

```sh
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"
```

The selected programs import NumPy, pandas, SciPy, scikit-learn, CatBoost, affine, pyproj, requests, urllib3 and, for the saved-fit figure/synthesis script, Matplotlib. Several fitters explicitly require CatBoost 1.2.10. Original settings/receipts record the run-specific environment. The byte-identical `pyproject.toml` and `uv.lock` are supplied solely as historical source-identity baselines, because original producer guards reference them. They include unrelated research dependencies: **do not install the whole project or synchronize that broad lock just to inspect this selection**. Use the separate `requirements-native-inspection.txt` for a prospective inspection environment; it omits TabM and Torch. The inspection requirements pin the versions used for this release's import/help and synthetic-component checks. They are not a recovered per-fit historical lock and do not certify native model reproduction. These scripts use Unix `resource` accounting and were developed on macOS; resource guard units need review before native refitting on another platform.

A prospective inspection environment can be prepared separately from the historical baseline:

```sh
python3 -m venv .venv-inspection
.venv-inspection/bin/python -m pip install -r requirements-native-inspection.txt
```

Use that environment's Python for the examples below. The release checks used an existing environment with the listed versions; a fresh installation was not performed locally. With numerical dependencies installed, `--help` describes an entry point without fitting or downloading, for example:

```sh
python3 scripts/run_hammondii_working_comparison.py --help
python3 scripts/run_south_survey_water.py --help
```

The following are **conditional native refit commands**, not commands that work from the distributed files alone. They require all original parent assets, relative paths and identity checks to be restored in a separate authorized workspace. Use fresh output IDs; the programs intentionally refuse existing output directories.

```sh
python3 scripts/run_hammondii_working_comparison.py --run-id local_controlled_refit
python3 scripts/run_hammondii_reciprocal_transfer.py --run-id local_transfer_refit
python3 scripts/run_seasonal_weighting_interaction.py --run-id local_matched_seasonal_refit
python3 scripts/run_south_water_benchmark.py --config config/south_water_benchmark_v1.json --run-id local_occurrence_water_refit
python3 scripts/run_south_survey_water.py --config config/south_survey_water_v1.json --run-id local_survey_water_refit
```

These commands use the original archived parents named in their configurations. The new `local_*` outputs are not automatically wired into later stages. For example, reciprocal transfer still expects the original controlled-reference run, and the seasonal synthesis still expects the original three saved-fit source runs. Replacing those parents is a separately versioned reconstruction requiring updated provenance and checks; do not bypass failed hash checks.

The historical final seasonal synthesis command is:

```sh
python3 scripts/summarize_seasonal_bioclim_reference.py --config config/seasonal_bioclim_reference_synthesis_v1.json
```

It requires restored parent scores, manifests, model metadata and cutoff identity receipts. Its configuration fixes output locations and it refuses existing destinations. Preserve the original outputs; choose distinct destinations in a separately saved configuration if reconstructing it in the full archive. This saved-fit synthesis is distinct from the public standard-library aggregate replay, which requires none of these native assets.

## Omitted inputs and access

The inventory's omitted paths are archive-relative identities, not public download URLs or a promise that a new request will return the same bytes. Native registries and producer manifests define physical arrays, masks, quadrature inclusion probabilities, source support, held geography and fitted-output identities. Small audit support files are supplied: the six original test modules referenced by producers, the separate North/South split preparer, two source/review catalogs, and original environment baseline files. These are identified separately from the 455 omitted native data/run assets in the inventory. Inspected third-party source-code snapshots and native audit receipts remain omitted. Original guards are unchanged; absence or identity mismatch in a required archived artifact will stop execution even if the numerical arrays have been restored. Historical review-gate statuses describe their original workflow dates and are not a current publication-authorization statement.

Public primary-source entry points include CHELSA 2.1 normals ([EnviDat](https://doi.org/10.16904/envidat.228)), the USGS occurrence/survey release ([data](https://doi.org/10.5066/P9O7MC7M)) and the official USGS integrated-model code ([software](https://doi.org/10.5066/P9YYKW1H)). These identify sources, not an automated reconstruction pipeline. The GBIF portion was acquired through cached live API pages, with per-page query/date/hash receipts and per-record provider terms; no GBIF download DOI was obtained. Fresh API results are not the frozen extraction.

The accepted occurrence preparation also contains published and author-source contributions. Those source-linked records and adjudications are not supplied here. Access depends on each original source's terms and permissions; this release does not promise their availability or authorize redistribution. First identify the exact omitted asset and its source identity in [Sources and access](../SOURCES_AND_ACCESS.md) and the [code inventory](code_selection.json). Obtain public materials from their providers; discuss access to unavailable derived/frozen inputs with the author without posting sensitive localities publicly. No public-issue request should include exact sensitive locations.

Upstream preprocessing is deliberately distinguished from fitting. The inventory names the original accepted-cell preparation, spatial design, native reference geography, source-support refresh, no-S35 membership audit and separate North/South split preparer. The North/South split preparer is included as a small source-snapshot dependency; the other upstream preparers are indexed rather than bundled. Their saved outputs are prerequisites here; rerunning those upstream procedures on new sources constitutes a new preparation, not byte-identical restoration. Their full source/acquisition dependency chains are outside this small code selection.

## Output and licensing boundaries

Running the native programs can create exact coordinates, native cell identifiers, source membership, full survey histories, arrays and model snapshots. Keep generated `runs/` and record-level inputs private unless each artifact has been separately reviewed for release. This code release does not reconstruct obscured survey coordinates and does not authorize doing so.

The code selection includes original Python implementations, not downloaded USGS R source, third-party palette source or manuscripts. Preserve the companion's original-software license, NOTICE and method attributions, including enmSdmX distance weighting and the spThin context. Provider records, climate products, maps and publications retain their own terms. Aggregate arithmetic checks, source hashes and AI-assisted code review do not constitute independent ecological replication or peer review.
