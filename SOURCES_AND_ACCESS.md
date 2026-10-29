# Sources and access

This companion makes coordinate-free saved results, their arithmetic and the three main manuscript figures reproducible. It does not distribute the complete occurrence-to-model research archive. The machine-readable [source inventory](metadata/public_sources.json) records provider identities, exact queries/object URLs, saved file hashes, source versions, access evidence and the role of each input class. Its full hashes identify the inspected receipts; hashes of native data are inherited from those receipts, not a new raw-data audit.

This document was prepared from saved source receipts on 29 October 2026. It reports the documented acquisition state. No underlying datasets were reacquired for this release. A working URL is a retrieval route, not a guarantee that today's payload matches the original bytes.

## Reproduction levels and material limits

| Level | Access and scope |
| --- | --- |
| Saved aggregate arithmetic and three main figures | The companion supplies compact scores, contrasts, source lineage and replay code. No record coordinates or native models are needed for this level. |
| Exact native occurrence-to-fit replay | Requires omitted frozen response pages, source adjudications, accepted cell assignments, native arrays/masks, quadrature, fitted models and run-specific settings. Some inputs are nonpublic or have unresolved redistribution rights. Access would require a separately reviewed arrangement; availability is not promised. |
| Fresh provider reacquisition | The recorded APIs, file IDs and object URLs permit a new preparation. Provider updates, live pagination and unpublished curation prevent treating it automatically as an exact recreation. |

Three limits are consequential. The GBIF base is a live API extraction, not a DOI-bearing atomic download. Accepted author-local additions and curation cannot be regenerated from public aggregate counts. The water analyses use distances supplied in the USGS tables; the complete versioned pond/stream geometry inventory was not acquired. The public background grid also supplies no documented survey-cell coordinate crosswalk.

## Occurrence inputs and curation

### GBIF base

The saved acquisition used [GBIF's occurrence API v1](https://techdocs.gbif.org/en/openapi/v1/occurrence). The 292 authoritative response pages were retrieved on **20 September 2026, 21:20:20–21:46:10 UTC**; extraction receipts completed at approximately 21:47 UTC. They represent global queries before project eligibility, deduplication and cell aggregation:

| Queried name | Recorded taxonKey | Returned unique GBIF IDs | Pagination |
| --- | ---: | ---: | --- |
| *Spea hammondii* | 2429783 | 6,435 | Species query exhausted |
| *S. intermontana* | 2429791 | 5,120 | Species query exhausted |
| *S. bombifrons* | 2429789 | 13,769 | 93 dataset partitions exhausted |
| *S. multiplicata* | 2429785 | 14,839 | 98 dataset partitions exhausted |

The base request was `https://api.gbif.org/v1/occurrence/search?taxonKey=<key>&limit=300&offset=<offset>`. For the two partitioned species, the saved facet query used `limit=0&facet=datasetKey&facetLimit=10000`; each returned `datasetKey` was then paginated separately. The inventory preserves every authoritative page URL, request parameters, timestamp and SHA-256, the partition counts, and the saved metadata identities of 139 provider datasets. Earlier abandoned broad-query pages are not additional observations.

**There is no GBIF download DOI for this extract.** Provider dataset DOIs recorded in the inventory identify those provider datasets, not a collective frozen download. The original response bytes remain in the research archive. Repeating the API calls today will query a changed live service; a new [authenticated GBIF predicate download](https://techdocs.gbif.org/en/data-use/api-downloads) could supply a new DOI-bearing snapshot, not retroactively identify this one. No account or download request was created during preparation of this companion.

GBIF record/provider licenses vary. The inventory reports the saved provider licenses individually; the research archive retains record-level attribution and terms. Neither the companion's license nor public API accessibility confers a blanket right to redistribute every underlying row. Exact locations and protected/obscured records are not exported or reconstructed here.

### Published and recovered public tables

The occurrence workflow also retained source-linked additions and cross-provider comparisons from:

| Source | Exact saved identity and retrieval route | Access/license evidence |
| --- | --- | --- |
| Neal et al. 2018 specimen supplement, S01 | [Article DOI](https://doi.org/10.1007/s10592-018-1066-7), file `10592_2018_1066_MOESM4_ESM.xlsx`; original workbook SHA-256 in the inventory | Public publisher supplement acquired; supplement redistribution license not stated in the saved receipt. |
| Legacy occurrence CSV, S03 | [Archived repository at commit `1f59a0693d9290325af175d1b7eeb5cc6ecd28c7`](https://github.com/kevinmneal-archive/spadefoot_machinelearning_R/tree/1f59a0693d9290325af175d1b7eeb5cc6ecd28c7), `speahammondii_combined_presence_points.csv` | Public bytes acquired; CSV license not stated in the saved receipt. This is an inherited source, not independent validation. |
| Dryad occurrence deposit, S05 | [doi:10.5061/dryad.f4qrfj6td](https://doi.org/10.5061/dryad.f4qrfj6td), version 2 / API version ID 76921, 11 August 2020; occurrence archive file ID 388844 | Saved metadata records CC0-1.0. The original archive/member names are `occurance/occuance_dataset.csv`, preserving the deposit's spelling. |
| Neal et al. Figshare coordinate tables, S30 | [doi:10.6084/m9.figshare.12721019.v1](https://doi.org/10.6084/m9.figshare.12721019.v1); `OCall208_85clust_popcoords_20190205.csv`, file 24079856, and `OCreduced148_85clust_popcoords_20190205.csv`, file 24079862 | Exact file URLs, SHA-256 and provider MD5 recorded; saved file metadata states CC BY 4.0. The tables overlap and must not be counted as independent animals. |

The exact normalizations and identity/source links matter. Published points are not automatically precise, dated or independent of GBIF entries. Conflicting assertions, unknown datum/uncertainty, date roles, rejected records and source versions remain distinct in the frozen preparation. The public source inventory identifies original bytes without distributing coordinate rows.

### Author-local additions and accepted decisions

The frozen genus preparation continued through dated source-review updates. Version 6 joined the audited public source base with recovered published assertions and curated alias/datum/record decisions. Version 7 added **95 author-source assertions and 26 existing-record enrichments**, rather than 121 new independent localities. Version 8 added **seven report-derived assertions representing four additional geographic supports**. These counts describe source processing; they are not the final fitting-cohort sizes.

The underlying author-local collection tables, the 2010 report/map georeference casework, exact coordinates, accepted record decisions and final source-to-cell assignments are retained in the nonpublic research archive. No verified general public download route or blanket redistribution permission is established for that complete curated set. Exact access is **request-dependent and subject to rights and location-sensitivity review**; it is not promised by this release. Provider retrieval can recover public portions, but cannot recover undocumented local decisions by inference. The inventory includes hashes of the dated preparation receipts, without publishing raw casework, personal storage paths or protected observation links.

The focal North/South comparison retains its **138/453-cell** cohort. The later **535-cell** South/Baja historical follow-up is a separate dataset and is not substituted. Likewise the earlier nonlinear development sensitivity used **136/453 cells** and is not a second reproduction of the primary comparison.

## CHELSA climate and environmental domains

### Contemporary climate

All three contemporary strands that use native climate draw on **CHELSA 2.1, 1981–2010 normals**, [dataset doi:10.16904/envidat.228](https://doi.org/10.16904/envidat.228), via the [official climatologies catalogue](https://www.chelsa-climate.org/datasets/chelsa_climatologies). The inventory lists 43 exact object identities: BIO01–BIO19 plus twelve monthly mean-temperature (`tas`) and twelve precipitation (`pr`) layers.

The saved object pattern is `https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/bioclim/bioNN/1981-2010/CHELSA_bioNN_1981-2010_V.2.1.tif`, with analogous `climatologies/tas/1981-2010/` and `climatologies/pr/1981-2010/` monthly paths. Each field has its recorded ETag and source byte size; monthly receipts also retain Last-Modified metadata. Recorded packing/conversion rules and inherited physical-array hashes identify the analysis inputs. A multipart object ETag is not asserted to be a whole-file MD5.

Saved raster-source evidence records **CC0-1.0**. This does not automatically assign identical rights to every supporting technical document. Regional source windows, retained byte ranges, aligned physical arrays, source grid and derived masks remain in the local archive. Matching a filename or BIO variable name alone is insufficient for exact replay: source object identity, unit conversion, grid alignment, period and missing-data support must match. The full-BIOCLIM reference contains thermal information absent from monthly mean tas/pr; it is not a nested timing-only ablation. Paleoclimate products are outside this contemporary source inventory.

### CEC and Natural Earth masks

The availability domains use **CEC (2021), Terrestrial Ecoregions Level III, Edition 2.0**, from the [official catalogue](https://www.cec.org/north-american-environmental-atlas/terrestrial-ecoregions-level-iii/). The two recorded official downloads are [source shapefile ZIP](https://www.cec.org/files/atlas_layers/1_terrestrial_ecosystems/1_06_3_terr_ecoregions_iii/terr_ecoregions_v2_level_iii_shapefile.zip) and [metadata ZIP](https://www.cec.org/files/atlas_layers/1_terrestrial_ecosystems/1_06_3_terr_ecoregions_iii/metadata.zip). The saved shapefile archive SHA-256 is `a99b4967e37b309651dfe117435573542d31576dafd779826f80874e015b3338`; bundled metadata/terms record **CC BY 4.0**.

The source inventory retains the selected object IDs as well as Level-III codes. Domain construction also requires the declared projected 25-km buffer, mainland clipping before and after buffering, pixel-center rasterization, climate validity mask and full native physical cell areas. The recorded CEC-to-WGS84 transformation was ballpark with unknown accuracy. Public source polygons are not measured dispersal boundaries or a substitute for the frozen masks.

The land support uses a **20 September 2026 hash-pinned Natural Earth `ne_10m_land.geojson` snapshot**, [recorded official repository URL](https://raw.githubusercontent.com/nvkelso/natural-earth-vector/master/geojson/ne_10m_land.geojson), SHA-256 `1ac90796408bc6ad6911d69448485d3c4dbf2190370080368a09976e1c9f7416`. Saved terms state **public domain**. The receipt does **not** establish a numbered land release or an immutable Git commit; the master URL is mutable. The retained original bytes and continental-component rule identify the actual support.

Figure S1's display-only administrative lines are a different Natural Earth acquisition: **Admin-0 and Admin-1 boundary lines, version 5.1.0**, with their exact URLs and hashes in the inventory. These are public-domain cartographic context, not model predictors. Exact S1 rebuilding still needs omitted native support and coordinate-bearing display inputs; the public aggregate replay checks the 3,914 displayed-cell count, not native geometry.

## Rose/USGS occurrence, water and survey inputs

### Released tables and public background

The source is **S35**, [doi:10.5066/P9O7MC7M](https://doi.org/10.5066/P9O7MC7M), [ScienceBase item `679916e7d34ea8c18376e9a5`](https://www.sciencebase.gov/catalog/item/679916e7d34ea8c18376e9a5). Saved metadata dates the data release to **25 February 2025**, with ScienceBase lastUpdated **26 February 2025, 20:27:19 UTC**. No numbered data-release edition is documented; the related code's version 1.0 is not the data version.

| File | Role and retained identity |
| --- | --- |
| `presence_only_records.csv` | 409 public occurrences, including dates and mapped-water distances. Primary water benchmark uses 324 records from 1995–2016; the 85 later records form a guarded secondary comparison. SHA-256 `7f5a7dd78d528e12d8ef121a8257392eb8ce83c2b7d78311f52a87e012c78306`. |
| `site_occupancy_data.csv` | 376 survey-cell histories and terrestrial/water covariates. SHA-256 `2923e2e610cf9c09d20190a9dbd6f8832facc7a17209f7ff5523fc278beeccaf`. |
| `so_detection_covariates.csv` | Annual trapping effort and precipitation, joined by explicit source IDs. SHA-256 `5ce7fe0c0b6f1864286f8819d9ed032ff672eb4d274bb6ec8066e5b37f3411c9`. |
| `spea_hammondii_southern_SDM.xml` | Original FGDC definitions and source projection. SHA-256 `e38d7b1cbe06b1c625336ced4c6f6427f8daa1c41d97845ab4fdd29e44e720de`. |
| `quadrature_points_1km.csv` | 28,931 public background points and supplied covariates/distances; benchmark retains 28,880 on established native land. SHA-256 `88655f1ccb2cefc474042b73f61a7b482956a5088cddeee650b5e3fe83064dd9`. |

Exact public ScienceBase file URLs, retrieval times and provider checksum comparisons are in the inventory. The three small data tables/XML were acquired in September; the public background grid was acquired on **2 October 2026**. Their original bytes remain local. Saved catalog evidence labels the release **public-domain-label-1.0**, with no access or use constraints; source-specific attribution and rights still apply. Raw occurrence points and survey histories are not bundled in the aggregate companion.

The survey model uses 2,133 surveyed cell-years and preserves 6,139 no-trapping years as missing, not nondetections. The released background `Id` and survey `cellID` have distinct documented roles. No explicit crosswalk or survey-cell coordinates were established. The survey comparison therefore contains no reconstructed coordinates or native BIOCLIM extraction.

### What the water variables do and do not provide

The tested `vp.pond.dist` and `ephem.dist` fields come directly from those released tables. [Rose et al. 2022, section 2.3](https://www.sciencedirect.com/science/article/pii/S2351989421004947), doi:10.1016/j.gecco.2021.e01944, documents separate distance calculations for presence-only and background points and covariate aggregation at survey arrays. This explicit methods description qualifies ambiguous trap-array wording retained in the original presence-only XML.

The article names regional vernal-pool inventories, NWI seasonal/temporary palustrine wetlands, NHD playas/intermittent ponds, and NHD ephemeral streams. **Individual layer versions and the complete geometry inventory were not acquired.** Reacquiring the released tables can recover the supplied distance values. Reconstructing the original geometry-derived distances is not established by this archive and would require separate source/version/access recovery. Inventory citation years are not pool observation dates; geometry rights remain source-specific and unresolved here. These distances are not hydroperiod, predator exclusion or exhaustive temporary-water availability.

### Official code as a methods reference

The inspected USGS software release is **S94**, [doi:10.5066/P9YYKW1H](https://doi.org/10.5066/P9YYKW1H), **version 1.0**, commit `d84c6fe0972a0d26a42e418b7396a47795e789cc`, code.usgs.gov project **11655**. The inventory gives exact commit-specific API URLs and SHA-256/Git-blob checks for `SDM_functions.R`, `prep_data.R` and `fit_models.R`. Saved license evidence is **CC0-1.0 / US public domain**.

These files informed the observation-process review; they were not executed as the manuscript model. The manuscript's latent-state diagnostic uses a separate implementation and is not a reproduction of the published integrated model. This code revision does not establish which version produced the 2022 paper.

## Omitted derived inputs and code

The [native code guide](native/README.md) identifies the original execution entry points. The [machine-readable code inventory](native/code_selection.json) lists every supplied implementation file and the omitted direct input/audit dependencies with saved identity hashes.

The public source routes do not by themselves regenerate the frozen scientific analysis. Exact native replay additionally needs accepted source-to-cell assignments and geographic rules; available and held masks; uncertainty-aware guards; sampled quadrature and inclusion probabilities; aligned physical arrays and cell areas; training transformations; fitted models and model-specific normalizers; and source-ID survey joins. Those derived assets are retained locally with original run manifests, but are not separate provider datasets with invented DOIs.

The public aggregate files and their replay are the intentionally portable reproducibility layer. Native access, if feasible, requires source-specific rights and sensitivity review. New downloads, new curation or changed software would require new receipts and validation rather than silently replacing the original inputs. The repository license applies only to the materials it actually licenses; it does not relicense third-party datasets or grant access to omitted private data.
