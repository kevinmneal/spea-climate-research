# Methods and numerical definitions

Public snapshot v0.1.0. Definitions for the original comparison and subsequent sensitivities; the later matched comparison is specified below.

## Observation and climate support

The baseline observation is one accepted recorded native climate cell, rather than one source assertion, animal or pond. The dated preparation retained source identifiers, event dates, coordinate-support metadata and duplicate/event relationships. Its event-year scope was 1970–2026; this does not imply that every fitting unit was sampled throughout that interval. Coordinates and source assertions remain in the frozen local preparation. No new historical or Baja records enter these comparisons. Northern hammondii and the southern unit use the inherited 138/453-cell cohorts, not the later 535-cell southern preparation.

All original seasonal columns derive from the same 24 physical CHELSA 2.1 layers: January–December mean temperature in °C and precipitation totals in mm/month, representing 1981–2010 climatological normals. These are not occurrence-year weather, interannual moments or daily extremes. The provider's model and versioned dataset citations are Karger et al. (2017, 2021); see [references.bib](../metadata/references.bib). Aligned layers were not resampled. Physical arrays were float32; feature arithmetic used float64 before casting model columns to float32. Each feature row required finite values, nonnegative precipitation and positive annual precipitation; invalid values raised an error rather than silently dropping cells.

The global source grid is 20,880 × 43,200 at approximately 1/120 degree, EPSG:4326. The cached window is 4,801 × 4,561, with global row/column offsets 3,479/6,480. Exact affine transforms are in [mask_specification.json](../metadata/mask_specification.json). Global identity is `source_row*43200+source_col`; working flat identity is `local_row*4561+local_col`. The two identifiers are not interchangeable. Native cell areas are WGS84 physical areas indexed by row.

## Available geography and geographic assignment

The contemporary primary available domains use the frozen CEC (2021) North American Environmental Atlas Level III, Edition 2.0 source objects. The compact mask specification lists every selected CEC object ID and its Level III code; codes alone would not identify all source polygons. The unions were clipped to the independent Natural Earth mainland polygon containing the fixed continental anchor (−105°, 40°), buffered by 25 km in `+proj=aea +lat_1=22 +lat_2=48 +lat_0=35 +lon_0=-107 +datum=WGS84 +units=m +no_defs`, and clipped to the mainland again. Buffer geometry used 16 segments per quarter circle. The processing box [−140, 5, −75, 65] was checked not to truncate those geometries. Source islands were removed before buffering. The CEC source-to-WGS84 transformation is recorded as a ballpark transformation with unknown accuracy.

Raster membership used pixel centres (`all_touched=False`), followed by complete cached BIOCLIM availability. Included pixels contributed their full native physical area; partial-polygon areas were not substituted. Saved native masks are authoritative for replay. The 100-km variant belongs to the earlier background sensitivity and is not the primary mask in this seasonal comparison. These pragmatic regional choices were informed by documented geography and cohort coverage, rather than an occurrence hull, a fitted climate threshold or an optimized score. They are not measured accessible areas or LGM dispersal envelopes. These choices remain assumptions of the developmental analysis.

The North/South design inherits group-filtered domains and occurrence assignments from its `all_1km` working design. Its published-map-guided geographic proxy has projection EPSG:5070 and signed distance

```math
s(x,y)=\big((x,y)-(-2049132.5303109284,1544994.498919337)\big)
\cdot(-0.2614952118263048,0.9652047731916352).
```

Positive distance points north. The saved 10-km review half-band is an analyst review aid, not a measured contact-zone width or coordinate error. Individual reference anchors are not exported in this preparation. Replay uses the frozen species-cell assignments and masks; applying older footprint eligibility rules anew would create a different cohort. Neal et al. (2018) supplies the biological premise for treating the units separately, not a genetic identification of each occurrence.

## Spatial folds, support and integration

Intermontana, bombifrons and multiplicata retain five globally anchored lattice folds. Coordinates for their block calculations use `+proj=aea +lat_1=20 +lat_2=50 +lat_0=35 +lon_0=-105 +datum=WGS84 +units=m +no_defs`. With origin (0,0), the held fold is

```math
\left(\lfloor x/150000\rfloor+2\lfloor y/150000\rfloor\right)\bmod5.
```

Training and held geography retain 50-km clearance after source uncertainty and source-to-cell displacement. Unknown uncertainty uses the inherited 1-km training scenario; held cells belong to the common core eligible under the 5-km scenario. Known coordinate-support radii are not capped or set to zero. Background guards apply to whole native pixel footprints and account for cache boundaries; source-linked/event exclusions also remain in place. Native membership and record-support decisions are in the saved spatial arrays, not reconstructed from the displayed coarse map.

North/South retain three folds formed from whole occupied 50-km blocks ordered by northing, with cuts at the one-third/two-thirds quantiles of the occupied block-centre distribution within each group. Their inherited unknown-radius scenario is 1 km on both sides, and the guard is 30 km. These are different tasks from the five-fold lattice. The 50-km evaluation blocks also have a different role from the 150-km fold lattice. Following spatial-validation principles does not imply that blockCV generated these custom splits (Valavi et al. 2019).

Each inner split intersects outer training presences/domains with another saved fold's training membership. Its validation intersects outer training with that other fold's held membership; the inherited guards remain. Eligibility requires at least 20 training cells and three recorded validation blocks. Intermontana/bombifrons/multiplicata have four eligible inner splits per outer fold; North has [2,2,2], South [2,1,1]. Two southern candidate inner datasets with 17 training cells remain ineligible.

Training backgrounds are integration quadrature over available physical area, including recorded cells. Every training presence is forced into quadrature with inclusion probability one. Other eligible native cells are sampled without replacement, up to eight per 25-km stratum in the five-fold units and sixteen in North/South. If a stratum has $`N_h`$ eligible nonpresence cells and budget $`k`$, each nonpresence cell has $`\pi_h=\min(k,N_h)/N_h`$. Its quadrature weight is $`a_i/\pi_i`$, using physical native area $`a_i`$. Presence neighbourhoods are not removed from this integration. Source labels distinguish the recorded sample from the available-area sample; a background row is not a verified absence.

Exact saved quadrature rows, inclusion probabilities and weights take precedence over a new random draw. Genus outer draws used base 20260924 plus `100*original_species_index+fold`; original species order is hammondii/intermontana/bombifrons/multiplicata. Genus inner draws used that base plus `1000*index+10*outer+inner`. North/South outer quadrature is inherited from the earlier working comparison, whose quadrature base was 2026092749. Their inner quadrature base is 2026093079 plus `1000*NorthSouth_index+10*outer+inner`. A North/South model seed does not replace the inherited quadrature seed.

## Original 2/4/28-column representations

Let $`T_m`$ and $`P_m`$ be the twelve physical monthly values, $`\bar T=\sum T_m/12`$, $`A=\sum P_m`$, $`\bar P=A/12`$, $`t_m=T_m-\bar T`$, and $`q_m=P_m/A`$.

| Family | Exact columns | Dimension |
|---|---|---:|
| Annual | $`\bar T`$; $`\log(1+A)`$ | 2 |
| Variability | Annual columns; $`\max T_m-\min T_m`$; $`\sqrt{\sum(P_m-\bar P)^2/12}/\bar P`$ | 4 |
| Monthly | Variability columns; twelve $`t_m`$; twelve $`q_m`$, January–December | 28 |

The precipitation CV is a ratio with population denominator, not percentage BIO15. Monthly anomalies sum to zero and fractions sum to one; those redundant columns were retained. Calendar position and profile shape enter together. Latitude/longitude, environmental-cluster deduplication, PCA and UMAP are not ecological predictors or preprocessing steps here.

## Fitting, selection and effective damping

CatBoost fits used weighted binary Logloss/BCE. With $`N`$ training presences, positive rows each initially weighed one. Landscape quadrature weights were normalized to total $`N`$, giving total object mass $`2N`$. The same location can have a recorded-source row and an integration row with different source labels. Effective leaf damping was $`\lambda=rN`$; native `l2_leaf_reg` was $`\lambda`$ divided by mean object weight. An exactly profiled unpenalized intercept was not used.

| Candidate, in tie order | Trees | Depth | Learning rate | Damping coefficient $`r`$ |
|---|---:|---:|---:|---:|
| Baseline | 300 | 4 | 0.05 | 0.005 |
| Shallow | 150 | 2 | 0.03 | 0.03 |
| Moderate | 300 | 3 | 0.05 | 0.01 |
| Damped | 300 | 4 | 0.05 | 0.05 |

Quantization used 128 GreedyLogSum cutoffs from unweighted training object rows (presences plus quadrature), with one shared 28-column bank; smaller families received the identical cutoffs for common prefix columns. No held values supplied cutoffs or standardization. CPU fitting used two threads, Plain boosting, symmetric trees, no bootstrap, zero random strength, Newton leaf updates with ten estimation iterations and AnyImprovement backtracking, `boost_from_average=False`, `use_best_model=False`, `nan_mode=Min`, and `one_hot_max_size=255`. Predictors were numeric. Configured tree counts were completed without statistical early stopping. CatBoost's general software citation is Prokhorenkova et al. (2018); its categorical/ordered-boosting innovations are not an explanation established for these fits.

Selection maximized equal-inner-fold occupied-block log relative density independently within each unit/family/outer fold. Every candidate had to complete every eligible inner split, and at least two eligible splits were required. Ties within 1e−10 followed the candidate order above. With insufficient support, the declared baseline was retained. All eligible selections chose shallow; South outer folds 1/2 used baseline for all three families. Thus 126 fixed/selected family/fold slots correspond to 120 distinct outer models and six exact baseline aliases.

Original genus model seed base is 2026093051 with original hammondii/intermontana/bombifrons/multiplicata indices; North/South base is 2026093071 with North/South indices. Outer seed is `base+1000*index+fold`; inner seed is `base+1000*index+10*outer+inner`. Individual fitted models and their native input indices are retained in the working archive and are not distributed here. Saved floating-point parameters can differ slightly from intended decimal settings. The [software snapshot](../metadata/software_snapshot.json) records frozen lock versions; those locks are not a recovered per-fit installed-environment inventory.

## Score and diagnostic weights

For each held target, raw model margins $`\eta_j`$ were normalized over every complete saved native target cell with physical area $`a_j`$:

```math
\ell_i=\eta_i-\log\!\left(\frac{\sum_j a_j e^{\eta_j}}{\sum_j a_j}\right).
```

Zero is uniform recording density per physical area. The primary target selects a fold uniformly, an occupied 50-km block uniformly within that fold, then a recorded cell uniformly in that block. The secondary target selects a fold uniformly, then a recorded cell uniformly in that fold. Fold-specific block counts mean that blocks need not have equal weight across the whole study. Baseline training weights and evaluation weights have different roles. Each held target is normalized separately, so these scores assess allocation within held geography rather than total abundance between regions.

Source ROC AUC and average precision use block-weighted presences and area-weighted landscape, each with total mass one. AP therefore assumes an artificial 1:1 source mixture. Neither diagnostic, a sigmoid margin, nor the normalized score estimates occupancy. A held cell exceeding at least one training-background marginal feature range was flagged for support description; it was not excluded. The any-column rule gives larger families more chances for a flag, does not test joint support and does not explain a poor score by itself.

Within each reported run, inner choices were frozen before its outer evaluation. These same outer geographies had already contributed to wider project development and inspection. Their reuse prevents treating this manuscript as independent confirmation or treating fold signs as significance tests.

## Additional summaries and sampling sensitivities

The later fixed shallow comparison crossed 19 BIOCLIM variables, BIOCLIM plus fourteen summaries (33), compact variability plus summaries (18), monthly (28), and monthly plus summaries (42). BIOCLIM amount variables 12/13/14/16/17/18/19 used `log1p`; BIO15 was unchanged. A shared 61-column cutoff union (19 + 28 + 14) was constructed from the original training rows before thinning. These feature families omit the original annual-only and variability-only controls.

The fourteen added columns are below. Linear quantiles use NumPy's `method="linear"`; skewness is the population standardized third moment, zero for exactly constant input. Numerical rescaling by maximum absolute deviation does not change its mathematical definition. Set $`\theta_m=2\pi(m-1)/12`$, with January at angle zero.

| Added column | Definition / units |
|---|---|
| Temperature median offset | $`\operatorname{median}(T)-\bar T`$, °C |
| Temperature IQR | $`Q_{.75}(T)-Q_{.25}(T)`$, °C |
| Temperature skewness | Population third standardized moment, dimensionless |
| Precipitation median ratio | $`\operatorname{median}(P)/\bar P`$ |
| Precipitation IQR ratio | $`(Q_{.75}(P)-Q_{.25}(P))/\bar P`$ |
| Precipitation skewness | Population third standardized moment, dimensionless |
| Precipitation entropy | $`-\sum q_m\log q_m/\log12`$, with $`0\log0=0`$ |
| Wettest three-month share | Maximum consecutive circular sum of three $`q_m`$ |
| Driest three-month share | Minimum consecutive circular sum of three $`q_m`$ |
| Precipitation sine harmonic | $`\sum q_m\sin\theta_m`$ |
| Precipitation cosine harmonic | $`\sum q_m\cos\theta_m`$ |
| Temperature sine harmonic | $`\sum t_m\sin\theta_m/6`$, °C |
| Temperature cosine harmonic | $`\sum t_m\cos\theta_m/6`$, °C |
| Precipitation–temperature offset | $`\sum q_m t_m`$, °C |

All are transformations of the same monthly normals. The five-family equal-cell/5-km/10-km comparison produced 315 fixed shallow fits. Its model seed was 2026100314 plus `100*species_index+fold`, with order intermontana/bombifrons/multiplicata/North/South. Training-only hard thinning used WGS84 geodesic cell-centre distance. Conflicts had strict distance less than the stated radius; equality was allowed. Each attempt repeatedly removed a point of maximum current neighbour degree, with seeded random ties. Ten attempts per radius used `2026100320+1000*species_index+20*fold+attempt`; the first maximum-cardinality set was frozen before fitting. This heuristic does not guarantee a globally largest set. Slots still required 20 training cells and three training blocks. The spThin manual supplies algorithm context; the R package was not executed. The later 2-km follow-up fitted only BIOCLIM/monthly, adding 42 fits.

Hard thinning retained the unchanged quadrature, including discarded recorded locations, and the original held targets and guards. Effective damping scaled with retained training-cell count. No size-matched random-deletion control separates loss of sample size from change in geographical distribution. No radius was adopted because of its outer score.

The separate 42-fit soft-weighting comparison retained every original presence. With WGS84 distance $`d_{ij}`$, raw weights were $`u_i=1/(1+\sum_{j\ne i}\max(1-d_{ij}/10\text{ km},0))`$, then $`w_i=Nu_i/\sum u_i`$. Original positive mass, background mass, damping, cutoffs and targets remained matched. The inspected enmSdmX `weightByDist` formula was translated to Python; the R package was not run. Kish effective sample size describes weight concentration, not independent pond count. Neither taper nor hard thinning recovers unknown sampling effort. These initial interventions did not include matched annual/variability families. The later matched soft-weight comparison below supplies those controls for one fixed procedure; hard thinning remains outside it.

## Later matched representation and weighting comparison

The later comparison fixes the shallow learner (150 trees, depth 2, learning rate 0.03, leaf damping 0.03 times presence mass) in every eligible outer fold. It crosses annual (2 columns), variability (4), monthly (28) and full BIOCLIM (19) with equal-cell or soft 10-km training weights. This yields 168 family/arm/fold slots across the same 21 unit/fold tasks. Annual/variability fits and exact monthly reuses comprise 126 slots; the subsequent synthesis adds 42 existing BIOCLIM slots. These are reused developmental geographies, not new test sets.

All original presence cells, quadrature and held targets remain fixed. Shared columns retain common training-only quantization cutoffs; source alignment was checked in the original synthesis. The model seed is 2026100314 + 100 × species index + fold, in intermontana, bombifrons, multiplicata, North, South order. Every arm has the same fixed learner settings, and training weights retain total presence mass. This differs from the original selected/fallback procedure and its 300-tree depth-four fixed baseline.

The 19 BIOCLIM fields are provider rasters on the same CHELSA 2.1 1981–2010 grid. Variables 12/13/14/16/17/18/19 are transformed with log1p. BIOCLIM includes diurnal and minimum/maximum thermal descriptors unavailable in monthly mean-temperature inputs. Its comparison is non-nested and cannot isolate calendar timing, compression, or a causal role for extremes.

Each contrast subtracts metrics within the same unit, fold and weight arm, then averages folds equally. Minimum, maximum and positive-fold counts describe that finite set; no significance test or independent replication is implied. [All values and their definitions](../data/README.md) accompany the [aggregate replay](../REPRODUCIBILITY.md).

## Replay boundary

The portable replay reconstructs aggregate table arithmetic from published fold scores and checks the saved summaries. It does not rebuild source records, native physical-area normalizers, predictions, or fitted models. The [reproducibility statement](../REPRODUCIBILITY.md) separates these levels. The intermontana geographic and historical work is a separate exploratory module.
