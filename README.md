# Spea climate and biogeography

**Recorded climatic differences and geographic prediction across Spea** is a working manuscript by Kevin M. Neal, with its evidence and analysis code. Public companion **v0.3.0, 29 October 2026**, brings the completed North/South, climate-representation and mapped-water comparisons into one paper. It supersedes the earlier seasonal manuscript; [release tags preserve that text](VERSION_HISTORY.md).

[Read the paper (PDF)](manuscript/Spea_Comparative_Manuscript.pdf) · [Full text](manuscript/spea_comparative_manuscript.md) · [Evidence index](manuscript/evidence_index.md)

## What the paper finds

Recorded climatic differences, extra environmental detail and habitat proxies do not necessarily improve geographic prediction. The paper examines that question through three complementary analyses:

- **North/South:** recorded climates differed, but adding separate climatic responses did not consistently improve prediction beyond the specified gradual geographic reference. Source and regularization sensitivity remain visible.
- **Climate representation:** monthly profiles improved on simple annual and variability controls in some units, especially intermontana. Full BIOCLIM had the higher primary mean in all five units under both training weights; some secondary scores reversed that comparison.
- **Mapped water:** distances improved the southern occurrence-density comparison but slightly worsened the separate survey-history comparison in every block. The designs differ in observations, controls and likelihoods, so their contrast does not isolate observation effort.

The results concern recording patterns and geographic prediction. They establish no new population boundary, physiological specialization or historical reconstruction. The five fitting units and reused development folds are not independent species-level replications. All primary folds, adverse results and consequential sensitivities remain available.

![Matched climate representation comparisons across the five Spea analysis units](manuscript/figures/figure_2_climate_representations.png)

*Monthly profiles compared with annual climate, variability summaries and full BIOCLIM under equal and soft spatial training weights. Points are retained folds and larger markers are equal-fold means; spread is descriptive. Full definitions and the other two comparisons are in the manuscript.*

## Run the reported arithmetic

After cloning, the integrated replay uses Python 3.13 and its standard library only:

```sh
git clone https://github.com/kevinmneal/spea-climate-research.git
cd spea-climate-research
python3 manuscript/replay/replay_aggregates.py --output replay_output/manuscript_receipt.json
```

It checks 3,664 aggregate arithmetic and source-identity comparisons across all three strands. The [reproduction guide](REPRODUCIBILITY.md) also provides figure/document rebuild commands and the preserved seasonal/environmental replays.

The [native implementation](native/README.md) exposes original fitters, transformations, fixed configurations and their project dependencies. It requires archived record preparation, grids, masks, parent runs and audit assets that are not all distributed here. The [source and access guide](SOURCES_AND_ACCESS.md) supplies precise provider/query/version identities and explains how public retrieval differs from exact archive restoration. In particular, the GBIF input has no collective download DOI, and some accepted source records and curation remain request-dependent.

## Navigate the companion

| Material | Contents |
| --- | --- |
| [Manuscript package](manuscript/README.md) | PDF, Markdown source, supporting information and all main figure inputs |
| [Reproducibility](REPRODUCIBILITY.md) | Runnable aggregate checks; figure and document rebuilding; native limits |
| [Native code](native/README.md) | Original fitting and transformation implementation, fixed configs and input inventory |
| [Sources and access](SOURCES_AND_ACCESS.md) | Source versions, queries, file identities, retrieval routes and omitted inputs |
| [Earlier aggregate modules](data/README.md) | Preserved seasonal sensitivities, decompositions and environmental-organization data |
| [Environmental organization](docs/environmental_organization.md) | Separate completed exploratory screen, including its negative result and sensitivity discrepancy |
| [Version history](VERSION_HISTORY.md) | Current manuscript and preserved predecessors |
| [Change log](CHANGELOG.md) | Changes by release |

The [broader research agenda](docs/research_agenda.md) is background to the completed paper. Historical, genetic and larval/weather investigations are not added to its main argument. Exact source records, survey-site crosswalks, downloaded papers and large native assets are outside this export. The manuscript's coverage map retains a 25-km display resolution; it is not an exact-locality release.

## Attribution and reuse

Maintained by **Kevin M. Neal**. Cite the release or commit used; [CITATION.cff](CITATION.cff) supplies repository citation metadata. This working manuscript has not been submitted to bioRxiv or presented as independently peer reviewed.

OpenAI Codex assisted with implementation, source triage, numerical checks, figures and drafting. Automated and AI-assisted review establish only the specific properties checked. They do not supply independent biological validation. Corrections and methodological discussion are welcome through repository issues; please keep sensitive locality information out of public discussions.

Original software is [MIT licensed](LICENSE). Original research prose, figures and aggregate tables are [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Third-party datasets and publications retain their own terms; see [NOTICE](NOTICE.md) and the source/access guide.
