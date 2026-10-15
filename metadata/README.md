# Preserved scientific metadata

`mask_specification.json`, `method_specification.json`, `software_snapshot.json`, `references.bib`, and `references.json` are exact copies from the 4 October manuscript-preparation packet. Their original hashes are recorded in the root provenance file. Archive-relative `runs/`, `reports/`, and `src/spea/` identifiers inside these JSON records document original lineage; they do not point to downloadable assets in this repository.

The original method specification precedes the completed matched weighting and BIOCLIM synthesis. The [public Methods](../docs/methods.md#later-matched-representation-and-weighting-comparison) supplies that later design. The [working manuscript's references](../docs/manuscript.md#references) also include the two provider-definition web pages added with the 5 October BIOCLIM synthesis; they are not entries in the older 13-item BibTeX/JSON list.

Mask specifications describe domains, source polygon identifiers, grid geometry and a published-map-guided North/South assignment rule. They do not contain individual occurrence anchors or the accepted native-cell membership arrays. The North/South line is an operational geographic proxy, not a measured contact-zone boundary. Recorded dependency locks identify intended frozen software versions rather than a recovered installed environment for every historical fit.
