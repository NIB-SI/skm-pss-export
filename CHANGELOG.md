# Changelog

Versions are tagged `v<version>` (e.g. `v0.2.0`); the SKM web app pins a tag.

## 0.2.0 (2026-10-07)

All PSS exports of the SKM web app (skm.nib.si) are made by pss-export.

- New formats: the network exports (reaction graph, interaction network, gene network per species; extended SIF)
  (#14), the reaction graph as JSON, the data of the PSS Explorer (#20), FAIDARE (#18), BoolNet (#21), SBGN-ML and
  its Newt variant (#22). A format registry (`pss_export.formats`) with titles, descriptions and files.
- Models (SBML, TabularQual, BoolNet, SBGN) in two variants, without and with model fixes (#15); model fixes add
  Boolean rules for transport reactions (#12); `nodes_to_ignore` for the models only (#14, #17).
- Species filter for all exports (#14); gene network genes get their functional cluster's class (#19).
- SBML: species annotations, notes and model metadata (#13); annotations from identifiers.org CURIEs (#16); species
  ids start with `s_` (#11).
- Export correctness fixes; `PSSAdapter` owns the Neo4j connection (#10).

## 0.1.0

First version.
