# Changelog

Versions are tagged `v<version>` (e.g. `v0.2.0`); the SKM web app pins a tag.

## 0.2.1 (2026-11-09)

Boolean models (TabularQual, BoolNet): corrected rules, and a new model fix.

- Boolean rules: several activations are bracketed before an inhibition, so the inhibition applies to all of them
  (it applied to the last activation only: 14 of the public model's 156 inhibited rules, e.g. EIN3).
- Boolean rules: a rule per target. A binding inhibition switches each partner off by the other partners, no
  longer by itself as well; a regulated process is no input of its own rule (Thermotolerance, Phenolics production
  and O2 deficiency could never switch on); a rule used by several reactions appears once.
- New model fix (`fix-form:activation-from-inactive`, "with model fixes" variant only): an active form that is only
  switched off (deactivated, degraded, bound) gets an activation from its inactive form when that is in the model.
  Without it, the Boolean models never switch it on again (15 activations in the public model, e.g. ARF, PP2C, JAM).
- Model fixes are deterministic: the location fixes are added in the same order every run.
- Tests for the Boolean rules of each reaction type.

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
