# Changelog

Versions are tagged `v<version>` (e.g. `v0.2.0`); the SKM web app pins a tag.

## Unreleased

- SBML with `--include-genes` (`include_genes=True`): the gene of a transcription/translation is a modifier with
  the role template (SBO:0000645), after the regulators, no longer a reactant (it was consumed by its own
  transcription). Breaking for tools that read the gene as the reactant (e.g. rate laws with `S` = the gene).
  `Reaction.templates` holds the genes; the default export (without genes) is unchanged.
- Model fixes: genes are never transported between compartments (a modelling rule; a gene in two compartments is
  a data error).
- SBML: the PSS pathways as groups (SBML Level 3 Groups package, not required): a group per pathway with the
  species of its entities as members; the notes `pathway` and `all pathways` on each species. Optional:
  `create_sbml(pathway_groups=False)`, CLI `--no-pathway-groups` (the notes stay).
- SBML: compartment containment as `bqbiol:isPartOf` annotations (config `compartment_part_of`); the apoplast is
  part of extracellular. The Level 2 `outside` call, which Level 3 ignores, is removed.
- Several species: `collect_reactions(species=['ath', 'stu'], species_match='any'|'all')`, CLI `--species ath,stu
  --species-match any|all`: the reactions of any of the species (default; the single-species exports together) or of
  all of them (gene clusters with genes in each); the network node files then have a `<code>_homologues` column per
  species instead of `genes`. Not for the gene network. `PSSAdapter.species` is now a tuple of codes (`('ath',)`)
  or None, instead of a code (a code is still accepted wherever a species is read).

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
