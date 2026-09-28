'''
Use https://github.com/sys-bio/TabularQual to export table

'''

from collections import defaultdict
from importlib import resources
from typing import Dict, Any, Tuple, List, Set, Protocol, Optional

from tabularqual.types import QualModel, ModelInfo
from tabularqual.types import Species as TabularQualSpecies

from tabularqual.types import Transition as TabularQualTransition
from tabularqual.types import InteractionEvidence as TabularQualInteractionEvidence
from tabularqual.types import Person as TabularQualPerson

from tabularqual.spreadsheet_writer import write_spreadsheet

from ..entity_classes import IDTracker, Species, SpeciesType, SpeciesReference, Reaction
from .boolean import reaction_rule_constructor, rule_composer

from ..annotations.annotation_manager import annotation_manager

class TabularQualAnnotationStrategy:
    """Formats node annotations specifically for TabularQual layouts independently."""
    def format_node(self, valid_records: List[Dict[str, Any]]) -> List[Tuple[str, str]]:
        tabular_entries = []
        for r in valid_records:
            raw_qualifier = r.get("qualifier", "bqbiol:isVersionOf")
            clean_qualifier = raw_qualifier.split(":", 1)[-1] if ":" in raw_qualifier else raw_qualifier
            db_identifier = f"{r['canonical_prefix']}:{r['local_id']}"
            tabular_entries.append((clean_qualifier, db_identifier))
        return tabular_entries

annotation_manager.register_export_strategy("tabularqual", TabularQualAnnotationStrategy())

#-------------------------------------
# TabularQual
#-------------------------------------

class TabularQual(IDTracker):

    def __init__(self, pss_adapter):
        '''


        '''

        IDTracker.__init__(self)

        self.pss_adapter = pss_adapter

        self.rules = defaultdict(lambda: {"activation":[], "inhibition":[]})
        self.rules_rx = defaultdict(list)

        self.species_dict = {}
        self.transitions = []
        self.interactions = []

    def to_QualModel(self):
        return QualModel(
            model=self.create_model_info(),
            species=self.species_dict,
            transitions=self.transitions,
            interactions=self.interactions
        )

    def write(self, filename):
        ''' Write TabularQual spreadsheet to file '''

        qualmodel = self.to_QualModel()

        # use the template in "resources" folder
        # resources/tabular_qual_pss_template.xlsx

        template_path = None #resources.files('pss_export.resources') / 'tabular_qual_pss_template.xlsx'

        # print(template)
        # print(template.exists())

        write_spreadsheet(qualmodel, filename, template_path=template_path)

    def create_model_info(self):

        """Prepare model-level information"""

        model_id = self.pss_adapter.model_id
        name = self.pss_adapter.model_name

        notes = [self.pss_adapter.model_description]
        versions = ["1.0.0"]

        source_urls = ["https://skm.nib.si"]
        described_by = [] # "Publication"
        derived_from = [] # "Origin_publication"
        biological_processes = []
        taxons = []
        created_iso = self.pss_adapter.export_datetime
        modified_iso = None
        creators = [TabularQualPerson(
            family_name=creator.family_name,
            given_name=creator.given_name,
            organization=creator.organization,
            email=creator.email
        ) for creator in self.pss_adapter.creators]
        contributors = []

        return ModelInfo(
            model_id=model_id,
            name=name,
            source_urls=source_urls,
            described_by=described_by,
            derived_from=derived_from,
            biological_processes=biological_processes,
            taxons=taxons,
            created_iso=created_iso,
            modified_iso=modified_iso,
            creators=creators,
            contributors=contributors,
            versions=versions,
            notes=notes
        )


    def get_tabularqual_species(self, species):
        """Prepare species-level information"""

        species_id, status =  self.get_species_id(species) # look up the id in the IDTracker

        if status == 1:
            species.set_id(species_id)
            return species_id

        # annotations... use fc id and "external_links" (list of <db.:<id>), parse to (qulaifier, db:id)

        # Collect raw links array from the adapter
        links = self.pss_adapter.node_annotations.get(species.name, {}).get("external_links", [])

        if links is None:
            links = []

        # If a functional cluster ID exists, inject it as a standard <db>:<id> reference
        functional_cluster_id = self.pss_adapter.node_annotations.get(species.name, {}).get("functional_cluster_id", None)
        if functional_cluster_id:
            links.append(f"skm:{functional_cluster_id}")

        # add tair from "ath_homologues" if exists, as link to "tair"
        ath_homologues = self.pss_adapter.node_annotations.get(species.name, {}).get("ath_homologues", [])
        if ath_homologues:
            for ath_homologue in ath_homologues:
                links.append(f"tair:{ath_homologue}")

        # Process the entire reference array using the pre-warmed TabularQual strategy
        # Unrecognized or malformed links will naturally fall into the 'skipped_links' array
        annotations, skipped_links = annotation_manager.process_node("tabularqual", links)

        # Optional: Print tracking alerts for your skipped items
        for skipped in skipped_links:
            print(f"TabularQual: warning, skipping or could not parse external link for species {species.name}: {skipped}")


        notes = []

        # "description" as note 1
        description = self.pss_adapter.node_annotations.get(species.name, {}).get("description", None)
        if description:
            notes.append((f"description:{description}"))

        # "form as note 2
        form = species.form
        if form:
            notes.append((f"species form:{form}"))

        # "additional_information" as note 3
        additional_information = self.pss_adapter.node_annotations.get(species.name, {}).get("additional_information", None)
        if additional_information:
            notes.append((f"additional_information:{additional_information}"))

        tabqual_species = TabularQualSpecies(
            species_id=species_id,
            name=species.label,
            compartment=species.compartment,
            constant=species.constant,
            initial_level=None,
            max_level=None,
            annotations=annotations,
            notes=notes,
        )

        self.set_species_id(species, species_id)
        species.set_id(species_id)

        self.species_dict[species_id] = tabqual_species

        # print(f"TabularQual: species id: {species.name} --> {species_id}", species.compartment, species.form, species.sbo_term)


        return species_id


    def add_reaction(self, reaction):

        # each reaction has ~4 nodes and ~3 arcs
        # substrate 1/+, product 1/+, process 1, modifier 1/+ nodes

        if reaction.reaction_type == 'unknown':
            # current_app.logger.info(f"{reaction_id}, undefined reaction type")
            print(f"TabularQual: {reaction.reaction_id}, unknown reaction type")

        # (1) create reaction object
        if reaction.reaction_id in self.reaction_ids:
            print(f"TabularQual: {reaction.reaction_id}, already exists")
            return -1

        # (2) substrate glyphs and arcs
        # (substrate)-[consumption]->(reaction)
        for species in reaction.substrates:
            # print("TabularQual: substrate species:", species.name)
            self.get_tabularqual_species(species)

        # (3) product glyphs and arcs
        # (reaction)-[production]->(product)
        for species in reaction.products:
            # print("TabularQual: product species:", species.name)
            self.get_tabularqual_species(species)

        # (4) modifier glyphs and arcs
        # (modifier)-[modifies]->(reaction)
        for species in reaction.modifiers:
            # print("TabularQual: modifier species:", species.name)
            self.get_tabularqual_species(species)

        rule_constructor = reaction_rule_constructor(reaction)
        if rule_constructor is None:
            print(f"TabularQual: {reaction.reaction_id}, could not construct reaction rule")
            return -1

        targets, reaction_rule = rule_constructor(reaction)

        # print("TabularQual: reaction rule:", reaction.reaction_id, "targets:", targets, "rule:", reaction_rule)

        if reaction_rule is None:
            print(f"TabularQual: {reaction.reaction_id}, no reaction rule generated")
            return

        for target in targets:
            self.rules[target][reaction.reaction_effect].append(reaction_rule)
            self.rules_rx[target].append(reaction.reaction_id)

    def create_transitions(self):
        ''' Create TabularQual Transitions from the collected rules '''

        for species_id, rule_dict in self.rules.items():

            # print("TabularQual: creating transition for species:", species_id)

            activation_rules = rule_dict["activation"]
            inhibition_rules = rule_dict["inhibition"]

            update_function = rule_composer(species_id, activation_rules, inhibition_rules)

            transition_id = f"tr_{species_id}"

            # reactions in this transition
            reactions = [self.pss_adapter.reactions[reaction_id] for reaction_id in self.rules_rx[species_id]]

            # annotations... use "external_links" (list of <db.:<id>), parce to (db, id)
            links = []
            for reaction in reactions:
                links.extend(reaction.external_links or [])
                # add reaction_id as link to skm
                links.append(f"skm:{reaction.reaction_id}")
            annotations, skipped_links = annotation_manager.process_node("tabularqual", links)

            # Optional: Print tracking alerts for your skipped items
            for skipped in skipped_links:
                print(f"TabularQual: warning, skipping or could not parse external link for species {species_id}: {skipped}")

            # notes -- generate a note based on how the transition was produced (e.g. if multiple reactions combined...)
            if len(reactions) > 1:
                notes = [f"Transition rule composed from {len(reactions)} reactions: {', '.join(self.rules_rx[species_id])}"]
            else:
                notes = [f"Transition rule generated from reaction: {self.rules_rx[species_id][0]}"]

            transition = TabularQualTransition(
                transition_id=transition_id,
                target=species_id,
                name=None,
                level=1,
                rule=update_function,
                annotations=annotations,
                notes=notes
            )

            self.transitions.append(transition)

