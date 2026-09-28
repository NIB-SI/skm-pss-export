'''

'''

import re

from .pss.config import pss_export_config
from .pss.pss_reaction_definitions import reaction_classes, participant_roles, reaction_subtypes

#-------------------------------------
#  Helper classes (for nodes and reactions)
#-------------------------------------

class Reaction:
    """
    This is the "Reaction" for SBML, or the "Transition" for SBMLqual/TabularQual
    """
    def __init__(self, reaction_id, reaction_type, reaction_properties, include_conditions=False, include_genes=False, export_notes=None):
        self.id = reaction_id
        self.reaction_id = reaction_id
        self.reaction_type = reaction_type

        # string attributes
        self.reaction_mechanism = reaction_properties.get('reaction_mechanism', None)
        self.evidence_sentence = reaction_properties.get('evidence_sentence', None)
        self.reaction_effect = reaction_properties.get('reaction_effect', None)

        # list attributes
        self.external_links = reaction_properties.get('external_links', [])

        # export attribute
        if export_notes is None:
            export_notes = ""
        self.export_notes = export_notes

        self.substrates = []
        self.products = []
        self.modifiers = []

        # Figure out the roles of the reaction participants based on the reaction type
        # (the subtype and SBO terms depend on the participants, see properties below)
        participant_roles.assign_roles(self)

        # settings for preparing reactions
        self.include_conditions = include_conditions
        self.include_genes = include_genes

    def add_edges(self, edge_list):

        for path in edge_list:

            edge = path.relationships[0]
            edge_type = edge.type # edges only have 1 type

            if edge_type in ['SUBSTRATE', 'TRANSLOCATE_FROM']:

                if self.reaction_type in reaction_classes.TRANSCRIPTIONAL_TRANSLATIONAL and not self.include_genes:
                    continue
                    # TODO handle transcription and translation reactions differently

                # (1) source is SUBSTRATE
                key = 'source'
                name = edge.start_node['name']
                location = edge[f'{key}_location']
                form = edge[f'{key}_form']
                self.add_substrate(Species(name, form, location))

            elif edge_type in ['PRODUCT', 'TRANSLOCATE_TO']:

                # (2) target is PRODUCT
                key = 'target'
                name = edge.end_node['name']
                location = edge[f'{key}_location']
                form = edge[f'{key}_form']
                self.add_product(Species(name, form, location))

            elif edge_type in  ['INHIBITS',  'ACTIVATES']:

                # (1) source is MODIFIER
                key = 'source'
                name = edge.start_node['name']
                location = edge[f'{key}_location']
                form = edge[f'{key}_form']
                if form == "condition" and not self.include_conditions:
                    continue
                self.add_modifier(Species(name, form, location))

            else:
                # current_app.logger.info(f"{edge_type}, {reaction_id}")
                print(f"Reaction: {edge_type}, {self.reaction_id}")

    # The subtype and SBO terms describe the reaction as exported: they are
    # computed from the current participants (after export settings and model
    # fixes), as downstream tools generate rate laws from them.

    @property
    def reaction_subtype(self):
        """ Subtype based on type, mechanism, effect and current participants. """
        return reaction_subtypes.classify(self)

    @property
    def reaction_type_sbo(self):
        """ SBO term for the reaction (int), or None if not configured. """
        return self._subtype_SBO_term('reaction_type_SBO')

    @property
    def kinetic_law_sbo(self):
        """ SBO term for the kinetic law (int), or None if not configured. """
        return self._subtype_SBO_term('kinetic_law_SBO')

    def _subtype_SBO_term(self, key):
        terms = pss_export_config.reaction_subtype_to_SBO.get(self.reaction_subtype)
        if terms is None:
            return None
        return int(terms[key])

    def add_substrate(self, substrate):
        self.substrates.append(substrate)

    def add_product(self, product):
        self.products.append(product)

    def add_modifier(self, modifier):
        self.modifiers.append(modifier)

    def has_modifiers(self):
        return len(self.modifiers) > 0

    def has_substrates(self):
        return len(self.substrates) > 0

    def __repr__(self):
        return (f"Reaction(id={self.id}, "
                f"reaction_type={self.reaction_type}, "
                f"substrates={self.substrates}, "
                f"products={self.products}, "
                f"modifiers={self.modifiers})")

class  SpeciesType:
    def __init__(self, name, form=None):
        '''
        Parameters
        ----------
        name: str
            e.g. 'WRKY11'
        form: str
            e.g. 'protein', 'gene', 'complex'
        '''

        self.name = name
        self.form = form.lower()

        self.set_SBO_term()

    def set_id(self, id_):
        ''' Make a unique, short species type id
        '''
        self.id = id_
        return id_

    def set_SBO_term(self):
        """ Set the SBO term for the species type based on its form. """
        if self.form in pss_export_config.node_form_to_SBO:
            self.sbo_term = int(pss_export_config.node_form_to_SBO[self.form])
        else:
            self.sbo_term = None

    def __repr__(self):
        return f"SpeciesType(name={self.name}, form={self.form})"


class Species:
    def __init__(self, name, form, compartment):
        '''
        Parameters
        ----------
        name: str
            e.g. 'WRKY11'
        form: str
            e.g. 'protein', 'gene', 'complex'
        compartment: str
            e.g. 'cytoplasm', 'nucleus', 'extracellular'
        '''

        self.name = name
        self.form = form

        self.label = self.name.split("[")[0]

        if ((compartment is None) or (compartment == 'unknown')):
            # in case of cellular location, all nodes not assigned
            # are put within cytoplasm
            compartment = 'cytoplasm'
        compartment = compartment.replace("putative:", "")

        self.compartment = compartment

        if self.form == "gene":
            self.constant = True
        else:
            self.constant = False

        self.set_SBO_term()

    def set_id(self, id_):
        self.id = id_
        return id_

    def set_SBO_term(self):
        """ Set the SBO term for the species type based on its form. """
        if self.form in pss_export_config.node_form_to_SBO:
            self.sbo_term = int(pss_export_config.node_form_to_SBO[self.form])
        else:
            self.sbo_term = None

    def __repr__(self):
        return f"Species(id={self.id}, name={self.name}, form={self.form}, compartment={self.compartment})"

class SpeciesReference():
    def __init__(self, species, stoichiometry=1, role=None):
        '''
        Parameters
        ----------
        species: Species
            The species object.
        stoichiometry: int
            The stoichiometry of the species in the reaction. // not relevant for modifiers.
            Default is 1.
        '''
        self.species = species
        self.stoichiometry = stoichiometry
        self.role = role

        self.set_SBO_term()

    def set_SBO_term(self):
        """ Set the SBO term for the species role. """
        if self.role in pss_export_config.node_role_to_SBO:
            self.sbo_term = int(pss_export_config.node_role_to_SBO[self.role])
        else:
            self.sbo_term = None

    def __repr__(self):
        return f"SpeciesReference(name={self.name}, form={self.form}, compartment={self.compartment}, stoichiometry={self.stoichiometry}, role={self.role})"


#-------------------------------------
#  Naming/identifier functions
#-------------------------------------

class IDTracker:
    """ A class to track and create IDs for species, species_types, and reactions,
    to ensure unique IDs are generated for each.
    """

    def __init__(self, location=True, verbose=False):
        '''
        Parameters
        ----------
        location: bool
            Whether to include location in species IDs.
            Default is True.
        '''

        self.verbose = verbose
        self.location = location

        self.species_ids = {}
        self.reaction_ids = {}

        self.species_types_ids = {}
        self.compartment_ids = {}

        self.counters = {
            'species_type': 0,
            'species':0,
            'reaction': 0,
            'compartment': 0,
        }


    def write_entities_table(self, filename, delim="\t"):
        """ Writes a table of all entities with their IDs and attributes. """

        with open(filename, 'w') as f:
            f.write(f"id{delim}type{delim}name{delim}form{delim}compartment\n")

            for (name, form, compartment), id_ in self.species_ids.items():
                s = delim.join([
                    id_,
                    'species',
                    name,
                    form,
                    compartment,
                ])
                f.write(f"{s}\n")

            for (name, form), id_ in self.species_types_ids.items():
                s = delim.join([
                    id_,
                    'species_type',
                    name,
                    form,
                    '',
                    ''
                ])
                f.write(f"{s}\n")

            for compartment, id_ in self.compartment_ids.items():
                s = delim.join([
                    id_,
                    'compartment',
                    '',
                    '',
                    compartment
                ])
                f.write(f"{s}\n")

            for reaction_id, id_ in self.reaction_ids.items():
                s = delim.join([
                    id_,
                    'reaction',
                    reaction_id,
                    '',
                    ''
                ])
                f.write(f"{s}\n")

            return filename

    def get_species_id(self, species):
        '''
        Returns the ID of a species.
        If the species already has an ID in the tracker, it returns that ID and a status of 1.
        If the species does not have an ID, it creates a new ID,
        and returns the new ID with a status of 0.
        '''

        if not self.location:
            # ignore compartment for species ID by setting all of them to 'none'
            compartment = 'none'
        else:
            compartment = species.compartment

        if (id_ := self.species_ids.get((species.name, species.form, compartment))) is not None:
            return id_, 1

        id_ = f"{self.remove_nonalphanum(self.get_display_label(species.name))}"\
              f"_{pss_export_config.compartment_to_short[compartment]}"\
              f"_{pss_export_config.node_form_to_short[species.form]}"

        # make sure ID does not exist in idtracker.species_ids
        while id_ in self.species_ids.values():

            if self.verbose:
                # print all the species details to see why it possible duplicated
                print(f"Duplicate species ID found: {id_} for species {species}")

                # list the duplicated species
                for (name, form, compartment), existing_id in self.species_ids.items():
                    if existing_id == id_:
                        print("Existing species with same ID:")
                        print(f" - species.name: {name}")
                        print(f" - species.form: {form}")
                        print(f" - species.compartment: {compartment}")

            id_ += '_1'

        return id_, 0

    def get_reaction_id(self, reaction):
        return self.reaction_ids.get(reaction.id)

    def get_species_type_id(self, species_type):

        if (id_ := self.species_types_ids.get((species_type.name, species_type.form))) is not None:
            return id_, 1

        id_ = f"s_{IDTracker.remove_nonalphanum(IDTracker.get_display_label(self.name))}"\
              f"_{pss_export_config.node_form_to_short[self.form]}"

        # make sure ID does not exist in idtracker.species_types_ids
        while id_ in self.species_types_ids.values():

            if self.verbose:
                # print all the species_type details to see why it possible duplicated
                print(f"Duplicate species type ID found: {id_} for species type {species_type}")

                # list the duplicated species types
                for (name, form), existing_id in self.species_types_ids.items():
                    if existing_id == id_:
                        print("Existing species type with same ID:")
                        print(f" - species_type.name: {name}")
                        print(f" - species_type.form: {form}")

            id_ += '_1'

        return id_, 0

    def get_compartment_id(self, compartment):
        '''
        Returns the ID of a compartment.
        If the compartment already has an ID, it returns that ID and a status of 1.
        If the compartment does not have an ID, it creates a new ID,
        and returns the new ID with a status of 0.
        '''
        if (id_ := self.compartment_ids.get(compartment)) is not None:
            return id_, 1
        # create a new compartment ID
        return self.create_compartment_id(compartment), 0

    def set_species_id(self, species, id_):

        if not self.location:
            # ignore compartment for species ID by setting all of them to 'none'
            compartment = 'none'
        else:
            compartment = species.compartment

        self.species_ids[(species.name, species.form, compartment)] = id_
        self.counters['species'] += 1

    def set_reaction_id(self, reaction, id_):
        self.reaction_ids[reaction.id] = id_
        self.counters['reaction'] += 1

    def set_species_type_id(self, species_type, id_):
        self.species_types_ids[(species_type.name, species_type.form)] = id_
        self.counters['species_type'] += 1

    def set_compartment_id(self, compartment, id_):
        self.compartment_ids[compartment] = id_
        self.counters['compartment'] += 1

    @staticmethod
    def get_display_label(name):
        m = re.match(r"(.*)\[(.*)\]$", name)
        if m:
            new_name = m.groups()[0]
        else:
            new_name = name

        return new_name

    @staticmethod
    def remove_nonalphanum(name):
        ''' Remove non-alphanumeric characters from a string.
        Replace them with an underscore.
        This is used to create unique IDs.
        '''
        return re.sub('[^0-9a-zA-Z_]+', '', name)

    def create_compartment_id(self, compartment):
        ''' Make a unique, short compartment id
        '''

        id_ = f"c_{pss_export_config.compartment_to_short[compartment]}"

        # make sure ID does not exist in self.compartment_ids
        while id_ in self.compartment_ids.values():
            id_ += '_1'

        return id_


#---------------------
#  Other helper classes
#---------------------

class Person:

    def __init__(self, family_name=None, given_name=None, organization=None, email=None):

        self.family_name = family_name.strip() if family_name else None
        self.given_name = given_name.strip() if given_name else None
        self.organization = organization.strip() if organization else None
        self.email = email.strip() if email else None

    def __repr__(self):
        return f"Person(family_name={self.family_name}, given_name={self.given_name}, organization={self.organization}, email={self.email})"

