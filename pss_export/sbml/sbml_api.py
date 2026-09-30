from datetime import datetime
from html import escape

import libsbml
from libsbml import (SBMLDocument, writeSBMLToFile, writeSBMLToString,
                    LIBSBML_OPERATION_SUCCESS, OperationReturnValue_toString,
                    CVTerm, BIOLOGICAL_QUALIFIER, MODEL_QUALIFIER)

from ..entity_classes import IDTracker, Species, SpeciesType, SpeciesReference, Reaction
from ..annotations.annotation_manager import annotation_manager

SBML_LEVEL = 3
SBML_VERSION = 2
OUTSIDE_COMPARTMENTS = ['cytoplasm', 'extracellular']

#-------------------------------------
# libSBML helper stuff
#-------------------------------------

def check(value, message):
    """If 'value' is None, prints an error message constructed using
    'message' and then exits with status code 1.  If 'value' is an integer,
    it assumes it is a libSBML return status code.  If the code value is
    LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
    prints an error message constructed using 'message' along with text from
    libSBML explaining the meaning of the code, and exits with status code 1.
    """
    if value is None:
        raise Exception('sbml: LibSBML returned a null value trying to ' + message + '.')

    if type(value) is int:
        if value == LIBSBML_OPERATION_SUCCESS:
            return
        else:
            err_msg = 'sbml: Error encountered trying to ' + message + '.' \
            + 'LibSBML returned error code ' + str(value) + ': "' \
            + OperationReturnValue_toString(value).strip() + '"'
            raise Exception(err_msg)

    return

#-------------------------------------
# Annotations
#-------------------------------------

def to_cvterm(annotation):
    ''' The RDF CV term (qualifier and identifiers.org url) of an Annotation '''
    if annotation.namespace == "bqmodel":
        cv = CVTerm(MODEL_QUALIFIER)
        cv.setModelQualifierType(libsbml.ModelQualifierType_fromString(annotation.qualifier))
    else:
        cv = CVTerm(BIOLOGICAL_QUALIFIER)
        cv.setBiologicalQualifierType(libsbml.BiolQualifierType_fromString(annotation.qualifier))
    cv.addResource(annotation.url)
    return cv

#-------------------------------------
# SBML
#-------------------------------------

class SBML(SBMLDocument, IDTracker):

    def __init__(self, pss_adapter, kinetic_laws=False):
        '''
        Parameters
        ----------
        pss_adapter : PSSAdapter
            With the reactions collected: its node annotations and model
            metadata are used.
        kinetic_laws : bool
            Give each reaction an (empty) kinetic law with the SBO term of
            its rate law type.
        '''

        self.kinetic_laws = kinetic_laws
        self.pss_adapter = pss_adapter

        SBMLDocument.__init__(self, SBML_LEVEL, SBML_VERSION)
        IDTracker.__init__(self)

        # start the SBML model
        self.sbml_model = self.createModel()
        check(self.sbml_model, 'create model')
        check(self.sbml_model.setId(pss_adapter.model_id), 'set identifier on the Model object')
        self.add_model_metadata()

    def add_model_metadata(self):
        ''' Model name, notes (description, version, access, source, export
        date) and, if there are creators, the model history. '''

        adapter = self.pss_adapter
        model = self.sbml_model

        check(model.setName(adapter.model_name), 'set model name')
        check(model.setMetaId(f"metaid_{adapter.model_id}"), 'set model metaid')

        SBML.add_note(model, 'description', adapter.model_description)
        SBML.add_note(model, 'version', adapter.model_version)
        SBML.add_note(model, 'access', adapter.access)
        SBML.add_note(model, 'source', 'https://skm.nib.si')
        SBML.add_note(model, 'export date', adapter.export_datetime)

        if adapter.creators:
            history = libsbml.ModelHistory()
            for person in adapter.creators:
                creator = libsbml.ModelCreator()
                if person.family_name:
                    creator.setFamilyName(person.family_name)
                if person.given_name:
                    creator.setGivenName(person.given_name)
                if person.organization:
                    creator.setOrganization(person.organization)
                if person.email:
                    creator.setEmail(person.email)
                check(history.addCreator(creator), f'add model creator {person}')

            # libSBML needs a date in W3CDTF (seconds, time zone)
            export_date = datetime.fromisoformat(adapter.export_datetime or datetime.now().isoformat())
            date = libsbml.Date(export_date.astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")[:-2] + ":" + export_date.astimezone().strftime("%z")[-2:])
            check(history.setCreatedDate(date), 'set model created date')
            check(history.setModifiedDate(date), 'set model modified date')
            check(model.setModelHistory(history), 'set model history')

    def write(self, filename, replace_markup=True):
        if filename is None:
            return writeSBMLToString(self)
        return writeSBMLToFile(self, filename)

    def get_sbml_species_type(self, species_type):
        ''' Get a species type node in SBML model. Create if it does not exist.
        '''

        id_, status =  self.get_species_type_id(species_type)

        if status == 1:
            pass

        else:
            typ = self.sbml_model.createSpeciesType()
            check(typ, f'create species type {id_}\n')
            check(typ.setId(id_), f'set species type id {id_}\n')
            typ.setName(f'{species_type.name} {species_type.form}')
            # typ.setConstant(False)

            if species_type.sbo_term:
                typ.setSBOTerm(species_type.sbo_term)

            self.set_species_type_id(species_type, id_)
            species_type.set_id(id_)

        return id_

    def get_sbml_species(self, species):
        ''' Get a species node in SBML model. Create if it does not exist.

        Returns
        -------
            id: str
        '''

        species_id, status =  self.get_species_id(species) # look up the id in the IDTracker

        if status == 1:
            pass

        else:
            sp = self.sbml_model.createSpecies()
            sp.setId(species_id)
            sp.setMetaId(f"metaid_{species_id}")
            sp.setName(species.name)

            if (SBML_VERSION >= 2) and (SBML_LEVEL == 2):
            # TODO -- Error: Error: sbml: LibSBML returned a null value trying to create species type VPg_p.
                species_type = SpeciesType(name=species.name, form=species.form)
                specie_type_identifier = self.get_sbml_species_type(species_type)
                sp.setSpeciesType(specie_type_identifier)

            if species.sbo_term:
                sp.setSBOTerm(species.sbo_term)

            compartment_id = self.get_sbml_compartment(species.compartment)
            sp.setCompartment(compartment_id)

            sp.setHasOnlySubstanceUnits(False)
            sp.setBoundaryCondition(False)
            sp.setConstant(species.constant)

            self.add_species_annotations(sp, species)

            self.set_species_id(species, species_id)
            species.set_id(species_id)
            # print(f"SBML: species id: {species.name} --> {species_id}", species.compartment, species.form, species.sbo_term)

        return species_id

    def add_species_annotations(self, sp, species):
        ''' Database links (external links, functional cluster, Arabidopsis
        genes) as annotations, and description, form and additional
        information as notes. '''

        annotations, skipped_links = annotation_manager.process_node(
            self.pss_adapter.species_links(species.name))
        for skipped in skipped_links:
            print(f"SBML: warning, skipping or could not parse external link for species {species.name}: {skipped}")

        for annotation in annotations:
            check(sp.addCVTerm(to_cvterm(annotation)), f"add annotation {annotation.url} to species {species.name}")

        node_annotations = self.pss_adapter.node_annotations.get(species.name, {})
        SBML.add_note(sp, 'description', node_annotations.get("description"))
        SBML.add_note(sp, 'species form', species.form)
        SBML.add_note(sp, 'additional_information', node_annotations.get("additional_information"))

    def get_sbml_compartment(self, compartment):

        id_, status =  self.get_compartment_id(compartment)

        if status == 1:
            pass

        else:
            comp = self.sbml_model.createCompartment()
            comp.setId(id_)
            comp.setName(compartment)
            comp.setSize(1)
            comp.setConstant(True)

            if not (compartment in OUTSIDE_COMPARTMENTS):
                cyto_identifier = self.get_sbml_compartment('cytoplasm')
                m = comp.setOutside(cyto_identifier)
                check(m, "Set 'outside' of compartment")

            self.set_compartment_id(compartment, id_)

        # print(f"SBML: compartment id: {compartment} --> {id_}")
        return id_

    def create_reactant_reference(self, species, role, reaction):
        ''' Create SBML "SpeciesReference" in the model
        Returns the SBML species reference object '''

        id_ = self.get_sbml_species(species)
        species_reference = SpeciesReference(species, stoichiometry=1, role=role)
        reactant = reaction.createReactant()
        reactant.setSpecies(id_)
        reactant.setStoichiometry(species_reference.stoichiometry)
        reactant.setConstant(species.constant)
        if species_reference.sbo_term:
            reactant.setSBOTerm(species_reference.sbo_term)

        return reactant


    def create_product_reference(self, species, role, reaction):
        ''' Create SBML "SpeciesReference" in the model
        Returns the SBML species reference object '''

        id_ = self.get_sbml_species(species)
        species_reference = SpeciesReference(species, stoichiometry=1, role=role)
        product = reaction.createProduct()
        product.setSpecies(id_)
        product.setStoichiometry(species_reference.stoichiometry)
        product.setConstant(species.constant)
        if species_reference.sbo_term:
            product.setSBOTerm(species_reference.sbo_term)

        return product

    def create_modifier_reference(self, species, role, reaction):
        ''' Create SBML "ModifierReference" in the model
        Returns the SBML modifier reference object '''

        id_ = self.get_sbml_species(species)
        species_reference = SpeciesReference(species, stoichiometry=None, role=role)
        modifier = reaction.createModifier()
        modifier.setSpecies(id_)
        if species_reference.sbo_term:
            modifier.setSBOTerm(species_reference.sbo_term)

        return modifier

    def create_sbml_reaction(self, reaction):
        ''' Create SBML "Reaction" in the model
        Returns the reaction itself, not the identifier '''

        rxn = self.sbml_model.createReaction()
        rxn.setId(reaction.reaction_id)

        self.set_reaction_id(reaction, reaction.reaction_id)

        rxn.setMetaId(f"metaid_skm_{reaction.reaction_id}")

        if reaction.reaction_type_sbo:
            rxn.setSBOTerm(reaction.reaction_type_sbo)

        rxn.setReversible(False)
        rxn.setFast(False)

        if self.kinetic_laws:
            kinetic_law = rxn.createKineticLaw()
            check(kinetic_law, f'create kinetic law for reaction {reaction.reaction_id}\n')
            if reaction.kinetic_law_sbo:
                kinetic_law.setSBOTerm(reaction.kinetic_law_sbo)

        links = [f"skm:{reaction.reaction_id}"] + list(reaction.external_links or [])
        annotations, skipped_links = annotation_manager.process_node(links)
        for skipped in skipped_links:
            print(f"SBML: warning, skipping or could not parse external link for reaction {reaction.reaction_id}: {skipped}")
        for annotation in annotations:
            check(rxn.addCVTerm(to_cvterm(annotation)), f"add annotation {annotation.url} to reaction {reaction.reaction_id}")

        if reaction.evidence_sentence:
            SBML.add_note(rxn, 'evidence_sentence', reaction.evidence_sentence)

        if reaction.reaction_mechanism:
            SBML.add_note(rxn, 'mechanism', reaction.reaction_mechanism)

        if reaction.export_notes:
            SBML.add_note(rxn, 'export_notes', reaction.export_notes)

        return rxn

    @staticmethod
    def add_note(node, prefix, note):
        ''' Add a note to a node (reaction or species) '''

        if not note:
            return

        note = f"<body xmlns='http://www.w3.org/1999/xhtml'><p>{escape(prefix)}:{escape(str(note))}</p></body>"
        if node.isSetNotes():
            status = node.appendNotes(note)
        else:
            status = node.setNotes(note)

        return check(status, "add note to node")

    def add_reaction(self, reaction):

        # each reaction has ~4 nodes and ~3 arcs
        # substrate 1/+, product 1/+, process 1, modifier 1/+ nodes

        if reaction.reaction_type == 'unknown':
            # current_app.logger.info(f"{reaction_id}, undefined reaction type")
            print(f"SBML: {reaction.reaction_id}, unknown reaction type")

        # (1) create reaction object
        if reaction.reaction_id in self.reaction_ids:
            print(f"SBML: {reaction.reaction_id}, already exists")
            return -1

        rxn = self.create_sbml_reaction(reaction)

        '''
                      (modifier)
                          |
                          |
                          |
                          o
        (substrate)----[process]---->(product)
        '''

        # (2) substrate glyphs and arcs
        # (substrate)-[consumption]->(reaction)
        for species in reaction.substrates:
            self.create_reactant_reference(species, role=reaction.substrate_role, reaction=rxn)

        # (3) product glyphs and arcs
        # (reaction)-[production]->(product)
        for species in reaction.products:
            self.create_product_reference(species, role=reaction.product_role, reaction=rxn)

        # (4) modifier glyphs and arcs
        # (modifier)-[modifies]->(reaction)
        for species in reaction.modifiers:
            self.create_modifier_reference(species, role=reaction.modifier_role, reaction=rxn)
