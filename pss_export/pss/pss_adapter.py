'''
Exports ..... formats ..... of PSS model
'''
# library imports
from .config import pss_export_config, pss_schema_config
from .collectors import PSSCollector

# internal imports
# (ModelFixer is imported in model_fixes(), as it needs optional dependencies)
from ..graph_db import GraphDB, resolve_connection_settings
from ..entity_classes import Person

# # SBGN
# from .sbgn_api import SBGN

# SBML
from ..sbml import SBML

from ..boolean import TabularQual

# # projection for DiNAR
# from .pss_dinar_translation import pss_dinar_translation


from datetime import datetime

################################################################################
# Handles all exports, and has paths for endpoints to use
################################################################################
class PSSAdapter():
    '''Exports for PSS model.s
    Exports (will) include:
        - SBML
        - *SBGN
        - *DiNAR projection
        - *JSON (for API)
    '''

    def __init__(self,
                    neo4j_uri=None,
                    neo4j_user=None,
                    neo4j_password=None,
                    model_id=None,
                    model_name=None,
                    model_description=None,
                    creator=None):
        '''
        Constructor for PSSAdapter class.

        The database is only connected to while collecting data
        (collect_reactions), and the connection is closed afterwards.

        Parameters
        ----------
        neo4j_uri, neo4j_user, neo4j_password : str, optional
            Database connection settings. Any setting not given is read from
            the environment (MY_NEO4J_URI, MY_NEO4J_USER, MY_NEO4J_PASSWORD),
            or from a .env file in the current folder.
        model_id, model_name, model_description : str, optional
            Model metadata to use in exports.
        creator : list of str, optional
            Model creators, each in the format of:
            familyName | givenName | organization | email

        Returns
        -------
        None

        '''
        # raises ValueError now if any setting is missing
        self.connection_settings = resolve_connection_settings(
            uri=neo4j_uri, user=neo4j_user, pwd=neo4j_password)

        self.reactions = {}
        self.reaction_ids = []
        self.include_genes = None

        self.node_annotations = {}
        self.reaction_pathways = {}

        self.additional_reactions = []

        self.model_id = model_id or "pss_exported_model"
        self.model_name = model_name or "PSS Exported Model"
        self.model_description = model_description or "Model exported from the Plant Stress Signalling knowledge graph (PSS) available at https://skm.nib.si using the skm-pss-export package."

        if creator:
            self.creators = [Person(*creator.split("|")) for creator in creator] # expects format of: familyName | givenName | organization | email
        else:
            self.creators = []

        self.export_datetime = None

    def collect_reactions(self, **kwargs):
        ''' Collect reactions and annotations from the database.
        Optionally limit to specific pathways OR reactions.

        Has to be done on the fly, to include any updates made to the database.

        Parameters
        ----------
        **kwargs : dict
            Additional keyword arguments to pass to the PSSCollector.

        Returns
        -------
        None

        '''

        # reset data structures in case of re-collection
        self.reactions = {}
        self.reaction_ids = []
        self.node_annotations = {}
        self.additional_reactions = []

        print("Collecting reactions and annotations from the database...")

        # connection is only open while collecting
        with GraphDB(**self.connection_settings) as graph_db:

            collector = PSSCollector(graph_db, **kwargs)

            # collect reactions
            self.reactions = collector.collect_reactions()
            self.reaction_ids = list(self.reactions.keys())
            print(f"Collected {len(self.reaction_ids)} reactions.")

            # collect node annotations
            self.node_annotations = collector.collect_node_annotations()

            # collect reaction pathways (for SBGN)
            self.reaction_pathways = collector.collect_reaction_pathways(self.reaction_ids)

        self.export_datetime = datetime.now().isoformat()

        # needed for model fixes
        self.include_genes = collector.include_genes

    def model_fixes(self, interactive=False, apply_fixes=True):
        ''' Identify model fixes to the collected reactions.
            1) Fix node 'form' issues by changing input/outputs to active forms.
            2) Add transport reactions for species in multiple compartments.
        '''
        try:
            from ..model_fixes import ModelFixer
        except ImportError as e:
            raise ImportError(
                "Model fixing needs optional dependencies: "
                "pip install 'pss-export[model-fixing]'") from e

        ModelFixer(self, apply_fixes=apply_fixes, interactive=interactive).identify_model_fixes()

    def create_sbml(self,
                    access='public',
                    filename=None,
                    entities_table=None,
                    kinetic_laws=True):

        sbml = SBML(self, kinetic_laws=kinetic_laws)

        for reaction_id in self.reaction_ids:
            sbml.add_reaction(self.reactions[reaction_id])

        for reaction_id in self.additional_reactions:
            print(reaction_id)
            sbml.add_reaction(self.reactions[reaction_id])

        print("-" * 40)
        print("Number of species in SBML: ", len(sbml.species_ids))
        print("Number of species types in SBML: ", len(sbml.species_types_ids))
        print("Number of compartments in SBML: ", len(sbml.compartment_ids))
        print("Number of reactions in SBML: ", len(sbml.reaction_ids))
        print("-" * 40)

        if entities_table:
            sbml.write_entities_table(entities_table)

        return sbml.write(filename)

    def create_tabularqual(self, filename=None):
        '''  '''

        tabqual = TabularQual(self)

        for reaction_id in self.reaction_ids:
            tabqual.add_reaction(self.reactions[reaction_id])

        for reaction_id in self.additional_reactions:
            print(reaction_id)
            tabqual.add_reaction(self.reactions[reaction_id])

        tabqual.create_transitions()

        print("-" * 40)
        print("Number of species in TabularQual spreadsheet: ", len(tabqual.species_ids))
        print("Number of compartments in TabularQual spreadsheet: ", len(tabqual.compartment_ids))
        print("Number of transitions in TabularQual spreadsheet: ", len(tabqual.transitions))
        print("-" * 40)

        return tabqual.write(filename)
