'''
Exports ..... formats ..... of PSS model
'''
# library imports
from .config import pss_export_config, pss_schema_config
from .collectors import PSSCollector

# internal imports
# (ModelFixer is imported in model_fixes(), as it needs optional dependencies)
from ..entity_classes import Person

# # SBGN
# from .sbgn_api import SBGN

# SBML
from ..sbml import SBML

from ..boolean import TabluarQqual

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

    def __init__(self, graph_db,
                    model_id=None,
                    model_name=None,
                    model_description=None,
                    creator=None):
        '''
        Constructor for PSSAdapter class.

        Parameters
        ----------
        graph_db : GraphDB
            The graph database connection object.
        include_genes : bool, optional
            Whether to include gene information in the reactions. Default is False.
        nodes_to_ignore : list or str, optional
            Nodes to ignore during export. Can be a list of node names or a single string. Default is 'default', which ignores nodes defined in the configuration.
        implement_model_fixes : bool, optional
            Whether to apply model fixes after collecting data. Default is True.

        Returns
        -------
        None

        '''
        self.graph_db = graph_db

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

        collector = PSSCollector(self, **kwargs)

        # collect reactions
        self.reactions = collector.collect_reactions()
        self.reaction_ids = list(self.reactions.keys())
        print(f"Collected {len(self.reaction_ids)} reactions.")

        # collect node annotations
        self.node_annotations = collector.collect_node_annotations()

        # collect reaction pathways (for SBGN)
        self.reaction_pathways = collector.collect_reaction_pathways()

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

    def create_tabulrqual(self, filename=None):
        '''  '''

        tabqual = TabluarQqual(self)

        for reaction_id in self.reaction_ids:
            tabqual.add_reaction(self.reactions[reaction_id])

        for reaction_id in self.additional_reactions:
            print(reaction_id)
            tabqual.add_reaction(self.reactions[reaction_id])

        tabqual.create_transitions()

        print("-" * 40)
        print("Number of species in TabluarQqual spreadsheet: ", len(tabqual.species_ids))
        print("Number of compartments in TabluarQqual spreadsheet: ", len(tabqual.compartment_ids))
        print("Number of transitions in TabluarQqual spreadsheet: ", len(tabqual.transitions))
        print("-" * 40)

        return tabqual.write(filename)
