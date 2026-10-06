'''
Exports ..... formats ..... of PSS model
'''
# library imports
import copy

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

from ..boolean import TabularQual, create_boolnet

from .. import networks
from .. import faidare

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
                    model_version=None,
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
        model_version : str, optional
            Version of the exported model, e.g. the PSS version.
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

        self.nodes = {}
        self.reaction_pathways = {}

        self.additional_reactions = []
        self.model_fixes_applied = None     # number of model fixes applied (None: none run)
        self._with_model_fixes = None

        self.model_id = model_id or "pss_exported_model"
        self.model_name = model_name or "PSS Exported Model"
        self.model_description = model_description or "Model exported from the Plant Stress Signalling knowledge graph (PSS) available at https://skm.nib.si using the skm-pss-export package."
        self.model_version = model_version

        if creator:
            self.creators = [Person(*creator.split("|")) for creator in creator] # expects format of: familyName | givenName | organization | email
        else:
            self.creators = []

        self.export_datetime = None
        self.access = None
        self.species = None

    def collect_reactions(self, species='ath', **kwargs):
        ''' Collect reactions and annotations from the database.
        Optionally limit to specific pathways OR reactions, and to a species.

        Has to be done on the fly, to include any updates made to the database.

        Parameters
        ----------
        species : str or None
            Limit to the reactions whose functional clusters all have genes in this species (a code
            of pss_schema_config.species, default 'ath'; reactions without functional clusters are
            kept). The gene annotations from the homologue lists (tair: for ath) and the gene network
            are then for this species. None: no species filter (no gene network).
        **kwargs : dict
            Additional keyword arguments to pass to the PSSCollector.

        Returns
        -------
        None

        '''

        # reset data structures in case of re-collection
        self.model_fixes_applied = None
        self._with_model_fixes = None
        self.reactions = {}
        self.reaction_ids = []
        self.nodes = {}
        self.additional_reactions = []

        print("Collecting reactions and annotations from the database...")

        # connection is only open while collecting
        with GraphDB(**self.connection_settings) as graph_db:

            collector = PSSCollector(graph_db, species=species, **kwargs)
            self.species = collector.species

            # collect reactions
            self.reactions = collector.collect_reactions()
            # sorted, so that the exports are in the same order every time
            self.reaction_ids = sorted(self.reactions)
            print(f"Collected {len(self.reaction_ids)} reactions.")

            # collect the entities and their annotations
            self.nodes = collector.collect_nodes()

            # collect reaction pathways (for SBGN)
            self.reaction_pathways = collector.collect_reaction_pathways(self.reaction_ids)

        self.export_datetime = datetime.now().isoformat()
        self.access = collector.access
        self.nodes_to_ignore = collector.nodes_to_ignore     # models only

        # needed for model fixes
        self.include_genes = collector.include_genes

    @property
    def model_reaction_ids(self):
        ''' The collected reactions the models (SBML, TabularQual) include: not those that the nodes to
        ignore (pss_export_config.yaml nodes_to_ignore) leave without substrates or products. The networks
        include all collected reactions. '''
        return [r for r in self.reaction_ids if self.reactions[r].in_model]

    @property
    def species_description(self):
        ''' e.g. "Solanum tuberosum (stu)", or "all species" (no species filter) '''
        if self.species is None:
            return "all species"
        return f"{pss_schema_config.species[self.species]['name']} ({self.species})"

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

        fixes = ModelFixer(self, apply_fixes=apply_fixes, interactive=interactive).identify_model_fixes()
        if apply_fixes:
            self.model_fixes_applied = (self.model_fixes_applied or 0) + fixes
        return fixes

    def with_model_fixes(self):
        ''' A copy of the collected data with the model fixes applied (pss_export.model_fixes: translation
        products in the active form, transport reactions between compartments), for the "with model fixes"
        variant of the models; this adapter stays as collected (exactly PSS). Made once per collection. '''
        if self._with_model_fixes is None:
            fixed = copy.deepcopy(self)
            fixed.model_fixes(apply_fixes=True, interactive=False)
            self._with_model_fixes = fixed
        return self._with_model_fixes

    def create_sbml(self,
                    filename=None,
                    entities_table=None,
                    kinetic_laws=False):

        sbml = SBML(self, kinetic_laws=kinetic_laws)

        for reaction_id in self.model_reaction_ids:
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

    def boolean_model(self):
        ''' The Boolean model of the collected reactions (the models' reactions and any added by the model
        fixes): a TabularQual with its species and transitions (rules), for TabularQual and BoolNet '''
        tabqual = TabularQual(self)
        for reaction_id in self.model_reaction_ids + self.additional_reactions:
            tabqual.add_reaction(self.reactions[reaction_id])
        tabqual.create_transitions()
        return tabqual

    def create_boolnet(self, filename, nodes_file=None):
        ''' The Boolean model in the BoolNet format (the TabularQual model's rules; inputs keep their value),
        and optionally a node file. Returns the number of rules. '''
        return create_boolnet(self.boolean_model(), filename, nodes_file)

    def create_tabularqual(self, filename=None):
        ''' The Boolean model as a TabularQual spreadsheet '''

        tabqual = self.boolean_model()

        print("-" * 40)
        print("Number of species in TabularQual spreadsheet: ", len(tabqual.species_ids))
        print("Number of compartments in TabularQual spreadsheet: ", len(tabqual.compartment_ids))
        print("Number of transitions in TabularQual spreadsheet: ", len(tabqual.transitions))
        print("-" * 40)

        return tabqual.write(filename)

    def create_reaction_graph(self, edges_file=None, nodes_file=None):
        ''' The reaction graph (extended SIF): entities and reactions, one edge per participant.
        Returns the number of edges. '''
        return networks.create_reaction_graph(self, edges_file, nodes_file)

    def create_interaction_network(self, edges_file=None, nodes_file=None):
        ''' The interaction network (extended SIF): entity -> entity influences through the
        reactions (rules: interaction_rules in pss_export_config.yaml). Returns the number of edges. '''
        return networks.create_interaction_network(self, edges_file, nodes_file)

    def create_gene_network(self, edges_file=None, nodes_file=None):
        ''' The gene network (extended SIF, the interaction network's columns): the interaction
        network with functional clusters expanded into their genes in the species collected for
        (collect_reactions(species=...); not without a species). Returns the number of edges. '''
        return networks.create_gene_network(self, edges_file, nodes_file)

    def create_faidare(self, filename):
        ''' The FAIDARE data discovery file (JSON): an entry per gene and species of the gene clusters in the
        collected reactions; needs a collection without a species filter (collect_reactions(species=None)).
        Returns the number of entries. '''
        return faidare.create_faidare(self, filename)

    def export(self, format_key, model_fixes=False, **arguments):
        ''' Write a format of the registry (pss_export.formats) by its key, e.g.
        export("gene-network", edges_file="e.tsv", nodes_file="n.tsv"). The arguments are the
        format's file arguments. model_fixes: the "with model fixes" variant (the formats with
        ExportFormat.model_fixes), from with_model_fixes(); both variants come from one collection. '''
        from ..formats import get_format

        fmt = get_format(format_key)
        allowed = {f.argument for f in fmt.files}
        unknown = set(arguments) - allowed
        if unknown:
            raise ValueError(f"'{format_key}' takes {', '.join(sorted(allowed))}, not {', '.join(sorted(unknown))}")
        if self.access is not None and self.access not in fmt.access:
            raise ValueError(f"'{format_key}' is not available for access '{self.access}'")
        if model_fixes and not fmt.model_fixes:
            raise ValueError(f"'{format_key}' has no variant with model fixes")
        adapter = self.with_model_fixes() if model_fixes else self
        return getattr(adapter, fmt.method)(**arguments)
