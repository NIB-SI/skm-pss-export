'''
Network exports: the reaction graph, the interaction network and the gene network.

All three are tab-separated edge and node files, the edges in extended SIF (as Pathway
Commons): a header, and the first three columns are source, interaction, target.

- Reaction graph: bipartite, entities and reactions, one edge per participant (lossless).
- Interaction network: entity -> entity influences through reactions, by the rules in
  pss_export_config.yaml (interaction_rules).
- Gene network: the interaction network in one species, functional clusters expanded into
  their genes. The same edge columns as the interaction network.
'''

from itertools import product as pairs

from ..entity_classes import Node
from ..pss.config import pss_export_config, pss_schema_config
from ..pss.pss_reaction_definitions import edge_types, participant_roles
from ..utils import sbo, write_tsv

# the rule groups (interaction_rules) that are modifiers: their influence on their own entity is kept
MODIFIER_GROUPS = ['modifiers', 'activating_modifiers', 'inhibiting_modifiers']

# rank as in CKN (0: best supported, to 4): PSS reactions are curated from the literature
PSS_RANK = 0

INTERACTION_EDGE_COLUMNS = [
    'source', 'interaction', 'target', 'directed', 'rank',
    'source_entity', 'target_entity', 'source_type', 'target_type',
    'source_role', 'target_role',
    'source_location', 'target_location', 'source_location_putative', 'target_location_putative',
    'source_form', 'target_form',
    'reaction_id', 'reaction_type', 'reaction_effect', 'reaction_mechanism',
    'evidence_sentence', 'experimental_techniques',
    'reaction_sbo', 'reaction_mechanism_sbo', 'influence_sbo',
]

ENTITY_NODE_COLUMNS = ['id', 'node_type', 'display_label', 'short_name', 'synonyms', 'description', 'additional_information',
                       'pathway', 'all_pathways', 'mapman', 'functional_cluster_id', 'components', 'component_cluster_ids',
                       'external_links']

GENE_NODE_COLUMNS = ['id', 'node_type', 'species', 'display_label', 'short_name', 'synonyms', 'description',
                     'additional_information', 'pathway', 'all_pathways', 'mapman', 'functional_cluster_id',
                     'components', 'component_cluster_ids', 'external_links']

REACTION_GRAPH_EDGE_COLUMNS = ['source', 'role', 'target', 'directed', 'rank',
                               'reaction_id', 'edge_type', 'role_sbo', 'location', 'location_putative', 'form']

def gene_columns(species):
    ''' The gene columns of an entity node file: the genes in the species exported for, or with no
    species filter, the genes of each species '''
    return ['genes'] if species else [f'{code}_homologues' for code in pss_schema_config.species]


REACTION_GRAPH_NODE_COLUMNS = ENTITY_NODE_COLUMNS + [
    'reaction_type', 'reaction_effect', 'reaction_mechanism', 'reaction_sbo', 'reaction_mechanism_sbo',
    'evidence_sentence']


def node(nodes, name):
    return nodes.get(name) or Node(name)


def participant_groups(reaction):
    ''' The participants (as collected, without conditions) by the names used in the rules '''
    participants = [p for p in reaction.participants if p.form != 'condition']
    return {
        'substrates': [p for p in participants if p.edge_type in edge_types.INPUTS],
        'products': [p for p in participants if p.edge_type in edge_types.OUTPUTS],
        'modifiers': [p for p in participants if p.edge_type in edge_types.MODIFIERS],
        'activating_modifiers': [p for p in participants if p.edge_type == 'ACTIVATES'],
        'inhibiting_modifiers': [p for p in participants if p.edge_type == 'INHIBITS'],
    }


#-------------------------------------
# Interaction network
#-------------------------------------

def interaction_edges(reaction, nodes):
    ''' The entity-level interaction edges of a reaction, by the interaction rules: no
    self-loops (except a modifier on its own entity), and of two edges between the same pair the
    one to the product. '''

    rules = pss_export_config.interaction_rules.get(reaction.reaction_type, {}).get(reaction.reaction_effect)
    if rules is None:
        print(f"Interaction network: {reaction.reaction_id}, no interaction rules for "
              f"({reaction.reaction_type}, {reaction.reaction_effect})")
        return []

    groups = participant_groups(reaction)
    edges = {}
    for source_group, target_group, influence, directed in rules:
        for source, target in pairs(groups[source_group], groups[target_group]):

            # no self-loops between the reaction's own inputs and outputs (an inactive -> active form,
            # a gene template -> its product, a translocation's from -> to); a modifier's influence on its
            # own entity (autoregulation, e.g. a TF activating its own transcription) is kept
            if source.name == target.name and source_group not in MODIFIER_GROUPS:
                continue

            edge = {
                'source': source.name, 'target': target.name,
                'interaction': pss_export_config.influences[influence]['interaction'],
                'directed': directed == 'directed',
                'rank': PSS_RANK,
                'source_entity': source.name, 'target_entity': target.name,
                'source_type': node(nodes, source.name).type, 'target_type': node(nodes, target.name).type,
                'source_role': participant_roles.of(reaction, source),
                'target_role': participant_roles.of(reaction, target),
                'source_location': source.location, 'target_location': target.location,
                'source_location_putative': source.location_putative,
                'target_location_putative': target.location_putative,
                'source_form': source.form, 'target_form': target.form,
                'reaction_id': reaction.reaction_id,
                'reaction_type': reaction.reaction_type,
                'reaction_effect': reaction.reaction_effect,
                'reaction_mechanism': reaction.reaction_mechanism,
                'evidence_sentence': reaction.evidence_sentence,
                'experimental_techniques': reaction.experimental_techniques,
                'reaction_sbo': sbo(reaction.reaction_type_sbo),
                'reaction_mechanism_sbo': sbo(reaction.reaction_mechanism_sbo),
                'influence_sbo': sbo(pss_export_config.influences[influence]['sbo']),
            }

            # net effect: one edge per pair of entities, the one to the product if there are several
            pair = (source.name, target.name)
            if pair in edges and edges[pair]['target_role'] == reaction.product_role:
                continue
            edges[pair] = edge

    return list(edges.values())


def entity_node_rows(names, nodes, species):
    ''' A row per entity, with its genes in the species (or, with no species, per species) '''
    rows = []
    for name in sorted(names):
        n = node(nodes, name)
        row = {'id': name, 'node_type': n.type, 'display_label': n.display_label, 'short_name': n.short_name,
               'synonyms': n.synonyms,
               'description': n.description, 'additional_information': n.additional_information,
               'pathway': n.pathway, 'all_pathways': n.all_pathways, 'mapman': n.mapman, 'functional_cluster_id': n.functional_cluster_id,
               'components': n.components,
               # the functional cluster ids of the components that are gene clusters (to match them to genes;
               # abstract clusters, like metabolites, are matched by node id)
               'component_cluster_ids': [node(nodes, c).functional_cluster_id for c in n.components
                                         if node(nodes, c).is_gene_cluster],
               'external_links': n.external_links}
        if species:
            row['genes'] = n.genes(species)
        else:
            row.update({f'{code}_homologues': n.genes(code) for code in pss_schema_config.species})
        rows.append(row)
    return rows


def create_interaction_network(pss_adapter, edges_file=None, nodes_file=None):
    ''' Write the interaction network (entity level). Returns the number of edges. '''

    edges = []
    for reaction_id in pss_adapter.reaction_ids:
        edges += interaction_edges(pss_adapter.reactions[reaction_id], pss_adapter.nodes)

    if edges_file:
        write_tsv(edges_file, INTERACTION_EDGE_COLUMNS, edges)
    if nodes_file:
        names = {e['source'] for e in edges} | {e['target'] for e in edges}
        write_tsv(nodes_file, ENTITY_NODE_COLUMNS + gene_columns(pss_adapter.species),
                  entity_node_rows(names, pss_adapter.nodes, pss_adapter.species))
    return len(edges)


#-------------------------------------
# Gene network
#-------------------------------------

def create_gene_network(pss_adapter, edges_file=None, nodes_file=None):
    ''' Write the gene network: the interaction network with functional clusters expanded into
    their genes in the species collected for (pss_adapter.species; the collector keeps only the
    reactions whose functional clusters all have genes there). Needs a species (genes of different
    species would otherwise be paired). Returns the number of edges. '''

    nodes = pss_adapter.nodes
    species = pss_adapter.species
    if species is None:
        raise ValueError("A gene network needs a species: collect_reactions(species=...)")

    def genes(name, clusters):
        return node(nodes, name).genes(species) if name in clusters else [name]

    edges = []
    gene_clusters = set()           # the clusters expanded into genes
    for reaction_id in pss_adapter.reaction_ids:
        reaction = pss_adapter.reactions[reaction_id]
        clusters = reaction.gene_clusters(nodes)

        for edge in interaction_edges(reaction, nodes):
            sources, targets = genes(edge['source'], clusters), genes(edge['target'], clusters)
            if edge['source'] == edge['target']:
                # autoregulation: each gene -> itself (not the cluster's other genes)
                gene_pairs = [(gene, gene) for gene in sources]
            else:
                # a gene in both clusters: no self-loop
                gene_pairs = [(source, target) for source, target in pairs(sources, targets) if source != target]
            # a gene's type is its cluster's class (PlantCoding, PlantNonCoding), as for CKN genes
            gene_clusters.update(name for name in (edge['source'], edge['target']) if name in clusters)
            for source, target in gene_pairs:
                edges.append({**edge, 'source': source, 'target': target})

    if edges_file:
        write_tsv(edges_file, INTERACTION_EDGE_COLUMNS, edges)
    if nodes_file:
        write_tsv(nodes_file, GENE_NODE_COLUMNS, gene_node_rows(edges, nodes, species, gene_clusters))
    return len(edges)


def gene_node_rows(edges, nodes, species, expanded):
    ''' A row per gene (from its functional clusters only: no gene-level annotations in PSS)
    and per other node. `expanded`: the functional clusters expanded into genes. A gene's
    node_type is the class of its clusters; genes in clusters of different classes are an error. '''

    gene_ids, clusters, others = set(), set(), set()
    for edge in edges:
        for side in ('source', 'target'):
            if edge[f'{side}_entity'] in expanded:
                gene_ids.add(edge[side])
                clusters.add(edge[f'{side}_entity'])
            else:
                others.add(edge[side])

    # all of a gene's functional clusters in the network (also those whose edges to it were self-loops)
    gene_clusters = {gene: [] for gene in gene_ids}
    for cluster in sorted(clusters):
        for gene in node(nodes, cluster).genes(species):
            if gene in gene_clusters:
                gene_clusters[gene].append(node(nodes, cluster))

    rows = []
    for gene in sorted(gene_clusters):
        cluster_nodes = gene_clusters[gene]
        types = {c.type for c in cluster_nodes}
        if len(types) != 1:
            raise ValueError(f"gene {gene} is in functional clusters of different classes: "
                             + ', '.join(f'{c.name} ({c.type})' for c in cluster_nodes))
        rows.append({
            'id': gene, 'node_type': types.pop(), 'species': species,
            # one entry per cluster, in the same order (empty where a cluster has no value)
            'display_label': [c.display_label for c in cluster_nodes],
            'short_name': [c.short_name for c in cluster_nodes],
            'pathway': [c.pathway for c in cluster_nodes],
            'functional_cluster_id': [c.functional_cluster_id for c in cluster_nodes],
            # all of the clusters' values
            'synonyms': list(dict.fromkeys(s for c in cluster_nodes for s in c.synonyms)),
            'all_pathways': sorted({p for c in cluster_nodes for p in c.all_pathways}),
            'mapman': sorted({b for c in cluster_nodes for b in c.mapman}),
        })
    rows += entity_node_rows(others, nodes, species)
    return rows


#-------------------------------------
# Reaction graph
#-------------------------------------

def create_reaction_graph(pss_adapter, edges_file=None, nodes_file=None):
    ''' Write the reaction graph: entities and reactions, one edge per participant (as
    collected, conditions and gene templates included). Returns the number of edges. '''

    edges = []
    for reaction_id in pss_adapter.reaction_ids:
        reaction = pss_adapter.reactions[reaction_id]
        for participant in reaction.participants:
            role = participant_roles.of(reaction, participant)
            output = participant.edge_type in edge_types.OUTPUTS
            edges.append({
                'source': reaction_id if output else participant.name,
                'role': role,
                'target': participant.name if output else reaction_id,
                'directed': True,
                'rank': PSS_RANK,
                'reaction_id': reaction_id,
                'edge_type': participant.edge_type,
                'role_sbo': sbo(pss_export_config.node_role_to_SBO.get(role)),
                'location': participant.location,
                'location_putative': participant.location_putative,
                'form': participant.form,
            })

    if edges_file:
        write_tsv(edges_file, REACTION_GRAPH_EDGE_COLUMNS, edges)
    if nodes_file:
        entities = {p.name for r in pss_adapter.reaction_ids for p in pss_adapter.reactions[r].participants}
        rows = entity_node_rows(entities, pss_adapter.nodes, pss_adapter.species)
        for reaction_id in pss_adapter.reaction_ids:
            reaction = pss_adapter.reactions[reaction_id]
            rows.append({
                'id': reaction_id, 'node_type': 'reaction', 'display_label': reaction_id,
                'external_links': reaction.external_links,
                'reaction_type': reaction.reaction_type,
                'reaction_effect': reaction.reaction_effect,
                'reaction_mechanism': reaction.reaction_mechanism,
                'reaction_sbo': sbo(reaction.reaction_type_sbo),
                'reaction_mechanism_sbo': sbo(reaction.reaction_mechanism_sbo),
                'evidence_sentence': reaction.evidence_sentence,
            })
        write_tsv(nodes_file, REACTION_GRAPH_NODE_COLUMNS + gene_columns(pss_adapter.species), rows)
    return len(edges)
