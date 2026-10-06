'''
FAIDARE: a data discovery file (JSON) for the FAIDARE portal (https://urgi.versailles.inrae.fr/faidare/),
https://urgi.versailles.inrae.fr/faidare/join (Datadiscovery files).

One entry per gene and species: the genes of the gene clusters (their homologue lists) that take part in a
collected reaction. Each entry describes the gene's functional cluster (description, genes, pathway, its reactions
and their participants, synonyms, links), links to the cluster in the PSS Explorer, and has the cluster's MapMan
bins (GoMapMan 2, MapMan4) as annotations. Made from all species (collect with species=None); the access level
of the collection decides which reactions are described (public: public reactions only).
'''
import json
from collections import defaultdict

from ..entity_classes import Node
from ..pss.config import pss_schema_config
from ..pss.pss_reaction_definitions import edge_types

URL = 'https://skm.nib.si/pss/?functional_cluster_id={}'


def node(nodes, name):
    return nodes.get(name) or Node(name)


def cluster_reactions(pss_adapter):
    ''' {gene cluster name: {reaction type: {interactor labels}}}: per cluster, the types of its reactions and the
    reactions' other participants (not conditions, not the complex a reaction makes, not the cluster itself) '''
    nodes, found = pss_adapter.nodes, defaultdict(lambda: defaultdict(set))
    for reaction_id in pss_adapter.reaction_ids:
        reaction = pss_adapter.reactions[reaction_id]
        for name in reaction.gene_clusters(nodes):
            for p in reaction.participants:
                other = node(nodes, p.name)
                if p.name == name or p.form == 'condition' or \
                        (p.edge_type in edge_types.OUTPUTS and other.type == 'Complex'):
                    continue
                found[name][reaction.reaction_type].add(other.short_name or other.name)
    return found


def description(gene, cluster, reactions):
    genes = [g for code in pss_schema_config.species for g in cluster.genes(code)]
    parts = [f"{gene} belongs to the FunctionalCluster {cluster.short_name or 'unknown'} with description "
             f"'{cluster.description or ''}'. This FunctionalCluster includes the gene(s) {', '.join(genes)}. "
             f"In the Plant Stress Signalling model, it forms part of the '{cluster.pathway or 'unknown'}' pathway. ",
             f"{cluster.short_name} takes part in "
             + ' and '.join(f"{rtype} with {', '.join(sorted(others))}" for rtype, others in sorted(reactions.items()))
             + '. ']
    synonyms = [s for s in cluster.synonyms if s != cluster.short_name]     # the other names
    if synonyms:
        parts.append(f"Synonyms are: {', '.join(synonyms)}. ")
    if cluster.external_links:
        parts.append(f"Links are: {', '.join(cluster.external_links)}. ")
    return ''.join(parts)


def create_faidare(pss_adapter, filename):
    ''' Write the FAIDARE file; returns the number of entries. Needs a collection without a species filter. '''
    if pss_adapter.species is not None:
        raise ValueError("FAIDARE covers every species: collect_reactions(species=None)")
    reactions = cluster_reactions(pss_adapter)
    entries = []
    for code, species in pss_schema_config.species.items():
        for name in sorted(reactions):
            cluster = node(pss_adapter.nodes, name)
            for gene in cluster.genes(code):
                entries.append({
                    'name': gene,
                    'url': URL.format(cluster.functional_cluster_id),
                    'description': description(gene, cluster, reactions[name]),
                    'entryType': 'Gene',
                    'species': [species['name']],
                    'node': 'NIB',
                    'databaseName': 'SKM-PSS',
                    # MapMan bins "<code>_<full name>" as MapMan4:<code> (not GMM:, the classic MapMan codes)
                    'annotationId': [f"MapMan4:{b.split('_')[0]}" for b in cluster.mapman],
                    'annotationName': [f"{b.split('_', 1)[1]} (MapMan4:{b.split('_')[0]})" for b in cluster.mapman],
                })
    with open(filename, 'w', encoding='utf-8') as out:
        json.dump(entries, out, indent=2, ensure_ascii=False)
    return len(entries)
