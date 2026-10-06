from ..entity_classes import Reaction, Node
from .config import pss_export_config, pss_schema_config

INVENTED_REASON_ALLOWLIST = ["invented:harmonise-location"]

# accepted access values -> internal access level
# ('all' is kept as a legacy alias for 'restricted')
ACCESS_LEVELS = {
    'public': 'public',
    'restricted': 'restricted',
    'all': 'restricted',
}


class PSSCollector:
    """All the logic for deciding which reactions to gather from PSS"""

    def __init__(self,
                 graph_db,
                 reactions=None,
                 access='public',
                 pathways=None,
                 species=None,
                 include_genes=False,
                 nodes_to_ignore='default'):

        self.graph_db = graph_db

        # species: only the reactions whose functional clusters all have genes in it (their
        # <species>_homologues lists); None: no species filter
        if species is not None and species not in pss_schema_config.species:
            raise ValueError(
                f"Unknown species '{species}', one of: {', '.join(pss_schema_config.species)} (or None)")
        self.species = species

        if access not in ACCESS_LEVELS:
            raise ValueError(
                f"Invalid access '{access}', must be one of: {', '.join(ACCESS_LEVELS)}")
        self.access = ACCESS_LEVELS[access]

        if reactions is not None:
            self.REACTIONS = reactions
            self.reaction_filter = True
        else:
            self.reaction_filter = False
            self.REACTIONS = None

        if (pathways is not None) and (self.reaction_filter is False):
            self.PATHWAYS = pathways
            self.pathway_filter = True
        else:
            # pathways from reactions
            self.PATHWAYS = [
                f"{x} - {y}" for d in pss_schema_config.pathways
                for x in d for y in d[x]
            ]
            self.pathway_filter = False

        self.nodes_to_ignore = self._resolve_nodes_to_ignore(nodes_to_ignore)

        self.include_genes = include_genes

    def _resolve_nodes_to_ignore(self, nodes_to_ignore):
        ''' Resolve nodes to ignore during export. '''

        if nodes_to_ignore == 'default':
            # the configured list (pss_export_config.yaml)
            return pss_export_config.nodes_to_ignore
        elif nodes_to_ignore is None:
            # no nodes to ignore
            return []
        elif isinstance(nodes_to_ignore, str):
            # assume a single node
            return [nodes_to_ignore.strip()]
        elif not isinstance(nodes_to_ignore, list):
            raise ValueError(
                "nodes_to_ignore must be a list or a comma-separated string.")
        else:
            # assume a list of nodes
            return [n.strip() for n in nodes_to_ignore]

    def _build_where_clause(self):
        ''' Filters on the reaction (r) only: selected reactions are always
        collected whole, with all of their edges. '''
        cy_filters = []
        # the nodes to ignore (models only; empty for the networks) are as if they weren't there: they don't
        # select a reaction (pathways) or leave it out (species)
        arguments = {'ignore': self.nodes_to_ignore}

        if self.reaction_filter:
            arguments['reaction_ids'] = self.REACTIONS
            cy_filters.append("r.reaction_id IN $reaction_ids")
        if self.pathway_filter:
            # reactions with at least one participant in the pathways
            arguments['pathways'] = self.PATHWAYS
            cy_filters.append('''
                EXISTS {
                    MATCH (r)--(m)
                    WHERE NOT m.name IN $ignore
                      AND size(apoc.coll.intersection(m.all_pathways, $pathways)) > 0
                }
                ''')
        if self.species is not None:
            # all gene clusters of the reaction, and those among the components of its complexes, have
            # genes in the species (reactions without them are kept; abstract clusters (PlantAbstract)
            # have no genes and, like metabolites, don't decide)
            arguments['homologues_key'] = f"{self.species}_homologues"
            cy_filters.append('''
                all(fc IN [(r)--(n:FunctionalCluster) WHERE NOT n:PlantAbstract AND NOT n.name IN $ignore | n]
                          + apoc.coll.flatten([(r)--(c:Complex) WHERE NOT c.name IN $ignore |
                              [(x:FunctionalCluster)-[:COMPONENT_OF]->(c)
                               WHERE NOT x:PlantAbstract AND NOT x.name IN $ignore | x]])
                    WHERE size(coalesce(fc[$homologues_key], [])) > 0)
                ''')
        # (nodes_to_ignore only applies to the models, see _build_reaction)
        if self.access == 'public':
            # public = at least one source which is not 'other' or 'invented',
            #          or an allowlisted 'invented' reason
            cy_filters.append('''
                (
                    size([link IN r.external_links WHERE NOT (link =~ 'other:.*' OR link =~ 'invented:.*') | 1]) > 0
                    OR size([link IN r.external_links WHERE link IN $invented_reason_allowlist | 1]) > 0
                )
                ''')
            arguments['invented_reason_allowlist'] = INVENTED_REASON_ALLOWLIST

        if len(cy_filters) == 0:
            return "", arguments

        return "WHERE " + " AND ".join(cy_filters), arguments

    def collect_reactions(self):
        '''Collect reaction list and pathway annotations (all reusable between formats)
            Limit to reactions in the PATHWAYS attr
            Limit to reactions in REACTIONS attr
        '''

        where_clause, arguments = self._build_where_clause()

        def _collect_reactions(tx, where_clause, arguments):

            cy = f'''
                MATCH (r:Reaction)
                {where_clause}
                OPTIONAL MATCH p=(r)-[]-(n)
                RETURN  r.reaction_id AS reaction_id,
                        r AS reaction,
                        collect(p) AS path
                '''
            result = tx.run(cy, **arguments)
            return list(result)

        reaction_data = self.graph_db.run_query(
            _collect_reactions, where_clause, arguments)

        reactions = {}

        for reaction_dict in reaction_data:

            reaction_id = reaction_dict['reaction_id']
            reaction_paths = reaction_dict['path']

            if len(reaction_paths) == 0:
                print(f"No edges on reaction {reaction_id}")
                continue

            reaction_properties = reaction_dict['reaction']

            reactions[reaction_id] = self._build_reaction(
                reaction_id, reaction_properties, reaction_paths)

        return reactions

    def _build_reaction(self, reaction_id, reaction_properties, reaction_paths):
        ''' Build a Reaction from its paths. The nodes to ignore are left out of the models only
        (Reaction.ignore_in_model); the networks keep every participant. '''
        reaction = Reaction(reaction_id,
                            reaction_properties['reaction_type'],
                            reaction_properties,
                            include_genes=self.include_genes)
        reaction.add_edges(reaction_paths)
        if self.nodes_to_ignore and (reason := reaction.ignore_in_model(self.nodes_to_ignore)):
            print(f"Reaction {reaction_id} left out of the models: {reason}")
        return reaction

    def collect_nodes(self):
        ''' The entities (all but reactions and families), with their annotations: {name: Node} '''

        def _collect_nodes(tx):
            cy = '''
                MATCH (n)
                WHERE NOT ('Reaction' IN labels(n) OR 'Family' in labels(n) )
                RETURN n.name AS name,
                       labels(n) AS labels,
                       n.short_name AS short_name,
                       n.display_label AS display_label,
                       n.synonyms AS synonyms,
                       n.description AS description,
                       n.additional_information AS additional_information,
                       n.pathway AS pathway,
                       n.all_pathways AS all_pathways,
                       n.mapman AS mapman,
                       n.external_links AS external_links,
                       n.functional_cluster_id AS functional_cluster_id,
                       apoc.coll.sort([(x)-[:COMPONENT_OF]->(n) | x.name]) AS components,
                       [k IN keys(n) WHERE k ENDS WITH '_homologues'
                                     AND NOT k IN ['all_homologues', '_all_homologues'] | [k, n[k]]] AS homologues
                '''
            return [record.data() for record in tx.run(cy)]

        nodes = {}
        for record in self.graph_db.run_query(_collect_nodes):
            # {species: [gene ids]}, e.g. {"ath": ["AT1G64280"], "stu": [...]}
            record["homologues"] = {k[:-len("_homologues")]: list(v or []) for k, v in record["homologues"]}
            nodes[record["name"]] = Node(**record)
        return nodes

    def collect_reaction_pathways(self, reaction_ids):

        def _collect_reaction_pathways(tx, reaction_ids):
            cy = '''
                MATCH (r:Reaction)--(n)
                WHERE r.reaction_id IN $reaction_ids
                RETURN r.reaction_id AS reaction_id, collect(DISTINCT n.pathway) AS pathway
                '''
            result = tx.run(cy, reaction_ids=reaction_ids)
            return [(r["reaction_id"], r["pathway"]) for r in result]

        # {reaction id: [the pathways of its participants]}
        return dict(self.graph_db.run_query(_collect_reaction_pathways, reaction_ids))
