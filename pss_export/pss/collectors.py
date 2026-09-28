from ..entity_classes import Reaction
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
                 include_genes=False,
                 nodes_to_ignore='default'):

        self.graph_db = graph_db

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
            # ignore nodes that are not reactions or families
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
        cy_filters = []
        arguments = {}

        if self.reaction_filter:
            arguments['reaction_ids'] = self.REACTIONS
            cy_filters.append("r.reaction_id IN $reaction_ids")
        if self.pathway_filter:
            arguments['pathways'] = self.PATHWAYS
            cy_filters.append(
                "size(apoc.coll.intersection(n.all_pathways, $pathways)) > 0")
        # (nodes_to_ignore is applied per reaction, see _build_reaction)
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
                OPTIONAL MATCH p=(r)-[]-(n)
                {where_clause}
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

            reaction = self._build_reaction(
                reaction_id, reaction_properties, reaction_paths)
            if reaction is None:
                continue

            reactions[reaction_id] = reaction

        return reactions

    def _build_reaction(self, reaction_id, reaction_properties, reaction_paths):
        ''' Build a Reaction from its paths, leaving out edges to ignored nodes.

        The reaction itself is dropped (returns None) only if ignoring nodes
        removed all of its edges, all of its substrates, or all of its products. Removed
        modifiers, or sides which were already empty (e.g. genes not
        included), never drop the reaction.
        '''

        def _make(paths):
            reaction = Reaction(reaction_id,
                                reaction_properties['reaction_type'],
                                reaction_properties,
                                include_genes=self.include_genes)
            reaction.add_edges(paths)
            return reaction

        # paths are (r)-[]-(n), so the participant is the end node
        kept_paths = [p for p in reaction_paths
                      if p.end_node['name'] not in self.nodes_to_ignore]

        if not kept_paths:
            print(f"Reaction {reaction_id} dropped: all edges are to ignored nodes")
            return None

        reaction = _make(kept_paths)

        if len(kept_paths) < len(reaction_paths):
            full_reaction = _make(reaction_paths)
            if ((full_reaction.substrates and not reaction.substrates) or
                    (full_reaction.products and not reaction.products)):
                print(f"Reaction {reaction_id} dropped: all substrates or products are ignored nodes")
                return None

        return reaction

    def collect_node_annotations(self):

        # fetch annotations
        def _collect_node_annotations(tx):
            cy = '''
                MATCH (n)
                WHERE NOT ('Reaction' IN labels(n) OR 'Family' in labels(n) )
                RETURN n.name AS name,
                       n.description AS description,
                       n.additional_information AS additional_information,
                       n.pathway AS pathway,
                       n.external_links AS external_links,
                       n.functional_cluster_id AS functional_cluster_id,
                       n.ath_homologues AS ath_homologues,
                       labels(n) AS labels
                '''
            result = tx.run(cy)
            return [x for x in result]

        node_annotations = self.graph_db.run_query(
            _collect_node_annotations)
        return {d["name"]: dict(d) for d in node_annotations}

    def collect_reaction_pathways(self, reaction_ids):

        def _collect_reaction_pathways(tx, reaction_ids):
            cy = '''
                MATCH (r:Reaction)--(n)
                WHERE r.reaction_id IN $reaction_ids
                RETURN r.reaction_id AS reaction_id, collect(DISTINCT n.pathway) AS pathway
                '''
            result = tx.run(cy, reaction_ids=reaction_ids)
            return {r["reaction_id"]: r["pathway"] for r in result}

        return self.graph_db.run_query(
            _collect_reaction_pathways, reaction_ids)
