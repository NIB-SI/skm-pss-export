'''
For analysing connectivity of PSS
'''

from collections import defaultdict

import networkx as nx

from ..entity_classes import IDTracker

class Graph:

    def __init__(self, pss_adapter, location=True):
        self.pss_adapter = pss_adapter
        self.location = location
        self.graph = self._create_digraph()

    def _create_digraph(self):
        idtracker = IDTracker(location=self.location)
        G = nx.DiGraph()
        for reaction in self.pss_adapter.reactions.values():
            if not reaction.in_model:
                continue
            G.add_node(reaction.id,
                       label=f"{reaction.id}\n{reaction.reaction_type}",
                       type='reaction',
                       form='reaction',
                       reaction_type=reaction.reaction_type)
            for substrate in reaction.substrates:
                substrate_id, status = idtracker.get_species_id(substrate)
                if status == 0:
                    G.add_node(substrate_id,
                               label=substrate_id,
                               type='species',
                               name=substrate.name,
                               form=substrate.form,
                               compartment=substrate.compartment)
                    idtracker.set_species_id(substrate, substrate_id)
                G.add_edge(substrate_id, reaction.id, type='substrate')

            for product in reaction.products:
                product_id, status = idtracker.get_species_id(product)
                if status == 0:
                    G.add_node(product_id,
                               label=product_id,
                               type='species',
                               name=product.name,
                               form=product.form,
                               compartment=product.compartment)
                    idtracker.set_species_id(product, product_id)
                G.add_edge(reaction.id, product_id, type='product')

            for modifier in reaction.modifiers:
                modifier_id, status = idtracker.get_species_id(modifier)
                if status == 0:
                    G.add_node(modifier_id,
                               label=modifier_id,
                               type='species',
                               name=modifier.name,
                               form=modifier.form,
                               compartment=modifier.compartment)
                    idtracker.set_species_id(modifier, modifier_id)
                G.add_edge(modifier_id, reaction.id, type='modifier')

        return G

    def get_species_per_node(self):
        species_dict = defaultdict(list)
        for node, data in self.graph.nodes(data=True):
            if data.get('type') == 'species':
                name = data.get('name')
                species_dict[name].append(node)
        return species_dict

    def is_node_connected(self, species):
        reactions = self._neighbourhood(species)
        # all_reactions_species = reactions + self._neighbourhood(reactions)
        all_reactions_species = reactions + species
        subgraph = self.graph.subgraph(all_reactions_species)

        # TODO this will not actually find all the problematic nodes
        # (e.g. there is a path from A--> A' through some other nodes... see DELLA)
        return nx.is_weakly_connected(subgraph), subgraph

    def _neighbourhood(self, nodes):
        neighbors = set()
        for n in nodes:
            neighbors.update(
                self.graph.to_undirected(as_view=True).neighbors(n))
        return list(neighbors)

    def find_problematic_nodes(self):
        species_dict = self.get_species_per_node()
        multi_species_nodes = {
            k: v
            for k, v in species_dict.items() if len(v) > 1
        }
        problematic_nodes = {}
        for node, species in multi_species_nodes.items():
            connected, subgraph = self.is_node_connected(species)
            if not connected:
                problematic_nodes[node] = {
                    "subgraph": subgraph,
                    "species": species
                }
        return problematic_nodes

class GraphVisualizer:

    @staticmethod
    def plot_subgraph(subgraph, species):
        nx.draw(
            subgraph,
            with_labels=True,
            font_weight='bold',
            # use "label" attribute for labels
            labels=nx.get_node_attributes(subgraph, 'label'),
            # make the reaction nodes smaller
            node_size=[
                400 if subgraph.nodes[n].get('type') == 'reaction' else 900
                for n in subgraph.nodes()
            ],
            # make the species of interest red, reactions yellow, the rest blue
            node_color=[
                'red' if n in species
                else 'yellow' if subgraph.nodes[n].get('type') == 'reaction'
                else 'blue'
                for n in subgraph.nodes()
            ],
            # larger arrow heads
            arrowsize=20,
            # arrow style for substrate/product/modifier
            arrowstyle=[
                '-|>' if subgraph.edges[e].get('type') == 'modifier'
                else '->'
                for e in subgraph.edges()
            ],
            # layout
            pos=nx.nx_agraph.graphviz_layout(subgraph, prog="dot"))

        from .model_fixes import _plt
        _plt().show()