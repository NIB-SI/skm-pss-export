'''Figures of the network projections, collected in one HTML document (docs/projections/index.html):

1. projections: an example reaction of each (reaction type, effect), as exported to the reaction graph,
   the interaction network and the gene network (docs/projections/types/);
2. chains: for each (reaction type, effect), the reactions its products are most often an input or
   modifier of, with an example pair of each of the most common links, in the interaction network of the
   two reactions (docs/projections/chains/, all counts in chains.tsv);
3. skm-tools: what simplify_pss, remove_and_rewire and remove_duplicated_binding_edges (skm-tools) do to
   the interaction network, on examples from the data (docs/projections/skm-tools/).

Everything comes from pss-export: the reactions are collected from the database (as for the exports), the
examples are exported with PSSAdapter.export(), and the figures draw the files written (read with
skm-tools for section 3). So after a change to the projections (interaction_rules in
pss_export_config.yaml, networks/network_api.py, pss_reaction_definitions.py), running this again shows
the new projections. Comments on the figures are kept in docs/projection_notes.yaml (by figure name),
and added to the document.

Run (with the database of .env; skm-tools importable for section 3, e.g. PYTHONPATH=../skm-tools):
    python docs/projection_figures.py                       # everything, and the document
    python docs/projection_figures.py --type catalysis      # only the catalysis figures (no document)
    python docs/projection_figures.py --reaction rx00001    # these reactions (no document)
'''
import argparse
import base64
import copy
import csv
import datetime
import html
import math
import re
import subprocess
import tempfile
import textwrap
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import yaml
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import ArrowStyle
from matplotlib.path import Path as MPath
from matplotlib.textpath import TextPath

from pss_export import PSSAdapter
from pss_export.pss.config import pss_export_config
from pss_export.pss.pss_reaction_definitions import edge_types, participant_roles, reaction_classes

DOCS = Path(__file__).parent
OUT = DOCS / "projections"
NOTES = DOCS / "projection_notes.yaml"
NETWORKS = [("reaction-graph", "Reaction graph"),
            ("interaction-network", "Interaction network"),
            ("gene-network", "Gene network")]

# Styles of the PSS Explorer (pss_viz services/web/utils.py: NODE_STYLE, EDGE_STYLE, GROUP): nodes by
# node_type, reaction graph edges by edge_type
NODE_STYLE = {
    "Complex": ("box", "#9cd6e4", "#40b0cb"),
    "PlantCoding": ("circle", "#66CDAA", "#48c39a"),
    "PlantNonCoding": ("box", "#98FB98", "#057c05"),
    "PlantAbstract": ("box", "#3cb371", "#2d8755"),
    "ForeignCoding": ("box", "#f4a460", "#ef7a17"),
    "ForeignNonCoding": ("box", "#f5deb3", "#e3a228"),
    "ForeignAbstract": ("box", "#bc8f8f", "#a66a6a"),
    "ForeignEntity": ("box", "#cd853f", "#965e27"),
    "ForeignAbiotic": ("box", "#cd3f40", "#a62b2c"),
    "Metabolite": ("box", "#fff0f5", "#ff6799"),
    "Process": ("box", "#c4bcff", "#6e5aff"),
    "Reaction": ("box", "#4169e1", "#122a73"),
    "default": ("box", "white", "#6c7881"),
}
NODE_GROUP = {"reaction": "Reaction", "MetaboliteFamily": "Metabolite"}
EDGE_STYLE = {          # (colour, arrowhead)
    "ACTIVATES": ("#008040", "circle"),
    "INHIBITS": ("#cd0000", "bar"),
    "SUBSTRATE": ("#5d5d5d", "arrow"),
    "TRANSLOCATE_FROM": ("#5d5d5d", "arrow"),
    "PRODUCT": ("#5d5d5d", "arrow"),
    "TRANSLOCATE_TO": ("#5d5d5d", "arrow"),
    "default": ("#000000", "arrow"),
}
EDGE_LABEL_COLOR = "#808080"
# Influence edges (interaction and gene network; the Explorer has none): as the SKM Cytoscape style
# (skm-tools resources/skm-styles.xml), in the Explorer's colours; directed = False dashed, with an open
# diamond at the source
INFLUENCE_STYLE = {
    "positive-influence": ("#008040", "arrow"),
    "negative-influence": ("#cd0000", "bar"),
    "unknown-influence": ("#5d5d5d", "diamond"),
}
HIGHLIGHT = "#ff8c00"   # border of the entity shared by chained reactions (or removed by rewiring)

# All figures the same size: the text width of an A4 page (17 cm), 3 to a page; the networks are scaled
# to their panels, the text (in points) is not
FIG_SIZE = (17 / 2.54, 8 / 2.54)
FONT_SIZE = 6           # node labels; edge labels one smaller
TITLE_SIZE = 7.5
PANEL_TITLE_SIZE = 6.5
EDGE_WIDTH = 1.0
ARROW_SIZE = 6
GAP = 10                # between nodes (points; the layout is in points, then scaled to the panels)
matplotlib.rcParams["svg.hashsalt"] = "projections"     # the same SVG ids on every run


# ---------------------------------------------------------------------------
# Collecting and exporting
# ---------------------------------------------------------------------------

def read_tsv(path):
    return pd.read_csv(path, sep="\t", quoting=csv.QUOTE_NONE, dtype=str, keep_default_na=False)


def export_reactions(adapter, reaction_ids, tmp, keys=tuple(k for k, _ in NETWORKS)):
    ''' The network exports (keys) of only these reactions (all: None), written to tmp:
    {key: (edges file, nodes file)} '''
    one = copy.copy(adapter)
    if reaction_ids is not None:
        one.reaction_ids = list(reaction_ids)
    files = {}
    for key in keys:
        edges_file, nodes_file = Path(tmp) / f"{key}-edges.tsv", Path(tmp) / f"{key}-nodes.tsv"
        one.export(key, edges_file=edges_file, nodes_file=nodes_file)
        files[key] = (edges_file, nodes_file)
    return files


def to_graph(edges_file, nodes_file):
    ''' The export files as they are (all values strings) '''
    g = nx.MultiDiGraph()
    for row in read_tsv(nodes_file).to_dict("records"):
        g.add_node(row["id"], **row)
    for row in read_tsv(edges_file).to_dict("records"):
        g.add_edge(row["source"], row["target"], **row)
    return g


def example_reactions(adapter, only_type=None):
    ''' An example reaction per (reaction type, effect): the most kinds of participants (edge types), then
    gene clusters of at most 4 genes (to see them expanded, but small enough for the figure), then the
    fewest participants. Transcriptional / translational reactions without conditions first (they are
    common there, and not in the interaction and gene networks) '''
    pairs = {(t, e) for t, effects in pss_export_config.interaction_rules.items() for e in effects}
    by_pair = {}
    for rid in adapter.reaction_ids:
        r = adapter.reactions[rid]
        pair = (r.reaction_type, r.reaction_effect)
        pairs.add(pair)
        by_pair.setdefault(pair, []).append(r)

    def score(r):
        genes = [len(adapter.nodes[c].genes(adapter.species)) for c in r.gene_clusters(adapter.nodes)]
        conditions = (r.reaction_type in reaction_classes.TRANSCRIPTIONAL_TRANSLATIONAL
                      and any(p.form == "condition" for p in r.participants))
        return (conditions, -len({p.edge_type for p in r.participants}), sum(genes) == 0,
                max(genes, default=0) > 4, len(r.participants), abs(sum(genes) - 4), r.reaction_id)

    examples = {}
    for pair in sorted(pairs, key=lambda p: (str(p[0]), str(p[1]))):
        if only_type and pair[0] != only_type:
            continue
        candidates = by_pair.get(pair)
        examples[pair] = min(candidates, key=score).reaction_id if candidates else None
    return examples


# ---------------------------------------------------------------------------
# Layout
# ---------------------------------------------------------------------------

def gene_members(gene_edges_file):
    """ The gene network's nodes by entity (from the edges' *_entity columns): {node: entity}, {entity: [nodes]} """
    entity_of = {}
    for row in read_tsv(gene_edges_file).to_dict("records"):
        for end in ("source", "target"):
            entity_of.setdefault(row[end], row[f"{end}_entity"])
    members = {}
    for node, entity in sorted(entity_of.items()):
        members.setdefault(entity, []).append(node)
    return entity_of, members


def side(edge_type):
    """ Inputs left, outputs right, activating modifiers above, inhibiting below """
    return ("left" if edge_type in edge_types.INPUTS else "right" if edge_type in edge_types.OUTPUTS
            else "bottom" if edge_type == "INHIBITS" else "top")


def participant_sides(adapter, reaction_id, nodes=None, placed=None):
    """ {participant: side} of a reaction (several: the first), of the participants in nodes (default all;
    e.g. conditions are not in the interaction network). placed: {entity: True (on its product side) or
    False (its other side)}, e.g. for a template that is also the product """
    placed = placed or {}
    sides = {}
    for p in adapter.reactions[reaction_id].participants:
        if (nodes is None or p.name in nodes) and (
                p.name not in placed or (p.edge_type in edge_types.OUTPUTS) == placed[p.name]):
            sides.setdefault(p.name, side(p.edge_type))
    return sides


def node_box(text, node_type, bold=False):
    """ A node's shape and size (points): a circle or a box around its label """
    shape = node_style(node_type)[0]
    w, h = text_size(text, bold)
    if shape == "circle":
        return shape, math.hypot(w, h) + 3, math.hypot(w, h) + 3
    return shape, w + 8, h + 5


def halves(labels):
    """ {node: (half width, half height)} in points, of node_labels """
    half = {}
    for n, label in labels.items():
        _, w, h = node_box(*label)
        half[n] = (w / 2, h / 2)
    return half


def ring_radius(diameter, n):
    """ Radius of the circle of a cluster's n genes (0: one node, in its entity's place) """
    return 0.0 if n <= 1 else max(0.8 * diameter, 1.3 * n * diameter / (2 * math.pi))


def cluster_halves(members, gene_half, half):
    """ half (sizes of the entities) with the room taken by the genes' circles of the clusters """
    half = dict(half)
    for entity, genes in members.items():
        if len(genes) > 1:
            d = 2 * max(max(gene_half[g]) for g in genes)
            r = ring_radius(d, len(genes)) + d / 2
            half[entity] = (r, r)
    return half


def entity_positions(sides, reaction_id, half, vertical=False):
    """ Reaction in the middle, its participants on their sides ({entity: side}): inputs left, outputs right,
    activating modifiers above, inhibiting below (vertical: inputs above, outputs below, activating
    modifiers right, inhibiting left), spaced by their sizes (half: {node: (half width, half height)},
    points; the positions are in points too) """
    if vertical:      # laid out horizontally, then turned
        half = {n: (hh, hw) for n, (hw, hh) in half.items()}
    size = lambda n: half.get(n, (10, 6))
    rhw, rhh = size(reaction_id)
    pos = {reaction_id: (0.0, 0.0)}
    column = 0.0
    for where, sign in (("left", -1), ("right", 1)):
        group = [e for e, s in sides.items() if s == where]
        if not group:
            continue
        x = sign * (rhw + 3.5 * GAP + max(size(e)[0] for e in group))
        # two or more further apart, three or more on an arc (the middle ones further out), so that the
        # edges between them (e.g. binding partners) are long enough to see and don't pass behind others
        gap = 2.5 * GAP if len(group) >= 2 else GAP
        bulge = 3 * GAP if len(group) >= 3 else 0
        total = sum(2 * size(e)[1] for e in group) + gap * (len(group) - 1)
        y = total / 2
        for i, e in enumerate(group):
            out = bulge * (1 - (2 * i / (len(group) - 1) - 1) ** 2) if bulge else 0
            pos[e] = (x + sign * out, y - size(e)[1])
            y -= 2 * size(e)[1] + gap
        column = max(column, total / 2)
    for where, sign in (("top", 1), ("bottom", -1)):
        group = [e for e, s in sides.items() if s == where]
        if not group:
            continue
        total = sum(2 * size(e)[0] for e in group) + GAP * (len(group) - 1)
        # above the columns of inputs and outputs if wider than the reaction
        y = sign * ((column if total / 2 > rhw + GAP else rhh) + 3.5 * GAP + max(size(e)[1] for e in group))
        x = -total / 2
        for e in group:
            pos[e] = (x + size(e)[0], y)
            x += 2 * size(e)[0] + GAP
    return {n: (y, -x) for n, (x, y) in pos.items()} if vertical else pos


def gene_positions(members, pos, gene_half):
    """ Genes on a circle around their functional cluster's position; a single gene and other nodes at
    their entity's position """
    gpos = {}
    for entity, nodes in members.items():
        cx, cy = pos.get(entity, (0.0, 0.0))
        r = ring_radius(2 * max(max(gene_half[n]) for n in nodes), len(nodes))
        for i, node in enumerate(nodes):
            a = 2 * math.pi * i / len(nodes) + math.pi / 2
            gpos[node] = (cx + r * math.cos(a), cy + r * math.sin(a))
    return gpos


def flow_positions(graphs, half):
    ''' A layered layout of a network (or of several, e.g. before and after a change, for the same
    positions in each) from its own edges, left to right: each node one layer after the nodes with edges
    into it (nodes on a cycle, e.g. binding partners, in one layer: the condensation's topological order),
    and a node with no edges into it (e.g. an enzyme) one layer before the first it has an edge to; within a
    layer, ordered by the mean height of their predecessors (fewer crossings), spaced by their sizes (half,
    points), with room for edge labels between the layers '''
    dg = nx.DiGraph()
    for g in graphs if isinstance(graphs, (list, tuple)) else [graphs]:
        dg.add_nodes_from(g)
        dg.add_edges_from((u, v) for u, v in g.edges() if u != v)
    c = nx.condensation(dg)
    layer_of = {}
    for scc in nx.topological_sort(c):
        layer_of[scc] = max((layer_of[p] + 1 for p in c.predecessors(scc)), default=0)
    for scc in c:
        if c.in_degree(scc) == 0 and c.out_degree(scc):
            layer_of[scc] = min(layer_of[v] for v in c.successors(scc)) - 1
    lowest = min(layer_of.values(), default=0)
    layer_of = {scc: layer - lowest for scc, layer in layer_of.items()}
    layers = {}
    for scc, layer in layer_of.items():
        layers.setdefault(layer, []).extend(sorted(c.nodes[scc]["members"]))
    size = lambda n: half.get(n, (10, 6))
    pos, x = {}, 0.0
    for layer in sorted(layers):
        nodes = layers[layer]
        if pos:
            x += previous + max(size(n)[0] for n in nodes) + 7 * GAP
        previous = max(size(n)[0] for n in nodes)
        def height(n):     # the mean height of the predecessors placed (0: none)
            ys = [pos[u][1] for u in dg.predecessors(n) if u in pos]
            return sum(ys) / len(ys) if ys else 0.0
        nodes = sorted(nodes, key=lambda n: -height(n))
        total = sum(2 * size(n)[1] for n in nodes) + 3 * GAP * (len(nodes) - 1)
        y = total / 2
        for n in nodes:
            pos[n] = (x, y - size(n)[1])
            y -= 2 * size(n)[1] + 3 * GAP
    return pos


def spine_positions(g, half):
    ''' A layout of a network along its main flow (falling back to flow_positions if it has none):

    - the spine, on a horizontal line: the longest path through the nodes with edges both into and out of
      them (e.g. the metabolites of a pathway; a cycle, e.g. binding partners, counts as one step), extended
      by a node with nothing upstream (or downstream) whose only neighbour on the spine is its first (last)
      node;
    - the other nodes beside the spine, at the mean position of their neighbours on it (e.g. an enzyme between
      its substrate and its product): above it those with edges into the spine only (regulators), below it
      the others; pushed outwards where they would overlap '''
    dg = nx.DiGraph((u, v) for u, v in g.edges() if u != v)
    dg.add_nodes_from(g)
    through = [n for n in dg if dg.in_degree(n) and dg.out_degree(n)]
    if not through:
        return flow_positions(g, half)
    c = nx.condensation(dg.subgraph(through))
    spine = []
    for scc in nx.dag_longest_path(c):
        members = c.nodes[scc]["members"]
        # a cycle: the member with the most edges to the rest of the path
        spine.append(max(sorted(members), key=lambda n: sum(1 for m in nx.all_neighbors(dg, n) if m not in members)))
    neighbours = lambda n: set(nx.all_neighbors(dg, n))
    first = [n for n in dg if dg.in_degree(n) == 0 and neighbours(n) & set(spine) == {spine[0]}]
    last = [n for n in dg if dg.out_degree(n) == 0 and neighbours(n) & set(spine) == {spine[-1]}]
    if first:
        spine.insert(0, max(first, key=dg.degree))
    if last:
        spine.append(max(last, key=dg.degree))

    size = lambda n: half.get(n, (10, 6))
    pos, x = {}, 0.0
    for i, n in enumerate(spine):
        if i:
            x += size(spine[i - 1])[0] + size(n)[0] + 8 * GAP
        pos[n] = (x, 0.0)
    spine_height = max(size(n)[1] for n in spine)

    placed = {"above": [], "below": []}
    rest = [n for n in dg if n not in pos]
    while rest:
        # next: the node with the most neighbours placed
        n = max(rest, key=lambda m: (len(neighbours(m) & pos.keys()), str(m)))
        rest.remove(n)
        near = [pos[m][0] for m in neighbours(n) if m in pos]
        nx_ = sum(near) / len(near) if near else x + 8 * GAP
        into = any(dg.has_edge(n, m) for m in spine)
        out_of = any(dg.has_edge(m, n) for m in spine)
        side = "above" if into and not out_of else "below"
        sign = 1 if side == "above" else -1
        y = spine_height + 3 * GAP + size(n)[1]
        # outwards until clear of the nodes on this side
        while any(abs(nx_ - px) < size(n)[0] + size(m)[0] + GAP and abs(y - abs(py)) < size(n)[1] + size(m)[1] + GAP
                  for m, (px, py) in ((m, pos[m]) for m in placed[side])):
            y += 2 * size(n)[1] + GAP
        pos[n] = (nx_, sign * y)
        placed[side].append(n)
    return pos


# ---------------------------------------------------------------------------
# Drawing, with networkx (nodes as the Explorer's, edges with its heads)
# ---------------------------------------------------------------------------

def wrap(text, width=16):
    """ Lines of at most width characters (keeping the line breaks in text) """
    return "\n".join(line for part in str(text).split("\n")
                     for line in textwrap.wrap(part, width, break_on_hyphens=False) or [part])


def node_style(node_type):
    return NODE_STYLE.get(NODE_GROUP.get(node_type, node_type), NODE_STYLE["default"])


class Heads(ArrowStyle._Base):
    """ An arrowstyle with the Explorer / Cytoscape heads, for networkx edges (matplotlib has no circle or
    diamond heads): end at the target (arrow, circle, bar, diamond), optionally start at the source
    (open diamond). parts: "all", or only the "line" or the "heads" (a dashed edge has solid heads) """

    def __init__(self, end, start=None, parts="all"):
        self.end, self.start, self.parts = end, start, parts
        super().__init__()

    @staticmethod
    def head(kind, tip, direction, size):
        """ The head's path, whether it is filled, and where the line ends (back from the tip) """
        d = direction / np.hypot(*direction)
        n = np.array([-d[1], d[0]])
        if kind == "arrow":
            return MPath([tip, tip - size * d + 0.45 * size * n, tip - size * d - 0.45 * size * n, tip],
                         closed=True), True, 0.5 * size
        if kind == "circle":
            return MPath.circle(tip - 0.35 * size * d, 0.35 * size), True, 0.35 * size
        if kind == "bar":       # a little before the node, so its border doesn't hide it
            at = tip - 0.3 * size * d
            return MPath([at + 0.6 * size * n, at - 0.6 * size * n]), False, 0.3 * size
        # diamond, open diamond
        verts = [tip, tip - size / 2 * d + 0.4 * size * n, tip - size * d, tip - size / 2 * d - 0.4 * size * n, tip]
        return MPath(verts, closed=True), kind == "diamond", size

    def transmute(self, path, mutation_size, linewidth):
        v = path.vertices.copy()
        paths, fillable = [], []
        for kind, i, j in [(self.end, -1, -2), (self.start, 0, 1)]:
            if kind:
                direction = v[i] - v[j]
                head, fill, back = self.head(kind, v[i], direction, mutation_size)
                paths.append(head)
                fillable.append(fill)
                v[i] = v[i] - back * direction / np.hypot(*direction)     # the line ends at the head
        if self.parts == "heads":
            return paths, fillable
        if self.parts == "line":
            return [MPath(v, path.codes)], [False]
        return [MPath(v, path.codes)] + paths, [False] + fillable


def text_size(text, bold=False):
    """ Width and height of a (multi-line) label, in points """
    lines = text.split("\n")
    prop = FontProperties(size=FONT_SIZE, weight="bold" if bold else "normal")
    width = max(TextPath((0, 0), line or " ", prop=prop).get_extents().width for line in lines)
    return width, len(lines) * FONT_SIZE * 1.2


def panel_scale(ax):
    """ Points per data unit (the layout is in points), at most 1: if a layout is larger than its panel, its
    nodes, text and edges are drawn smaller with it """
    box = ax.get_position()
    x0, x1 = ax.get_xlim()
    return min(1.0, box.width * ax.figure.get_figwidth() * 72 / (x1 - x0))


def draw_nodes(ax, pos, nodes, highlight=()):
    """ nodes: {id: (label, node_type, bold)}, drawn with networkx as the Explorer's shapes and colours: a
    circle or a box around the label; highlight: nodes with a thick border.
    Returns {id: (shape, width, height)} in points """
    g = nx.Graph()
    g.add_nodes_from(nodes)
    geometry = {}
    scale = panel_scale(ax)
    for n, (text, node_type, bold) in nodes.items():
        shape, fill, border = node_style(node_type)
        shape, w, h = node_box(text, node_type, bold)
        w, h = w * scale, h * scale
        marker = "o" if shape == "circle" else MPath(
            [(-w / 2, -h / 2), (w / 2, -h / 2), (w / 2, h / 2), (-w / 2, h / 2), (-w / 2, -h / 2)],
            closed=True)        # (a marker is scaled to its longest side)
        geometry[n] = (shape, w, h)
        nx.draw_networkx_nodes(g, pos, nodelist=[n], ax=ax, node_size=max(w, h) ** 2, node_shape=marker,
                               node_color=fill, edgecolors=HIGHLIGHT if n in highlight else border,
                               linewidths=(2.2 if n in highlight else 0.8) * scale)
        nx.draw_networkx_labels(g, pos, labels={n: text}, ax=ax, font_size=FONT_SIZE * scale,
                                font_color="white" if node_type == "reaction" else "black",
                                font_weight="bold" if bold else "normal")
    return geometry


CONNECTIONS = [f"arc3,rad={0.15 + 0.15 * i}" for i in range(6)]   # parallel edges bend apart


def border(geometry, direction):
    """ Distance (points) from a node's centre to its border, in direction (data units, aspect equal) """
    shape, w, h = geometry
    if shape == "circle":
        return w / 2
    dx, dy = np.abs(direction) / np.hypot(*direction)
    return min(w / 2 / dx if dx else np.inf, h / 2 / dy if dy else np.inf)


def draw_edges(ax, g, pos, geometry, style_of, label_of=None):
    """ The edges of g with networkx, one at a time so that each ends at its nodes' borders (networkx clips
    at a circle): style_of(data) -> (colour, head); dashed with an open diamond at the source if
    directed = False; labelled by label_of(data) """
    count = Counter()
    scale = panel_scale(ax)
    edges = g.edges(keys=True, data=True) if g.is_multigraph() else ((u, v, None, d) for u, v, d in g.edges(data=True))
    for u, v, k, d in edges:
        color, kind = style_of(d)
        dashed = str(d.get("directed")) == "False"
        connection = CONNECTIONS[count[(u, v)] % len(CONNECTIONS)]    # as networkx numbers parallel edges
        count[(u, v)] += 1
        # the curve leaves u and enters v towards its arc3 control point
        p0, p1 = np.asarray(pos[u], float), np.asarray(pos[v], float)
        rad = float(connection.split("=")[1])
        control = (p0 + p1) / 2 + rad * np.array([p1[1] - p0[1], -(p1[0] - p0[0])])
        if u == v:
            sizes = {u: (geometry[u][2]) ** 2}
        else:
            sizes = {u: (2 * border(geometry[u], control - p0)) ** 2, v: (2 * border(geometry[v], control - p1)) ** 2}
        h = g.subgraph([u, v])
        edge = (u, v) if k is None else (u, v, k)
        edge_pos, edges, colors, extra = pos, [edge], [color], []
        if u == v:
            # networkx makes a self-loop 1/10 of the height of the edges drawn with it (from the node's
            # centre): an invisible spacer edge sets it to a little over the node's
            h = nx.MultiDiGraph(h) if h.is_multigraph() else nx.DiGraph(h)
            loop = 2.6 * geometry[u][2] / 2 / scale       # data units
            spacer = ("spacer:a", "spacer:b")
            h.add_edge(*spacer)
            edge_pos = {**pos, spacer[0]: pos[u], spacer[1]: (pos[u][0], pos[u][1] + 10 * loop)}
            edges, colors, extra = [edge, spacer], [color, (0, 0, 0, 0)], ["-"]
            sizes.update({spacer[0]: 0, spacer[1]: 0})
        common = dict(edgelist=edges, ax=ax, width=EDGE_WIDTH * scale, arrows=True, arrowsize=ARROW_SIZE * scale,
                      nodelist=list(sizes), node_size=list(sizes.values()), connectionstyle=connection)
        # dashed: the line, then its heads solid
        nx.draw_networkx_edges(h, edge_pos, edge_color=colors, style="dashed" if dashed else "solid",
                               arrowstyle=[Heads(kind, "open diamond" if dashed else None,
                                                 "line" if dashed else "all")] + extra, **common)
        if dashed:
            nx.draw_networkx_edges(h, edge_pos, edge_color=colors,
                                   arrowstyle=[Heads(kind, "open diamond", "heads")] + extra, **common)
        if label_of:
            nx.draw_networkx_edge_labels(h, pos, edge_labels={edge: label_of(d)}, ax=ax,
                                         font_size=(FONT_SIZE - 1) * scale, font_color=EDGE_LABEL_COLOR,
                                         connectionstyle=connection, nodelist=list(sizes),
                                         node_size=list(sizes.values()),
                                         bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.85))


def node_labels(g, label_of=None):
    """ {id: (label, node_type, bold)} as the Explorer: the display label; a reaction its id (bold; its type is
    in the title) """
    labels = {}
    for n, d in g.nodes(data=True):
        if d.get("node_type") == "reaction":
            labels[n] = (n, "reaction", True)
        else:
            labels[n] = (wrap(label_of(n) if label_of else d.get("display_label") or n, 12), d.get("node_type"), False)
    return labels


def influence_style(d):
    return INFLUENCE_STYLE[d["interaction"]]


def draw_influences(ax, g, pos, highlight=(), label_of=None):
    """ An interaction (or gene) network: nodes, and edges by their influence """
    geometry = draw_nodes(ax, pos, node_labels(g), highlight)
    draw_edges(ax, g, pos, geometry, influence_style, label_of)


def figure(panels, title, stacked=False):
    ''' A figure of FIG_SIZE with panels side by side (stacked: one above the other, for wide networks), and
    the title above them '''
    fig = plt.figure(figsize=FIG_SIZE)
    title = wrap(title, 115)
    fig.suptitle(title, fontsize=TITLE_SIZE, y=0.985)
    height = 0.86 - 0.045 * title.count("\n")
    if stacked:
        gap = 0.08      # (room for the panel titles)
        each = (height - gap * (panels - 1)) / panels
        axes = [fig.add_axes([0.01, 0.01 + (panels - 1 - i) * (each + gap), 0.98, each]) for i in range(panels)]
    else:
        gap, width = 0.01, (0.98 - 0.01 * (panels - 1)) / panels
        axes = [fig.add_axes([0.01 + i * (width + gap), 0.01, width, height]) for i in range(panels)]
    return fig, axes


def setup_panel(ax, title, pos, half, pad=(15, 15)):
    ''' A panel showing the nodes (at pos, of half sizes, plus pad; points) at the same scale in x and y, the
    panel's box fixed (the limits are widened to its shape) '''
    boxes = [(x, y, *half.get(n, (10, 6))) for n, (x, y) in pos.items()]
    x0, x1 = min(x - w for x, y, w, h in boxes) - pad[0], max(x + w for x, y, w, h in boxes) + pad[0]
    y0, y1 = min(y - h for x, y, w, h in boxes) - pad[1], max(y + h for x, y, w, h in boxes) + pad[1]
    ax.set_title(title, fontsize=PANEL_TITLE_SIZE, pad=3)
    for spine in ax.spines.values():
        spine.set_color("#dddddd")
    box = ax.get_position()
    box_ratio = box.width * ax.figure.get_figwidth() / (box.height * ax.figure.get_figheight())
    if (x1 - x0) / (y1 - y0) < box_ratio:     # widen x, or else y, to the shape of the box
        grow = ((y1 - y0) * box_ratio - (x1 - x0)) / 2
        x0, x1 = x0 - grow, x1 + grow
    else:
        grow = ((x1 - x0) / box_ratio - (y1 - y0)) / 2
        y0, y1 = y0 - grow, y1 + grow
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])


def empty_panel(ax, text):
    ax.text(0.5, 0.5, text, ha="center", va="center", transform=ax.transAxes, fontsize=PANEL_TITLE_SIZE,
            color="#7f7f7f")


def save(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path.with_suffix(".svg"), metadata={"Date": None})
    fig.savefig(path.with_suffix(".png"), dpi=200)
    plt.close(fig)
    return path.with_suffix(".svg")


def draw_legend(path):
    ''' One legend for all figures '''
    def line(color, marker, text, ls="-", mfc=None):
        return Line2D([], [], color=color, lw=1.5, ls=ls, marker=marker, markersize=5,
                      markerfacecolor=mfc or color, label=text)
    handles = [line(c, {"arrow": ">", "circle": "o", "bar": "|"}[k], f"{t.replace('_', ' ')} (reaction graph)")
               for t, (c, k) in EDGE_STYLE.items() if t in ("SUBSTRATE", "ACTIVATES", "INHIBITS")]
    handles += [line(c, {"arrow": ">", "bar": "|", "diamond": "D"}[k], t) for t, (c, k) in INFLUENCE_STYLE.items()]
    handles.append(line("#5d5d5d", "D", "both directions (undirected)", ls="--", mfc="white"))
    handles += [Line2D([], [], ls="", marker="o" if s == "circle" else "s", markersize=7, markerfacecolor=f,
                       markeredgecolor=b, label=t) for t, (s, f, b) in NODE_STYLE.items()
                if t in ("PlantCoding", "Complex", "Metabolite", "Reaction", "default")]
    handles[-1].set_label("other (e.g. condition)")
    handles.append(Line2D([], [], ls="", marker="o", markersize=7, markerfacecolor="white", markeredgecolor=HIGHLIGHT,
                          markeredgewidth=2, label="shared / removed entity"))
    fig = plt.figure(figsize=(FIG_SIZE[0], 0.75))
    fig.legend(handles=handles, loc="center", fontsize=FONT_SIZE, frameon=False, ncol=4)
    return save(fig, path)


# ---------------------------------------------------------------------------
# 1. Projections: a reaction of each type in the three networks
# ---------------------------------------------------------------------------

def rules_of(reaction_type, effect):
    return pss_export_config.interaction_rules.get(reaction_type, {}).get(effect) or []


def reaction_title(reaction):
    return (f"{reaction.reaction_id}: {reaction.reaction_type}, {reaction.reaction_effect}"
            + (f" ({reaction.reaction_mechanism})" if reaction.reaction_mechanism else ""))


def draw_reaction(adapter, reaction_id, path):
    reaction = adapter.reactions[reaction_id]
    with tempfile.TemporaryDirectory() as tmp:
        files = export_reactions(adapter, [reaction_id], tmp)
        graphs = {key: to_graph(*files[key]) for key, _ in NETWORKS}
        entity_of, members = gene_members(files["gene-network"][0])
    # genes by their id (their display_label is their clusters')
    def gene_label(n):    # an AGI on two lines (smaller circles)
        if entity_of.get(n) == n:
            return graphs["gene-network"].nodes[n].get("display_label") or n
        agi = re.fullmatch(r"(AT[0-9CM]G)(\d+)", n)
        return f"{agi[1]}\n{agi[2]}" if agi else n
    labels = node_labels(graphs["reaction-graph"])
    gene_labels = node_labels(graphs["gene-network"], gene_label)
    half, gene_half = halves(labels), halves(gene_labels)
    sides = participant_sides(adapter, reaction_id, set(graphs["reaction-graph"]))
    pos = entity_positions(sides, reaction_id, cluster_halves(members, gene_half, half), vertical=True)
    gpos = gene_positions(members, pos, gene_half)
    clusters = sorted({e for n, e in entity_of.items() if n != e})
    # the same scale in all panels: the nodes of all
    all_pos = {**{("rg", n): p for n, p in pos.items()}, **{("gn", n): p for n, p in gpos.items()}}
    all_half = {**{("rg", n): h for n, h in half.items()}, **{("gn", n): h for n, h in gene_half.items()}}

    fig, axes = figure(3, f"{reaction_title(reaction)}   [species: {adapter.species}]")
    for ax, (key, title) in zip(axes, NETWORKS):
        g = graphs[key]
        setup_panel(ax, f"{title}: {g.number_of_nodes()} nodes, {g.number_of_edges()} edges", all_pos, all_half)
        if g.number_of_edges() == 0:
            empty_panel(ax, "no edges")
            continue
        if key == "reaction-graph":
            geometry = draw_nodes(ax, pos, labels)
            draw_edges(ax, g, pos, geometry, lambda d: EDGE_STYLE.get(d["edge_type"], EDGE_STYLE["default"]),
                       lambda d: f"{d['edge_type'].replace('_', ' ')}\n({d['role']})")
        elif key == "interaction-network":
            draw_influences(ax, g, pos)
        else:
            for e in clusters:     # the functional cluster, behind its genes
                shape, fill, border_color = node_style(adapter.nodes[e].type if e in adapter.nodes else None)
                ax.text(*pos[e], wrap(adapter.nodes[e].display_label if e in adapter.nodes else e),
                        ha="center", va="center", fontsize=(FONT_SIZE + 1) * panel_scale(ax), color=border_color,
                        fontweight="bold",
                        zorder=0)
            geometry = draw_nodes(ax, gpos, gene_labels)
            draw_edges(ax, g, gpos, geometry, influence_style)
    return save(fig, path)


def autoregulation_reactions(adapter):
    ''' Reactions with a modifier that is also their product (autoregulation, e.g. a transcription factor
    activating its own transcription; the only self-loops of the networks): the first of each
    (type, effect, modifier edge type) '''
    found = {}
    for rid in adapter.reaction_ids:
        r = adapter.reactions[rid]
        products = {p.name for p in r.participants if p.edge_type in edge_types.OUTPUTS}
        for p in r.participants:
            if p.edge_type not in edge_types.INPUTS + edge_types.OUTPUTS and p.name in products:
                found.setdefault((r.reaction_type, r.reaction_effect, p.edge_type), rid)
    return found


def draw_autoregulation(adapter, out):
    records = []
    for (reaction_type, effect, edge_type), rid in sorted(autoregulation_reactions(adapter).items()):
        name = f"{file_name(reaction_type, effect)}-autoregulation"
        svg = draw_reaction(adapter, rid, out / name)
        print(f"autoregulation ({reaction_type}, {effect}, {edge_type}): {rid} -> {svg}")
        records.append({"name": name, "title": f"{reaction_type}, {effect}: autoregulation ({edge_type})",
                        "svg": svg, "rules": rules_of(reaction_type, effect),
                        "text": f"Example reaction: {rid}."})
    return records


def binding_partners_reaction(adapter, effect="activation"):
    ''' A binding with three or more interactors (substrates), to show the mutual edges between all of them:
    the fewest interactors from three, then gene clusters of at most 4 genes, then the fewest genes '''
    def score(r):
        inputs = {p.name for p in r.participants if p.edge_type in edge_types.INPUTS}
        genes = [len(adapter.nodes[c].genes(adapter.species)) for c in r.gene_clusters(adapter.nodes)]
        return (len(inputs) < 3, len(inputs), max(genes, default=0) > 4, sum(genes), r.reaction_id)
    candidates = [adapter.reactions[rid] for rid in adapter.reaction_ids
                  if adapter.reactions[rid].reaction_type == "binding/oligomerisation"
                  and adapter.reactions[rid].reaction_effect == effect]
    best = min(candidates, key=score, default=None)
    return best.reaction_id if best is not None and not score(best)[0] else None


def unmodified_translocation(adapter):
    ''' A translocation without modifiers (no transporter): the fewest participants, then genes '''
    def score(r):
        genes = [len(adapter.nodes[c].genes(adapter.species)) for c in r.gene_clusters(adapter.nodes)]
        return (len(r.participants), sum(genes), r.reaction_id)
    candidates = [r for r in (adapter.reactions[rid] for rid in adapter.reaction_ids)
                  if r.reaction_type == "translocation"
                  and not any(p.edge_type in edge_types.MODIFIERS for p in r.participants)]
    return min(candidates, key=score).reaction_id if candidates else None


def draw_unmodified_translocation(adapter, out):
    rid = unmodified_translocation(adapter)
    if rid is None:
        return []
    reaction = adapter.reactions[rid]
    name = f"{file_name(reaction.reaction_type, reaction.reaction_effect)}-no_modifiers"
    svg = draw_reaction(adapter, rid, out / name)
    print(f"translocation without modifiers: {rid} -> {svg}")
    return [{"name": name, "title": f"translocation, {reaction.reaction_effect}: without modifiers", "svg": svg,
             "rules": rules_of(reaction.reaction_type, reaction.reaction_effect),
             "text": f"Example reaction: {rid}."}]


def draw_binding_partners(adapter, out):
    rid = binding_partners_reaction(adapter)
    if rid is None:
        return []
    name = "binding_oligomerisation-activation-interactors"
    svg = draw_reaction(adapter, rid, out / name)
    print(f"binding with interactors: {rid} -> {svg}")
    return [{"name": name, "title": "binding/oligomerisation, activation: three or more interactors", "svg": svg,
             "rules": rules_of("binding/oligomerisation", "activation"),
             "text": f"Example reaction: {rid}."}]


def draw_types(adapter, out, only_type=None, reaction_ids=None):
    ''' The projection figures: an example per (type, effect), or these reactions. Returns their records
    for the document '''
    if reaction_ids:
        todo = {(adapter.reactions[r].reaction_type, adapter.reactions[r].reaction_effect, r): r for r in reaction_ids}
        names = {k: file_name(k[0], k[1]) + f"-{k[2]}" for k in todo}
    else:
        todo = example_reactions(adapter, only_type)
        names = {k: file_name(*k) for k in todo}
    records = []
    for key, reaction_id in todo.items():
        if reaction_id is None:
            print(f"{key}: no reaction of this type and effect in the collection")
            records.append({"name": names[key], "title": f"{key[0]}, {key[1]}", "rules": rules_of(*key[:2]),
                            "text": "No reaction of this type and effect in the collection."})
            continue
        svg = draw_reaction(adapter, reaction_id, out / names[key])
        print(f"{key}: {reaction_id} -> {svg}")
        records.append({"name": names[key], "title": f"{key[0]}, {key[1]}", "svg": svg,
                        "rules": rules_of(*key[:2]),
                        "text": f"Example reaction: {reaction_id}."})
    return records


# ---------------------------------------------------------------------------
# 2. Chains: a product of one reaction as an input or modifier of another
# ---------------------------------------------------------------------------

def chain_counts(adapter):
    ''' Each product of a reaction that is an input or modifier of another reaction:
    {(type, effect): {(downstream type, effect, role there): [(reaction, product, downstream reaction)]}} '''
    uses = {}
    for rid in adapter.reaction_ids:
        r = adapter.reactions[rid]
        for p in r.participants:
            if p.edge_type not in edge_types.OUTPUTS:
                uses.setdefault(p.name, []).append((r, participant_roles.of(r, p)))
    chains = {}
    for rid in adapter.reaction_ids:
        r = adapter.reactions[rid]
        for entity in sorted({p.name for p in r.participants if p.edge_type in edge_types.OUTPUTS}):
            for down, role in uses.get(entity, []):
                if down is not r:
                    chains.setdefault((r.reaction_type, r.reaction_effect), {}).setdefault(
                        (down.reaction_type, down.reaction_effect, role), []).append((rid, entity, down.reaction_id))
    return chains


def example_chain(adapter, links):
    ''' The example of the links (reaction, product, downstream reaction): one shared entity, the most kinds
    of participants (edge types), then the fewest participants '''
    def score(link):
        r1, r2 = adapter.reactions[link[0]], adapter.reactions[link[2]]
        shared = {p.name for p in r1.participants} & {p.name for p in r2.participants}
        kinds = len({p.edge_type for p in r1.participants}) + len({p.edge_type for p in r2.participants})
        return (len(shared), -kinds, len(r1.participants) + len(r2.participants), link)
    return min(links, key=score)


def reaction_label(d):
    """ An edge's reaction: its id and type (broken after a '/' too) """
    return f"{d['reaction_id']}\n" + wrap(d["reaction_type"].replace("/", "/ "), 16).replace("/ ", "/")


def label(adapter, name):
    """ A node's display label """
    return adapter.nodes[name].display_label if name in adapter.nodes else name


def chain_title(adapter, link, role):
    r1, entity, r2 = link
    a, b = adapter.reactions[r1], adapter.reactions[r2]
    return (f"{r1}: {a.reaction_type}, {a.reaction_effect}  →  {r2}: {b.reaction_type}, {b.reaction_effect}; "
            f"{adapter.nodes[entity].display_label if entity in adapter.nodes else entity} as {role}")


def draw_chain(adapter, link, role, path):
    r1, entity, r2 = link
    with tempfile.TemporaryDirectory() as tmp:
        g = to_graph(*export_reactions(adapter, [r1, r2], tmp, ["interaction-network"])["interaction-network"])
    half = halves(node_labels(g))
    pos = flow_positions(g, half)
    fig, (ax,) = figure(1, chain_title(adapter, link, role))
    setup_panel(ax, f"Interaction network of {r1} and {r2}: {g.number_of_nodes()} nodes, "
                    f"{g.number_of_edges()} edges", {n: pos[n] for n in g}, half)
    if g.number_of_edges() == 0:
        empty_panel(ax, "no edges")
    draw_influences(ax, g, pos, highlight={entity}, label_of=reaction_label)
    return save(fig, path)


def draw_chains(adapter, out, top, only_type=None):
    ''' The chain figures (the top links of each type) and chains.tsv. Returns (the records of the figures
    for the document, the table of all links, the links drawn: {figure name: (link, role)}) '''
    chains = chain_counts(adapter)
    out.mkdir(parents=True, exist_ok=True)
    rows, records, drawn = [], [], {}
    for up in sorted(chains, key=lambda p: (str(p[0]), str(p[1]))):
        ranked = sorted(chains[up].items(), key=lambda kv: (-len(kv[1]), str(kv[0])))
        for rank, (down, links) in enumerate(ranked, 1):
            link = example_chain(adapter, links)
            rows.append({"type": up[0], "effect": up[1], "downstream_type": down[0], "downstream_effect": down[1],
                         "role": down[2], "links": len(links), "example": " ".join(link)})
            if rank > top or (only_type and up[0] != only_type):
                continue
            name = f"{file_name(*up)}--{file_name(down[0], down[1])}-{down[2]}"
            svg = draw_chain(adapter, link, down[2], out / name)
            print(f"{up} -> {down} ({len(links)} links): {' '.join(link)} -> {svg}")
            drawn[name] = (link, down[2])
            r1, r2 = (adapter.reactions[r] for r in (link[0], link[2]))
            records.append({"name": name, "group": f"{up[0]}, {up[1]}",
                            "title": f"→ {down[0]}, {down[1]}, as {down[2]} ({len(links)} pairs)", "svg": svg,
                            "rules": rules_of(r1.reaction_type, r1.reaction_effect),
                            "rules2": rules_of(r2.reaction_type, r2.reaction_effect),
                            "text": f"Example: {link[0]} makes {label(adapter, link[1])}, which {link[2]} uses."})
    table = pd.DataFrame(rows)
    table.to_csv(out / "chains.tsv", sep="\t", index=False)
    return records, table, drawn


# ---------------------------------------------------------------------------
# 3. skm-tools: simplify_pss, remove_and_rewire, remove_duplicated_binding_edges
# ---------------------------------------------------------------------------

def reaction_ids_label(d):
    ids = d.get("reaction_id")
    ids = ids if isinstance(ids, list) else [ids]
    return wrap(", ".join(ids), 20) + ("\n(rewired)" if str(d.get("note") or "").startswith("rewired") else "")


def load_interaction_network(adapter, reaction_ids, tmp):
    ''' The interaction network export of these reactions (all: None), read with skm-tools '''
    from skm_tools.pss import pss_interaction_network_to_networkx
    edges_file, nodes_file = export_reactions(adapter, reaction_ids, tmp, ["interaction-network"])["interaction-network"]
    return pss_interaction_network_to_networkx(edges_file, nodes_file)


SIMPLIFY_CASES = [
    # (name, title, which node pairs (score: lowest first; None: not this case), text)
    ("same_influence", "parallel edges of the same influence",
     lambda n, directions, influences, pair: None if directions > 1 or influences > 1 or n < 2 or n > 4
     else (abs(n - 3), pair),
     "Parallel edges in the same direction and with the same influence become one edge, labelled with all "
     "their reaction ids."),
    ("different_influences", "parallel edges of different influences",
     lambda n, directions, influences, pair: None if directions > 1 or influences < 2 or n > 4
     else (-influences, -n, pair),
     "Parallel edges with different influences (e.g. one positive, one negative) become one edge of unknown "
     "influence. To keep one edge per influence instead, use <code>split_on_attrs='interaction'</code>."),
    ("both_directions", "edges in both directions",
     lambda n, directions, influences, pair: None if directions < 2 or n > 4 else (-min(influences, 3), -n, pair),
     "Edges are merged per direction, so two nodes linked both ways (e.g. binding partners) keep two edges, "
     "one each way. They are drawn as two curves, but they are not leftover parallel edges."),
]


def draw_simplify(adapter, out):
    ''' simplify_pss on the edges between two nodes, for each case of SIMPLIFY_CASES: the node pair that
    fits it best (at most 4 edges: readable) '''
    from skm_tools.pss import simplify_pss
    with tempfile.TemporaryDirectory() as tmp:
        g = load_interaction_network(adapter, None, tmp)
    pairs = {}
    for u, v, k, d in g.edges(keys=True, data=True):
        if u != v:
            pairs.setdefault(frozenset((u, v)), []).append((u, v, k))
    records = []
    for case, title, score, text in SIMPLIFY_CASES:
        scored = []
        for pair, edges in pairs.items():
            directions = len({e[:2] for e in edges})
            influences = len({g.edges[e]["interaction"] for e in edges})
            key = score(len(edges), directions, influences, sorted(pair))
            if key is not None:
                scored.append((key, pair, edges))
        if not scored:
            continue
        _, pair, edges = min(scored)
        sub = g.edge_subgraph(edges).copy()
        u, v = sorted(pair)
        pos = {u: (-70.0, 0.0), v: (70.0, 0.0)}
        half = halves(node_labels(sub))
        panels = [("Interaction network", sub),
                  ("simplify_pss(g)", simplify_pss(sub)),
                  ("split_on_attrs='interaction'", simplify_pss(sub, split_on_attrs="interaction"))]
        fig, axes = figure(3, f"simplify_pss, {title}: {label(adapter, u)} and {label(adapter, v)}")
        for ax, (panel_title, h) in zip(axes, panels):
            setup_panel(ax, f"{panel_title}: {h.number_of_edges()} edges", pos, half, pad=(15, 40))
            draw_influences(ax, h, pos, label_of=reaction_ids_label)
        name = f"simplify_pss-{case}"
        svg = save(fig, out / name)
        print(f"simplify_pss ({case}): {u}, {v} -> {svg}")
        records.append({"name": name, "title": f"simplify_pss: {title}", "svg": svg, "html": True,
                        "text": f"{text} Example: the {len(edges)} edges between {html.escape(label(adapter, u))} "
                                f"and {html.escape(label(adapter, v))}."})
    if records:
        records[0]["text"] = ("The interaction network has one edge per reaction, so two nodes can be linked by "
                              "several parallel edges. simplify_pss merges them into one edge per "
                              "direction. " + records[0]["text"])
    return records


def draw_rewire(adapter, drawn, out):
    ''' remove_and_rewire of the shared entity of chain examples: the first (in the order of the chain
    figures) whose rewiring gives a positive edge, a negative edge, an unknown one, or nothing (removing a
    binding partner) '''
    from skm_tools.pss import remove_and_rewire, simplify_pss
    found = {}
    for name, (link, role) in drawn.items():
        r1, entity, r2 = link
        with tempfile.TemporaryDirectory() as tmp:
            g = simplify_pss(load_interaction_network(adapter, [r1, r2], tmp))
        if entity not in g:
            continue
        added = remove_and_rewire(g, entity, dry_run=True)
        kind = ("nothing rewired" if not added else
                "negative" if any(d["interaction"] == "negative-influence" for *_, d in added) else
                "positive" if all(d["interaction"] == "positive-influence" for *_, d in added) else "unknown")
        # nothing rewired: through a binding partner (the mutual binding edges are not rewired through), not
        # e.g. an entity with nothing upstream
        if kind != "nothing rewired" or role == participant_roles.INTERACTOR:
            found.setdefault(kind, (name, link, role, g))
    records = []
    for kind in ("positive", "negative", "unknown", "nothing rewired"):
        if kind not in found:
            continue
        name, link, role, g = found[kind]
        r1, entity, r2 = link
        after = g.copy()
        reason = remove_and_rewire(after, entity)[entity]
        half = halves(node_labels(g))
        pos = flow_positions([g, after], half)
        fig, axes = figure(2, f"remove_and_rewire of {g.nodes[entity].get('display_label') or entity}: "
                              f"{chain_title(adapter, link, role)}", stacked=True)
        for ax, (title, h, highlight) in zip(axes, [("simplify_pss(g)", g, {entity}),
                                                     (f"remove_and_rewire(g, '{entity}'): {reason}", after, ())]):
            setup_panel(ax, f"{title}; {h.number_of_edges()} edges", {n: pos[n] for n in g}, half)
            draw_influences(ax, h, pos, highlight, label_of=reaction_ids_label)
        file = f"remove_and_rewire-{kind.replace(' ', '_')}"
        svg = save(fig, out / file)
        print(f"remove_and_rewire ({kind}): {' '.join(link)} -> {svg}")
        records.append({"name": file, "title": f"remove_and_rewire: {kind}", "svg": svg,
                        "text": f"Removing {label(adapter, entity)}, which links reactions {r1} and {r2}."})
    if records:
        records[0]["text"] = ("remove_and_rewire removes a node and connects each node upstream of it to each "
                              "node downstream, with the combined influence (e.g. activating an inhibitor gives a "
                              "negative edge). The network is simplified first, as remove_and_rewire needs one "
                              "edge per direction. " + records[0]["text"])
    return records


METABOLITE_TYPES = ("Metabolite", "MetaboliteFamily")


def metabolite_pathway(adapter, length=3):
    ''' A path of catalysis reactions through metabolites (each one's product a substrate of the next), each
    with an enzyme, whose last product acts (as a modifier, interactor, ...; not a substrate) in a reaction of
    another type: [r1, (metabolite, r2), ..., (metabolite, other reaction)]. The longest (up to length
    catalyses), then the fewest participants '''
    def is_metabolite(name):
        return name in adapter.nodes and adapter.nodes[name].type in METABOLITE_TYPES

    def has_enzyme(rid):
        return any(p.edge_type in edge_types.MODIFIERS and p.form != "condition" and not is_metabolite(p.name)
                   for p in adapter.reactions[rid].participants)
    following, ends = {}, {}
    for (down_type, down_effect, role), links in chain_counts(adapter).get(("catalysis", "activation"), {}).items():
        for r1, entity, r2 in links:
            if not is_metabolite(entity) or not has_enzyme(r1):
                continue
            if down_type == "catalysis" and role == "substrate":
                if has_enzyme(r2):
                    following.setdefault(r1, []).append((entity, r2))
            elif down_type != "catalysis" and role != "substrate":
                ends.setdefault(r1, []).append((entity, r2))

    paths = []
    def extend(path, reactions):
        last = path[-1] if len(path) == 1 else path[-1][1]
        paths.extend(path + [end] for end in ends.get(last, []))
        if len(reactions) < length:
            for entity, r in following.get(last, []):
                if r not in reactions:
                    extend(path + [(entity, r)], reactions + [r])
    for r1 in sorted(following.keys() | ends.keys()):
        extend([r1], [r1])

    def score(path):
        reactions = [path[0]] + [r for _, r in path[1:]]
        return (-len(reactions), sum(len(adapter.reactions[r].participants) for r in reactions), reactions)
    return min(paths, key=score) if paths else None


def draw_metabolite_pathway(adapter, out):
    ''' remove_and_rewire of all metabolites of a metabolic pathway and the reaction its product acts in '''
    from skm_tools.pss import remove_and_rewire, simplify_pss
    path = metabolite_pathway(adapter)
    if path is None:
        return []
    reactions = [path[0]] + [r for _, r in path[1:]]
    with tempfile.TemporaryDirectory() as tmp:
        g = simplify_pss(load_interaction_network(adapter, reactions, tmp))
    metabolites = [n for n in g if g.nodes[n].get("node_type") in METABOLITE_TYPES]
    after = g.copy()
    remove_and_rewire(after, metabolites)
    half = halves(node_labels(g))
    pos = spine_positions(g, half)
    end = adapter.reactions[reactions[-1]]
    fig, axes = figure(2, f"remove_and_rewire of the metabolites of {' → '.join(reactions[:-1])}, and "
                          f"{reactions[-1]} ({end.reaction_type})", stacked=True)
    for ax, (title, h, highlight) in zip(axes, [("simplify_pss(g)", g, set(metabolites)),
                                                 (f"remove_and_rewire(g, metabolites)", after, ())]):
        setup_panel(ax, f"{title}; {h.number_of_edges()} edges", {n: pos[n] for n in g}, half)
        draw_influences(ax, h, pos, highlight, label_of=reaction_ids_label)
    name = "remove_and_rewire-metabolites"
    svg = save(fig, out / name)
    print(f"remove_and_rewire (metabolites): {' '.join(reactions)} -> {svg}")
    return [{"name": name, "title": "remove_and_rewire: the metabolites of a pathway", "svg": svg,
             "text": f"Removing all metabolites, e.g. to show transcriptomics data on the genes without clutter. "
                     f"Example: the pathway "
                     f"{' → '.join(reactions[:-1])}, whose product {label(adapter, path[-1][0])} acts in "
                     f"{reactions[-1]} ({end.reaction_type}); removed: "
                     f"{', '.join(label(adapter, m) for m in metabolites)}."}]


def draw_duplicated_binding(adapter, binding_reaction, out):
    ''' remove_duplicated_binding_edges on the interaction network of a binding reaction '''
    from skm_tools.pss import remove_duplicated_binding_edges
    with tempfile.TemporaryDirectory() as tmp:
        g = load_interaction_network(adapter, [binding_reaction], tmp)
    after = g.copy()
    removed = remove_duplicated_binding_edges(after)
    half = halves(node_labels(g))
    pos = flow_positions(g, half)
    fig, axes = figure(2, f"remove_duplicated_binding_edges: {reaction_title(adapter.reactions[binding_reaction])}")
    for ax, (title, h) in zip(axes, [("Interaction network", g),
                                     (f"remove_duplicated_binding_edges(g): {len(removed)} removed", after)]):
        setup_panel(ax, f"{title}; {h.number_of_edges()} edges", {n: pos[n] for n in g}, half)
        draw_influences(ax, h, pos)
    svg = save(fig, out / "remove_duplicated_binding_edges")
    print(f"remove_duplicated_binding_edges: {binding_reaction} -> {svg}")
    return {"name": "remove_duplicated_binding_edges", "title": "remove_duplicated_binding_edges", "svg": svg,
            "text": "Binding partners are linked in both directions (dashed). "
                    "remove_duplicated_binding_edges keeps only one of the two edges, for a less cluttered "
                    "drawing; the edges to the complex stay. Use it for drawing only, not before path or "
                    "neighbourhood analysis, where the removed direction would be missing."}


# ---------------------------------------------------------------------------
# The document
# ---------------------------------------------------------------------------

def file_name(reaction_type, effect):
    return f"{reaction_type}-{effect}".replace("/", "_").replace(" ", "_")


def git_version():
    try:
        return subprocess.run(["git", "-C", str(DOCS), "describe", "--always", "--dirty"], capture_output=True,
                              text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def svg_img(path, alt):
    data = base64.b64encode(Path(path).read_bytes()).decode()
    return f'<img src="data:image/svg+xml;base64,{data}" alt="{html.escape(alt)}">'


def rules_html(rules, caption="rules"):
    if not rules:
        return f'<p class="rules">No {caption} for this type and effect.</p>'
    rows = "".join(f"<tr><td>{html.escape(s)}</td><td>→</td><td>{html.escape(t)}</td><td>{html.escape(i)}</td>"
                   f"<td>{html.escape(d)}</td></tr>" for s, t, i, d in rules)
    return f'<table class="rules"><caption>{caption}</caption>{rows}</table>'


def figure_html(record, notes):
    parts = [f'<figure id="{record["name"]}"><figcaption><h3>{html.escape(record["title"])}</h3>']
    if record.get("text"):
        # (text with markup, e.g. <code>, is marked "html")
        parts.append(f"<p>{record['text'] if record.get('html') else html.escape(record['text'])}</p>")
    if "rules" in record:
        parts.append(rules_html(record["rules"], "rules" + (" (upstream reaction)" if "rules2" in record else "")))
    if "rules2" in record:
        parts.append(rules_html(record["rules2"], "rules (downstream reaction)"))
    if record["name"] in notes:
        parts.append(f'<div class="note">{notes[record["name"]]}</div>')
    parts.append("</figcaption>")
    if record.get("svg"):
        parts.append(svg_img(record["svg"], record["title"]))
    parts.append("</figure>")
    return "\n".join(parts)


CSS = """
:root { --fg: #1d1d1f; --muted: #6e6e73; --line: #e2e2e6; --note: #fff8e6; --note-line: #f0c36d; }
body { font-family: system-ui, -apple-system, "Segoe UI", sans-serif; color: var(--fg); background: #fff;
       margin: 0; line-height: 1.45; font-size: 15px; }
main { max-width: 17cm; margin: 2em auto; padding: 0 16px; }
nav.toc { font-size: 0.82em; padding: 1em 16px; border-bottom: 1px solid var(--line); }
nav.toc ol { padding-left: 1.2em; margin: 0; } nav.toc ul { padding-left: 1em; margin: 0.2em 0 0.6em; }
nav.toc li { margin: 0.15em 0; } nav.toc a { color: #2a4fb8; text-decoration: none; }
nav.toc a:hover { text-decoration: underline; } nav.toc > ol > li > a { font-weight: 600; }
.toc-title { font-weight: 600; color: var(--muted); margin: 0 0 0.5em; text-transform: uppercase;
             letter-spacing: 0.05em; font-size: 0.9em; }
@media (min-width: 1100px) {
  nav.toc { position: fixed; top: 0; left: 0; bottom: 0; width: 270px; overflow-y: auto; overflow-x: hidden; border-bottom: none;
            border-right: 1px solid var(--line); background: #fafafb; box-sizing: border-box; }
  main { margin-left: max(300px, calc(50vw - 8.5cm + 135px)); }
}
h2, h3, figure { scroll-margin-top: 1em; }
h1 { font-size: 1.7em; margin-bottom: 0.2em; } h2 { margin-top: 2.2em; border-bottom: 1px solid var(--line); }
h3 { font-size: 1em; margin: 0 0 0.3em; }
h3.group { font-size: 1.3em; margin: 2.2em 0 0.8em; padding: 0.35em 0.7em; background: #eef2fb;
           border-left: 5px solid #4169e1; } .meta, caption { color: var(--muted); font-size: 0.85em; }
figure { margin: 1.6em 0; break-inside: avoid; page-break-inside: avoid; }
figure img { width: 100%; height: auto; border: 1px solid var(--line); }
figcaption p { margin: 0.2em 0; }
table { border-collapse: collapse; font-size: 0.85em; margin: 0.4em 0; }
td, th { padding: 2px 8px; border-bottom: 1px solid var(--line); text-align: left; vertical-align: top; }
table.rules td { font-family: ui-monospace, monospace; border: none; padding: 0 6px; }
caption { text-align: left; }
.note { background: var(--note); border-left: 3px solid var(--note-line); padding: 0.4em 0.8em; margin: 0.4em 0; }
nav ol { padding-left: 1.2em; } details { margin: 0.6em 0; }
@media print { nav.toc { display: none; } main { margin: 0; max-width: none; } body { font-size: 10pt; }
               h2 { break-before: page; } }
"""


def group_id(group):
    return "".join(c if c.isalnum() else "_" for c in group)


def side_panel(types, chain_records, skm_records):
    ''' The links to the sections and figures, for the side panel '''
    def link(target, text):
        return f'<li><a href="#{target}">{html.escape(text)}</a></li>'
    groups = list(dict.fromkeys(r["group"] for r in chain_records))
    return ('<nav class="toc"><p class="toc-title">Contents</p><ol>'
            '<li><a href="#projections">Projections of each reaction type</a><ul>'
            + "".join(link(r["name"], r["title"]) for r in types) + '</ul></li>'
            '<li><a href="#chains">How reactions chain</a><ul>'
            + "".join(link(f"from-{group_id(g)}", f"from {g}") for g in groups) + '</ul></li>'
            '<li><a href="#skm-tools">skm-tools: simplifying and rewiring</a><ul>'
            + "".join(link(r["name"], r["title"]) for r in skm_records) + '</ul></li></ol></nav>')


def write_html(path, adapter, legend, types, chain_records, chain_table, skm_records, top):
    notes = yaml.safe_load(NOTES.read_text()) if NOTES.exists() else {}
    notes = notes or {}
    sections = []

    sections.append('<h2 id="projections">1. Projections of each reaction type</h2>')
    sections.append(f"<p>One example reaction for each reaction type and effect, drawn in the three networks: "
                    f"the reaction graph, the interaction network and the gene network (species: "
                    f"{html.escape(adapter.species)}). Nodes are in the same place in each: inputs above, outputs "
                    f"below, activators right, inhibitors left. In the gene network, the genes of a gene cluster "
                    f"are drawn around the cluster's name.</p>")
    sections.append("<p>The table above each figure lists the rules that make the interaction network edges "
                    "for that reaction type: from which participants to which, with which influence. "
                    "<em>substrates</em> are the inputs, <em>products</em> the outputs and <em>modifiers</em> "
                    "the activators and inhibitors (<em>activating_modifiers</em> and "
                    "<em>inhibiting_modifiers</em> when the two are treated differently). <em>directed</em> "
                    "edges go one way; <em>mutual</em> edges link each pair of participants in both directions."
                    "</p>")
    sections += [figure_html(r, notes) for r in types]

    sections.append('<h2 id="chains">2. How reactions chain in the interaction network</h2>')
    sections.append(f"<p>Reactions connect when an entity made by one reaction (upstream) is used by another "
                    f"(downstream), as an input or a modifier. This section shows, for each reaction type, its "
                    f"{top} most common downstream connections, each with one example pair of reactions drawn in "
                    f"the interaction network. The shared entity is outlined in orange, and edges are labelled "
                    f"with their reaction id and type.</p>"
                    f"<p>The table lists all connections: the upstream reaction type, the downstream reaction "
                    f"type, the role of the shared entity in the downstream reaction, the number of reaction "
                    f"pairs connected this way, and an example (upstream reaction, shared entity, downstream "
                    f"reaction). The counts are also in <code>chains/chains.tsv</code>.</p>")
    summary = chain_table.groupby(["type", "effect"], sort=True)
    rows = []
    for (t, e), group in summary:
        for i, row in enumerate(group.itertuples()):
            rows.append(f"<tr><td>{html.escape(t) + ', ' + html.escape(e) if i == 0 else ''}</td>"
                        f"<td>{html.escape(row.downstream_type)}, {html.escape(row.downstream_effect)}</td>"
                        f"<td>{html.escape(row.role)}</td><td>{row.links}</td><td>{html.escape(row.example)}</td></tr>")
    sections.append("<details><summary>All links (" + str(len(chain_table)) + "; click to expand)</summary><table><tr><th>upstream"
                    "</th><th>downstream</th><th>role downstream</th><th>pairs</th><th>example</th></tr>" + "".join(rows)
                    + "</table></details>")
    group = None
    for r in chain_records:
        if r["group"] != group:
            group = r["group"]
            sections.append(f'<h3 class="group" id="from-{group_id(group)}">From {html.escape(group)}</h3>')
        sections.append(figure_html(r, notes))

    sections.append('<h2 id="skm-tools">3. skm-tools: simplifying and rewiring</h2>')
    sections.append("<p>The skm-tools functions that simplify the interaction network for analysis or "
                    "drawing, each shown on an example from the data, before and after. Edge labels are reaction "
                    "ids.</p>")
    sections += [figure_html(r, notes) for r in skm_records]

    intro = notes.get("introduction", "")
    doc = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>PSS Network Projections</title><style>{CSS}</style></head><body>
{side_panel(types, chain_records, skm_records)}
<main>
<h1>PSS network projections</h1>
<p class="meta">Generated by <code>docs/projection_figures.py</code> (skm-pss-adapters {html.escape(git_version())})
on {datetime.date.today().isoformat()}, from {len(adapter.reaction_ids)} reactions (access: {html.escape(adapter.access)},
species: {html.escape(adapter.species)}).</p>
{f'<div class="note">{intro}</div>' if intro else ''}
<figure>{svg_img(legend, "Legend")}</figure>
{chr(10).join(sections)}
</main>
</body></html>
"""
    path.write_text(doc)
    return path


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--type", help="only the figures of this reaction type, e.g. catalysis (no document)")
    parser.add_argument("--reaction", nargs="+", help="only figures of these reactions (no document)")
    parser.add_argument("--species", default="ath", help="species of the gene network (default ath)")
    parser.add_argument("--access", default="public")
    parser.add_argument("--top", type=int, default=3, help="chain figures of the top N links per type (default 3)")
    parser.add_argument("--out", type=Path, default=OUT)
    args = parser.parse_args()

    adapter = PSSAdapter()
    adapter.collect_reactions(access=args.access, species=args.species)
    adapter.access = getattr(adapter, "access", None) or args.access

    if args.reaction or args.type:
        draw_types(adapter, args.out / "types", args.type, args.reaction)
        if args.type:
            draw_chains(adapter, args.out / "chains", args.top, args.type)
        return

    for folder in ("types", "chains", "skm-tools"):      # no figures left from earlier examples
        for old in (args.out / folder).glob("*.[ps][nv]g"):
            old.unlink()
    legend = draw_legend(args.out / "legend")
    types = draw_types(adapter, args.out / "types")
    # more examples, each after the example of its type: a binding with three or more interactors, a
    # translocation without modifiers
    for after, draw in (("binding_oligomerisation-inhibition", draw_binding_partners),
                        ("translocation-activation", draw_unmodified_translocation)):
        at = next(i for i, r in enumerate(types) if r["name"] == after) + 1
        types[at:at] = draw(adapter, args.out / "types")
    types += draw_autoregulation(adapter, args.out / "types")
    chain_records, chain_table, drawn = draw_chains(adapter, args.out / "chains", args.top)
    skm_records = draw_simplify(adapter, args.out / "skm-tools")
    skm_records += draw_rewire(adapter, drawn, args.out / "skm-tools")
    skm_records += draw_metabolite_pathway(adapter, args.out / "skm-tools")
    # on the binding with three or more interactors (mutual edges between each pair)
    binding = binding_partners_reaction(adapter)
    if binding is not None:
        skm_records.append(draw_duplicated_binding(adapter, binding, args.out / "skm-tools"))
    path = write_html(args.out / "index.html", adapter, legend, types, chain_records, chain_table, skm_records,
                      args.top)
    print(f"-> {path}")


if __name__ == "__main__":
    main()
