'''Draws docs/pss-export-overview.svg (how pss-export works). Run: python docs/make_overview_svg.py'''
from pathlib import Path
from html import escape

W, H = 1400, 750
FILL = {"use": "#fdeee8", "in": "#e8f0fa", "cfg": "#fdf6e3", "data": "#eef6ee", "write": "#ffffff", "out": "#f3eefa"}
parts = []


def box(x, y, w, h, title, lines=(), kind="write", dashed=False, note=False):
    style = 'stroke-dasharray="5,3"' if dashed else ""
    if note:
        parts.append(f'<path d="M{x},{y} h{w-12} l12,12 v{h-12} h-{w} z" fill="{FILL[kind]}" class="box"/>'
                     f'<path d="M{x+w-12},{y} v12 h12" fill="none" class="box"/>')
    else:
        parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="{FILL[kind]}" class="box" {style}/>')
    ty = y + 18
    parts.append(f'<text x="{x + w/2}" y="{ty}" class="title">{title}</text>')
    for i, line in enumerate(lines):
        parts.append(f'<text x="{x + w/2}" y="{ty + 15 + i*13}" class="small">{escape(line)}</text>')
    return (x, y, w, h)


def cylinder(x, y, w, h, title, lines):
    parts.append(f'<path d="M{x},{y+8} v{h-16} a{w/2},8 0 0 0 {w},0 v-{h-16}" fill="{FILL["in"]}" class="box"/>'
                 f'<ellipse cx="{x+w/2}" cy="{y+8}" rx="{w/2}" ry="8" fill="{FILL["in"]}" class="box"/>')
    parts.append(f'<text x="{x + w/2}" y="{y + 32}" class="title">{title}</text>')
    for i, line in enumerate(lines):
        parts.append(f'<text x="{x + w/2}" y="{y + 47 + i*13}" class="small">{escape(line)}</text>')
    return (x, y, w, h)


def group(x, y, w, h, label, dashed=True):
    style = 'stroke-dasharray="6,4"' if dashed else ""
    parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" class="group" {style}/>'
                 f'<text x="{x + 12}" y="{y + 18}" class="group-label">{escape(label)}</text>')


def arrow(points, label=None, dotted=False, dashed=False, both=False, label_at=None):
    d = "M" + " L".join(f"{px},{py}" for px, py in points)
    style = ' stroke-dasharray="2,4"' if dotted else (' stroke-dasharray="5,3"' if dashed else "")
    start = ' marker-start="url(#arrow-start)"' if both else ""
    parts.append(f'<path d="{d}" class="edge"{style} marker-end="url(#arrow)"{start}/>')
    if label:
        lx, ly = label_at or ((points[0][0] + points[-1][0]) / 2, (points[0][1] + points[-1][1]) / 2 - 5)
        parts.append(f'<text x="{lx}" y="{ly}" class="edge-label">{escape(label)}</text>')


def right(b, dy=0): x, y, w, h = b; return (x + w, y + h/2 + dy)
def left(b, dy=0): x, y, w, h = b; return (x, y + h/2 + dy)
def bottom(b, dx=0): x, y, w, h = b; return (x + w/2 + dx, y + h)
def top(b, dx=0): x, y, w, h = b; return (x + w/2 + dx, y)


# columns
group(15, 35, 250, 420, "Used through")
group(285, 35, 230, 330, "Inputs")
group(535, 35, 290, 695, "1. Collect: collect_reactions(…)", dashed=False)
group(845, 35, 300, 695, "2. Write: create_…()", dashed=False)
group(1165, 35, 220, 695, "Outputs")
group(852, 58, 286, 370, "Networks: everything")
group(852, 440, 286, 280, "Models: nodes_to_ignore left out")

cli = box(30, 60, 220, 70, "CLI: pss-export", ["to-sbml, to-tabularqual, to-sbgn,", "to-reaction-graph[-json], to-faidare,", "to-interaction-network, to-gene-network"], "use")
reg = box(30, 150, 220, 70, "Format registry", ["pss_export.formats: titles, descriptions,", "files, access, model fixes", "(e.g. for the SKM web app)"], "use")
api = box(30, 245, 220, 90, "Python: PSSAdapter", ["collect_reactions()", "with_model_fixes()", "create_sbml(), create_…()", "export(key, files…, model_fixes)"], "use")
scr = box(30, 365, 220, 70, "scripts/export_pss.py", ["the reference exports in pss/:", "every format but FAIDARE,", "public, Arabidopsis"], "use")

pss = cylinder(300, 60, 200, 80, "PSS (Neo4j)", ["reactions, participants,", "entities, annotations"])
cfg = box(300, 160, 200, 90, "pss_export_config.yaml", ["species, nodes_to_ignore", "SBO: forms, roles,", "subtypes, mechanisms", "interaction_rules"], "cfg", note=True)
defs = box(300, 270, 200, 80, "pss_reaction_definitions.py", ["reaction and edge types,", "participant roles,", "subtypes (keys of the SBO tables)"], "cfg", note=True)

col = box(550, 60, 260, 100, "PSSCollector (Cypher)", ["filters: access (public / restricted),", "reactions or pathways, species (all", "functional clusters have genes);", "nodes_to_ignore: models only"], "in")
dat = box(550, 190, 260, 125, "Collected data", ["Reaction: type, effect, mechanism, links", "· participants (as collected)", "· substrates / products / modifiers", "· roles, reaction SBO, mechanism SBO", "node annotations: links, homologues", "metadata: version, access, species"], "data")
fix = box(550, 345, 260, 85, "Model fixes (opt-in)", ["with_model_fixes(): a fixed copy for the", "\"with model fixes\" variant of the models:", "active forms, transport reactions,", "activation from the inactive form"], "write", dashed=True)
ann = box(550, 560, 260, 70, "annotation_registry.yaml", ["databases → qualifier,", "identifiers.org prefix"], "cfg", note=True)

rg   = box(862, 82, 266, 50, "Reaction graph", ["entities and reactions, one edge per participant"])
rgj  = box(862, 150, 266, 45, "Reaction graph (JSON)", ["the same graph: the PSS Explorer's data"])
inet = box(862, 215, 266, 60, "Interaction network", ["entity → entity influences (positive /", "negative / unknown), by interaction_rules"])
gnet = box(862, 295, 266, 50, "Gene network", ["functional clusters → genes in the species"])
fd   = box(862, 365, 266, 50, "FAIDARE", ["an entry per gene and species (all species)"])
sbml = box(862, 465, 266, 60, "SBML", ["species per compartment, SBO terms,", "RDF annotations, notes, metadata"])
bm   = box(862, 545, 266, 70, "Boolean model", ["a rule per species: any activation", "and no inhibition (TabularQual, BoolNet)"])
sbgn = box(862, 635, 266, 70, "SBGN-ML", ["Process Description, Graphviz layout,", "compartments, annotations;", "a variant for Newt"])

o_rg  = box(1180, 82, 190, 50, "edges.tsv + nodes.tsv", ["extended SIF"], "out", note=True)
o_rgj = box(1180, 148, 190, 50, "graph.json", ["nodes and edges"], "out", note=True)
o_in  = box(1180, 220, 190, 50, "edges.tsv + nodes.tsv", ["extended SIF"], "out", note=True)
o_gn  = box(1180, 295, 190, 50, "edges.tsv + nodes.tsv", ["extended SIF, per species"], "out", note=True)
o_fd  = box(1180, 365, 190, 50, "entries.json", ["FAIDARE data discovery"], "out", note=True)
o_sb  = box(1180, 470, 190, 50, "model.xml", ["SBML L3V2"], "out", note=True)
o_tq  = box(1180, 530, 190, 45, "model.xlsx", ["TabularQual"], "out", note=True)
o_bn  = box(1180, 585, 190, 45, "model.bnet + nodes.tsv", ["BoolNet"], "out", note=True)
o_sg  = box(1180, 645, 190, 50, "map.sbgn", ["SBGN-ML (or for Newt)"], "out", note=True)

arrow([bottom(cli), top(reg)], "calls", label_at=(170, 145))
arrow([bottom(reg), top(api)], "export(key)", label_at=(180, 239))
arrow([top(scr), bottom(api)], "export(key)", label_at=(180, 357))
arrow([right(api, -35), (522, right(api, -35)[1]), (522, 130), (550, 130)], "collect_reactions()", label_at=(400, 266))
arrow([right(pss), (550, right(pss)[1])])
arrow([right(cfg), (538, right(cfg)[1]), (538, 150), (550, 150)], dotted=True)
arrow([bottom(col), top(dat)])
arrow([right(defs), (530, right(defs)[1]), (530, 285), (550, 285)], dotted=True)
arrow([bottom(dat), top(fix)], dashed=True)
for w in (rg, inet, fd, sbml, bm, sbgn):
    arrow([right(dat), (832, right(dat)[1]), (832, left(w, -12)[1]), left(w, -12)])
arrow([right(fix), (840, right(fix)[1]), (840, 446), (852, 446)], dashed=True)
arrow([bottom(rg), top(rgj)])
arrow([bottom(inet), top(gnet)], "same edges, expanded", label_at=(1060, 288))
for w in (sbml, bm, sbgn):
    arrow([right(ann), (848, right(ann)[1]), (848, left(w, 20)[1]), left(w, 20)], dotted=True)
for w, o in ((rg, o_rg), (rgj, o_rgj), (inet, o_in), (gnet, o_gn), (fd, o_fd), (sbml, o_sb), (sbgn, o_sg)):
    arrow([right(w), (o[0], right(w)[1])])
arrow([right(bm, -12), (1155, right(bm, -12)[1]), (1155, right(o_tq)[1] - 0), (1180, right(o_tq)[1])])
arrow([right(bm, 12), (1155, right(bm, 12)[1]), (1155, right(o_bn)[1]), (1180, right(o_bn)[1])])

svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}" font-family="Helvetica, Arial, sans-serif">
<title>How pss-export works</title>
<defs>
  <marker id="arrow" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M0,0 L10,5 L0,10 z" fill="#5a6b7b"/></marker>
  <marker id="arrow-start" viewBox="0 0 10 10" refX="1" refY="5" markerWidth="7" markerHeight="7" orient="auto"><path d="M10,0 L0,5 L10,10 z" fill="#5a6b7b"/></marker>
  <style>
    .box {{ stroke: #5a6b7b; stroke-width: 1.2; }}
    .group {{ fill: #f8fafb; stroke: #9aa7b3; stroke-width: 1; }}
    .group-label {{ font-size: 12px; fill: #44525e; }}
    .title {{ font-size: 12.5px; font-weight: bold; text-anchor: middle; fill: #1f2a33; }}
    .small {{ font-size: 10px; text-anchor: middle; fill: #2f3b45; }}
    .edge {{ fill: none; stroke: #5a6b7b; stroke-width: 1.2; }}
    .edge-label {{ font-size: 10px; fill: #44525e; text-anchor: middle; }}
  </style>
</defs>
<rect width="{W}" height="{H}" fill="#ffffff"/>
{chr(10).join(parts)}
</svg>
'''
Path(__file__).with_name("pss-export-overview.svg").write_text(svg)
