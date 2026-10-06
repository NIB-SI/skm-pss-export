'''Draws docs/pss-export-overview.svg (how pss-export works). Run: python docs/make_overview_svg.py'''
from pathlib import Path
from html import escape

W, H = 1400, 600
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
group(15, 35, 250, 330, "Used through")
group(285, 35, 230, 330, "Inputs")
group(535, 35, 290, 540, "1. Collect: collect_reactions(…)", dashed=False)
group(845, 35, 300, 540, "2. Write: create_…()", dashed=False)
group(1165, 35, 220, 540, "Outputs")

cli = box(30, 60, 220, 70, "CLI: pss-export", ["to-sbml, to-tabularqual, to-reaction-graph,", "to-interaction-network, to-gene-network,", "formats"], "use")
reg = box(30, 150, 220, 70, "Format registry", ["pss_export.formats: titles, descriptions,", "files, access", "(e.g. for the SKM web app)"], "use")
api = box(30, 250, 220, 90, "Python: PSSAdapter", ["collect_reactions()", "create_sbml(), create_…()", "export(key, files…)"], "use")

pss = cylinder(300, 60, 200, 80, "PSS (Neo4j)", ["reactions, participants,", "entities, annotations"])
defs = box(300, 270, 200, 80, "pss_reaction_definitions.py", ["reaction and edge types,", "participant roles,", "subtypes (keys of the SBO tables)"], "cfg", note=True)
cfg = box(300, 160, 200, 90, "pss_export_config.yaml", ["species, nodes_to_ignore", "SBO: forms, roles,", "subtypes, mechanisms", "interaction_rules"], "cfg", note=True)

col = box(550, 60, 260, 100, "PSSCollector (Cypher)", ["filters: access (public / restricted),", "reactions or pathways, species (all", "functional clusters have genes);", "nodes_to_ignore: models only"], "in")
dat = box(550, 190, 260, 125, "Collected data", ["Reaction: type, effect, mechanism, links", "· participants (as collected)", "· substrates / products / modifiers", "· roles, reaction SBO, mechanism SBO", "node annotations: links, homologues", "metadata: version, access, species"], "data")
fix = box(550, 340, 260, 55, "model_fixes (optional)", ["active forms, transport reactions"], "write", dashed=True)
ann = box(550, 450, 260, 70, "annotation_registry.yaml", ["databases → qualifier,", "identifiers.org prefix"], "cfg", note=True)

rg   = box(860, 60, 270, 55, "Reaction graph", ["entities and reactions, one edge per participant"])
inet = box(860, 140, 270, 70, "Interaction network", ["entity → entity influences", "(positive / negative / unknown),", "by interaction_rules"])
gnet = box(860, 250, 270, 55, "Gene network", ["functional clusters → genes in the species"])
sbml = box(860, 375, 270, 70, "SBML", ["species per compartment, SBO terms,", "RDF annotations, notes, metadata"])
tq   = box(860, 470, 270, 70, "TabularQual", ["species, Boolean transition rules,", "model information"])

o_rg = box(1180, 60, 190, 55, "edges.tsv + nodes.tsv", ["extended SIF"], "out", note=True)
o_in = box(1180, 147, 190, 55, "edges.tsv + nodes.tsv", ["extended SIF"], "out", note=True)
o_gn = box(1180, 250, 190, 55, "edges.tsv + nodes.tsv", ["extended SIF, per species"], "out", note=True)
o_sb = box(1180, 382, 190, 55, "model.xml", ["SBML L3V2"], "out", note=True)
o_tq = box(1180, 477, 190, 55, "model.xlsx", ["TabularQual"], "out", note=True)

arrow([bottom(cli), top(reg)], "calls", label_at=(170, 145))
arrow([bottom(reg), top(api)], "export(key)", label_at=(180, 243))
arrow([right(api, -35), (522, right(api, -35)[1]), (522, 130), (550, 130)], "collect_reactions()", label_at=(400, 256))
arrow([right(pss), (550, right(pss)[1])])
arrow([right(cfg), (538, right(cfg)[1]), (538, 150), (550, 150)], dotted=True)
arrow([bottom(col), top(dat)])
arrow([right(defs), (530, right(defs)[1]), (530, 285), (550, 285)], dotted=True)
arrow([bottom(dat), top(fix)], dashed=True, both=True)
for w in (rg, inet, sbml, tq):
    arrow([right(dat), (832, right(dat)[1]), (832, left(w, -12)[1]), left(w, -12)])
arrow([bottom(inet), top(gnet)], "same edges, expanded", label_at=(1060, 232))
arrow([right(ann, -12), (848, right(ann, -12)[1]), (848, left(sbml, 16)[1]), left(sbml, 16)], dotted=True)
arrow([right(ann, 12), (848, right(ann, 12)[1]), (848, left(tq, 16)[1]), left(tq, 16)], dotted=True)
for w, o in ((rg, o_rg), (inet, o_in), (gnet, o_gn), (sbml, o_sb), (tq, o_tq)):
    arrow([right(w), (o[0], right(w)[1])])

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
