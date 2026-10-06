'''
SBGN-ML (Process Description) of the PSS model reactions, laid out with Graphviz.

The same reactions and participants as SBML (the model reactions, nodes_to_ignore left out, model fixes in the
"with model fixes" variant), with the same species ids. Glyphs:

- entity pool nodes by form (protein: macromolecule with an active/inactive state variable, metabolite: simple
  chemical, gene/mRNA/miRNA/ncRNA: nucleic acid feature, complex: complex with its components, process and foreign
  entities: unspecified entity, abiotic: perturbing agent), in their cellular compartment;
- a process node per reaction (association for binding, dissociation for dissociation), with an input and an
  output port;
- arcs by the participant's role: substrate/interactor: consumption, product: production, template (the gene of a
  transcription/translation): necessary stimulation, with a source and sink as the process's input, catalyst and
  transporter: catalysis, stimulator: stimulation, inhibitor: inhibition, modifier: stimulation or inhibition by its
  edge (ACTIVATES/INHIBITS); degradation/secretion produces a source and sink. Conditions are left out.

Annotations: MIRIAM (bqbiol/bqmodel RDF in each glyph's extension, identifiers.org URLs, from the CURIEs of the
database, as in SBML), and notes (descriptions, evidence). The Newt variant (newt=True) adds Newt's custom properties
(sio:SIO_000223) and its render information (colours by form and compartment).

Layout: Graphviz sfdp (force-directed), compartment by compartment, bottom-up: each compartment's content is laid
out on its own, then the compartment is one node in its parent's layout. (dot, layered, with the compartments as
clusters, lines the map up into a strip about 20 times wider than high; fdp, which supports clusters, places
connected glyphs barely closer than at random.) Coordinates are SBGN-ML's: the top left corner of each bounding box, y from top to bottom.
'''

import json
import shutil
import subprocess
from collections import defaultdict
from xml.etree import ElementTree as ET

from ..annotations.annotation_manager import annotation_manager
from ..entity_classes import IDTracker
from ..pss.config import pss_export_config
from ..pss.pss_reaction_definitions import participant_roles, reaction_classes
from ..utils import MODEL_FIXES_NOTE

SBGN_NS = 'http://sbgn.org/libsbgn/0.3'
PD_VERSION = 'http://identifiers.org/combine.specifications/sbgn.pd.level-1.version-1.3'
RDF_NS = 'http://www.w3.org/1999/02/22-rdf-syntax-ns#'
BQ_NS = {'bqbiol': 'http://biomodels.net/biology-qualifiers/', 'bqmodel': 'http://biomodels.net/model-qualifiers/'}
SIO_NS = 'http://semanticscience.org/resource/'
XHTML_NS = 'http://www.w3.org/1999/xhtml'
RENDER_NS = 'http://www.sbml.org/sbml/level3/version1/render/version1'

for prefix, uri in [('', SBGN_NS), ('rdf', RDF_NS), ('sio', SIO_NS), *BQ_NS.items()]:
    ET.register_namespace(prefix, uri)

# entity pool node class and state variable by participant form
FORM_GLYPH = {
    'metabolite': ('simple chemical', None),
    'protein': ('macromolecule', 'inactive'),
    'protein_active': ('macromolecule', 'active'),
    'gene': ('nucleic acid feature', None),
    'mrna': ('nucleic acid feature', None),
    'mirna': ('nucleic acid feature', None),
    'ncrna': ('nucleic acid feature', None),
    'complex': ('complex', 'inactive'),
    'complex_active': ('complex', 'active'),
    # unspecified entities have no state variables in SBGN PD: the state is in the label
    'process': ('unspecified entity', None),
    'process_active': ('unspecified entity', None),
    'foreign_entity': ('unspecified entity', None),
    'abiotic': ('perturbing agent', None),
}
ACTIVE_LABEL = {'process_active': 'active'}

# the class of a complex's component, by its PSS class
COMPONENT_GLYPH = {
    'Metabolite': 'simple chemical', 'MetaboliteFamily': 'simple chemical',
    'PlantCoding': 'macromolecule', 'ForeignCoding': 'macromolecule', 'PlantAbstract': 'macromolecule',
    'PlantNonCoding': 'nucleic acid feature', 'ForeignNonCoding': 'nucleic acid feature',
    'Complex': 'complex',
}

PROCESS_GLYPH = {'binding/oligomerisation': 'association', 'dissociation': 'dissociation'}

# arc class by the participant's role (modifier: by its edge)
ROLE_ARC = {
    participant_roles.SUBSTRATE: 'consumption',
    participant_roles.INTERACTOR: 'consumption',
    participant_roles.PRODUCT: 'production',
    participant_roles.TEMPLATE: 'necessary stimulation',
    participant_roles.CATALYST: 'catalysis',
    participant_roles.TRANSPORTER: 'catalysis',
    participant_roles.STIMULATOR: 'stimulation',
    participant_roles.INHIBITOR: 'inhibition',
}

# compartments in compartments (the others are at the top level: cytoplasm, extracellular, apoplast)
PARENT_COMPARTMENT = {
    'chloroplast': 'cytoplasm', 'endoplasmic reticulum': 'cytoplasm', 'golgi apparatus': 'cytoplasm',
    'mitochondrion': 'cytoplasm', 'nucleus': 'cytoplasm', 'peroxisome': 'cytoplasm', 'vacuole': 'cytoplasm',
    'nucleolus': 'nucleus',
}

# glyph sizes (SBGN-ML units)
SIZE = {
    'macromolecule': (120, 40), 'nucleic acid feature': (120, 40), 'simple chemical': (90, 40),
    'unspecified entity': (120, 40), 'perturbing agent': (120, 50), 'source and sink': (24, 24),
    'process': (20, 20), 'association': (20, 20), 'dissociation': (20, 20),
}
COMPONENT_SIZE = (110, 30)
COMPLEX_PADDING = 10
PORT_DISTANCE = 15
STATE_SIZE = (40, 16)

# Newt colours
FORM_COLOUR = {
    'metabolite': '#abb0d0', 'protein': '#7fffd4', 'protein_active': '#7fffd4',
    'gene': '#006b47', 'mrna': '#006b47', 'mirna': '#77c0c9', 'ncrna': '#77c0c9',
    'process': '#a239fa', 'process_active': '#a239fa', 'complex': '#ffefd5', 'complex_active': '#ffefd5',
    'foreign_entity': '#f2a1bd', 'abiotic': '#f2a1bd',
}
COMPARTMENT_COLOUR = {
    'cytoplasm': '#c8c831', 'chloroplast': '#228b22', 'endoplasmic reticulum': '#f16f36',
    'golgi apparatus': '#e69500', 'nucleus': '#8b2257', 'peroxisome': '#2200e6', 'vacuole': '#31c8c8',
    'mitochondrion': '#adff2f', 'apoplast': '#1b1b07', 'nucleolus': '#c8317d', 'extracellular': '#f5f5dc',
}
COMPONENT_COLOUR = {'simple chemical': '#abb0d0', 'macromolecule': '#7fffd4', 'nucleic acid feature': '#77c0c9',
                    'complex': '#ffefd5'}


NEWT_MAP_PROPERTIES = [
    ('compoundPadding', '0'), ('extraCompartmentPadding', '14'), ('extraComplexPadding', '10'),
    ('arrowScale', '1.25'), ('showComplexName', 'true'), ('dynamicLabelSize', 'regular'),
    ('inferNestingOnLoad', 'false'), ('fitLabelsToNodes', 'false'), ('fitLabelsToInfoboxes', 'false'),
    ('recalculateLayoutOnComplexityManagement', 'false'), ('rearrangeOnComplexityManagement', 'true'),
    ('animateOnDrawingChanges', 'true'), ('adjustNodeLabelFontSizeAutomatically', 'false'),
    ('enablePorts', 'true'), ('enableSIFTopologyGrouping', 'false'), ('allowCompoundNodeResize', 'true'),
    ('mapColorScheme', 'black_white'), ('mapColorSchemeStyle', 'solid'), ('mapName', 'PSS'),
    ('mapDescription', ''), ('highlightColor', '#ff0000'), ('extraHighlightThickness', '10'),
]


def compartment_path(name):
    ''' The compartment and the compartments it is in, outermost first '''
    path = [name]
    while path[0] in PARENT_COMPARTMENT:
        path.insert(0, PARENT_COMPARTMENT[path[0]])
    return path


def common_compartment(names):
    ''' The innermost compartment that has all of them (e.g. nucleus and cytoplasm: cytoplasm), or None '''
    paths = [compartment_path(name) for name in names]
    common = None
    for level in zip(*paths) if paths else []:
        if len(set(level)) > 1:
            break
        common = level[0]
    return common


def q(tag, ns=SBGN_NS):
    return f'{{{ns}}}{tag}'


class SBGN(IDTracker):
    ''' The SBGN map of the model reactions of a PSSAdapter (create_sbgn) '''

    def __init__(self, pss_adapter, newt=False):
        IDTracker.__init__(self)
        self.pss_adapter = pss_adapter
        self.newt = newt
        self.glyphs = {}            # id: {'class', 'label', 'compartment', ...}
        self.arcs = []              # {'id', 'class', 'source', 'target'}
        self.compartments = {}      # name: glyph id
        self.reaction_pathways = {}  # reaction id: pathway (Newt properties)
        self.arc_counter = 0

    # --- building ---------------------------------------------------------------------------------------------

    def node(self, name):
        return self.pss_adapter.nodes.get(name)

    def compartment(self, name):
        ''' The glyph id of a cellular compartment (and of its parents) '''
        if name not in self.compartments:
            parent = PARENT_COMPARTMENT.get(name)
            parent_id = self.compartment(parent) if parent else None
            id_ = f"c_{pss_export_config.compartment_to_short[name]}"
            self.compartments[name] = id_
            self.glyphs[id_] = {'class': 'compartment', 'label': name, 'compartment': parent_id, 'name': name}
        return self.compartments[name]

    def species_glyph(self, species):
        ''' The glyph of a participant (one per name, form and compartment, with the SBML species id) '''
        id_, known = self.get_species_id(species)
        if known:
            return id_
        self.set_species_id(species, id_)
        class_, state = FORM_GLYPH.get(species.form, ('unspecified entity', None))
        node = self.node(species.name)
        label = node.display_label if node else IDTracker.get_display_label(species.name)
        if species.form in ACTIVE_LABEL:
            label = f'{label} ({ACTIVE_LABEL[species.form]})'
        components = []
        if class_ == 'complex' and node:
            for component in node.components:
                component_node = self.node(component)
                components.append({
                    'label': component_node.display_label if component_node else IDTracker.get_display_label(component),
                    'class': COMPONENT_GLYPH.get(component_node.type if component_node else None, 'unspecified entity'),
                })
        self.glyphs[id_] = {
            'class': class_, 'label': label, 'state': state, 'components': components,
            'compartment': self.compartment(species.compartment), 'form': species.form,
            'location': species.compartment, 'name': species.name, 'node': node,
        }
        return id_

    def source_and_sink(self, id_, compartment):
        self.glyphs[id_] = {'class': 'source and sink', 'label': None, 'compartment': self.compartment(compartment)}
        return id_

    def arc(self, class_, source, target):
        self.arc_counter += 1
        self.arcs.append({'id': f'a{self.arc_counter}', 'class': class_, 'source': source, 'target': target})

    def add_reaction(self, reaction):
        ''' A process node with its participants (the model participants) '''
        reaction_id = reaction.reaction_id
        self.set_reaction_id(reaction, reaction_id)
        substrates = list(reaction.substrates)
        templates = []
        if reaction.reaction_type in reaction_classes.TRANSCRIPTIONAL_TRANSLATIONAL:
            templates = [s for s in substrates if s.form == 'gene']
            substrates = [s for s in substrates if s.form != 'gene']
            if not self.pss_adapter.include_genes:
                # the models leave the gene templates out: SBGN shows them, as necessary stimulation
                ignore = set(self.pss_adapter.nodes_to_ignore or [])
                templates = [p for p in reaction.participants if p.edge_type == 'SUBSTRATE' and p.form == 'gene'
                             and p.name not in ignore]
        participants = substrates + templates + list(reaction.products) + list(reaction.modifiers)
        compartments = [p.compartment for p in participants if p.form != 'condition']
        # laid out in the innermost compartment with all its participants, else (e.g. apoplast and cytoplasm) in
        # its first participant's
        home = common_compartment(compartments) or (compartments[0] if compartments else None)

        self.glyphs[reaction_id] = {
            'class': PROCESS_GLYPH.get(reaction.reaction_type, 'process'), 'label': None,
            # not in a compartment (SBGN; for Newt, which would infer one, it is), but laid out in one
            'layout_compartment': self.compartment(home) if home else None,
            'reaction': reaction, 'ports': True,
        }
        ins, outs = f'{reaction_id}.in', f'{reaction_id}.out'
        where = home or (substrates[0].compartment if substrates else 'cytoplasm')

        for species in substrates:
            self.arc('consumption', self.species_glyph(species), ins)
        if not substrates:
            # made from nothing (transcription/translation, or no substrate curated)
            self.arc('consumption', self.source_and_sink(f'{reaction_id}.source', where), ins)
        for species in templates:
            self.arc('necessary stimulation', self.species_glyph(species), reaction_id)
        for species in reaction.products:
            self.arc('production', outs, self.species_glyph(species))
        if not reaction.products:
            # degradation/secretion (or no product curated)
            self.arc('production', outs, self.source_and_sink(f'{reaction_id}.sink', where))
        for species in reaction.modifiers:
            if species.form == 'condition':
                continue
            role = participant_roles.of(reaction, species)
            arc_class = ROLE_ARC.get(role) or ('inhibition' if species.edge_type == 'INHIBITS' else 'stimulation')
            self.arc(arc_class, self.species_glyph(species), reaction_id)

    # --- layout -----------------------------------------------------------------------------------------------

    def size(self, glyph):
        if glyph['class'] == 'complex':
            n = max(len(glyph['components']), 1)
            return (COMPONENT_SIZE[0] + 2 * COMPLEX_PADDING,
                    n * COMPONENT_SIZE[1] + (n + 1) * COMPLEX_PADDING + STATE_SIZE[1])
        return SIZE[glyph['class']]

    def layout(self):
        ''' The bounding box (top left, w, h) of every glyph and compartment, bottom-up: the content of each
        compartment is laid out on its own (sfdp, force-directed), and the compartment then takes part in its
        parent's layout as one node of that size; arcs between compartments pull them together. '''
        if not shutil.which('sfdp'):
            raise RuntimeError("SBGN layout needs Graphviz (the 'sfdp' program)")
        children = defaultdict(list)
        for id_, glyph in self.glyphs.items():
            parent = glyph['layout_compartment'] if 'layout_compartment' in glyph else glyph.get('compartment')
            children[parent].append(id_)
        # which layout node stands for a glyph in a compartment's layout: itself, or the child compartment it is in
        container = {id_: (g['layout_compartment'] if 'layout_compartment' in g else g.get('compartment'))
                     for id_, g in self.glyphs.items()}
        ends = [(arc['source'].removesuffix('.in').removesuffix('.out'),
                 arc['target'].removesuffix('.in').removesuffix('.out')) for arc in self.arcs]

        def representative(id_, compartment):
            while id_ is not None and container.get(id_) != compartment:
                id_ = container.get(id_)
            return id_

        relative = {}       # id: (x, y) of its centre in its compartment's layout
        sizes = {}

        def lay_out(compartment):
            members = children[compartment]
            for member in members:
                if self.glyphs[member]['class'] == 'compartment':
                    lay_out(member)
                else:
                    sizes[member] = self.size(self.glyphs[member])
            edges = set()
            for source, target in ends:
                s, t = representative(source, compartment), representative(target, compartment)
                if s and t and s != t:
                    edges.add((s, t))
            lines = ['graph g {', '  graph [overlap=prism, sep="+12", splines=false];',
                     '  node [shape=box, fixedsize=true, label=""];']
            lines += [f'  "{m}" [width={sizes[m][0] / 72:.3f}, height={sizes[m][1] / 72:.3f}];' for m in members]
            lines += [f'  "{s}" -- "{t}";' for s, t in sorted(edges)]
            result = subprocess.run(['sfdp', '-Tjson'], input='\n'.join(lines + ['}']),
                                    capture_output=True, text=True, check=True)
            positions = {o['name']: tuple(float(v) for v in o['pos'].split(','))
                         for o in json.loads(result.stdout).get('objects', []) if 'pos' in o}
            # the box around the members, with a margin; y from top to bottom
            margin = 30 if compartment else 0
            left = min(positions[m][0] - sizes[m][0] / 2 for m in members) - margin
            right = max(positions[m][0] + sizes[m][0] / 2 for m in members) + margin
            bottom = min(positions[m][1] - sizes[m][1] / 2 for m in members) - margin
            top = max(positions[m][1] + sizes[m][1] / 2 for m in members) + margin
            for m in members:
                relative[m] = (positions[m][0] - left, top - positions[m][1])
            if compartment:
                sizes[compartment] = (right - left, top - bottom + (30 if margin else 0))
            return right - left, top - bottom

        width, height = lay_out(None)
        self.map_bbox = (0, 0, width, height)

        def place(compartment, x0, y0):
            for m in children[compartment]:
                w, h = sizes[m]
                x, y = x0 + relative[m][0] - w / 2, y0 + relative[m][1] - h / 2
                self.glyphs[m]['bbox'] = (x, y, w, h)
                if self.glyphs[m]['class'] == 'compartment':
                    place(m, x, y + 30)     # room for the compartment's label

        place(None, 0, 0)

    # --- writing ----------------------------------------------------------------------------------------------

    @staticmethod
    def bbox(parent, x, y, w, h):
        ET.SubElement(parent, q('bbox'), x=f'{x:.1f}', y=f'{y:.1f}', w=f'{w:.1f}', h=f'{h:.1f}')

    def annotation(self, element, id_, links, properties):
        ''' MIRIAM annotations (and, for Newt, its custom properties) in the element's extension '''
        annotations, skipped = annotation_manager.process_node(links)
        for link in skipped:
            print(f"SBGN: warning, skipping or could not parse external link of {id_}: {link}")
        properties = {k: v for k, v in properties.items() if v} if self.newt else {}
        if not annotations and not properties:
            return
        extension = ET.SubElement(element, q('extension'))
        # the namespaces again on rdf:RDF: Newt (libsbgn.js) parses it on its own
        rdf = ET.SubElement(ET.SubElement(extension, q('annotation')), q('RDF', RDF_NS), {
            'xmlns:rdf': RDF_NS, 'xmlns:bqbiol': BQ_NS['bqbiol'], 'xmlns:bqmodel': BQ_NS['bqmodel'],
            'xmlns:sio': SIO_NS})
        description = ET.SubElement(rdf, q('Description', RDF_NS), {q('about', RDF_NS): f'#{id_}'})
        by_qualifier = defaultdict(list)
        for annotation in annotations:
            by_qualifier[(annotation.namespace, annotation.qualifier)].append(annotation.url)
        for (namespace, qualifier), urls in by_qualifier.items():
            bag = ET.SubElement(ET.SubElement(description, q(qualifier, BQ_NS[namespace])), q('Bag', RDF_NS))
            for url in urls:
                ET.SubElement(bag, q('li', RDF_NS), {q('resource', RDF_NS): url})
        for key, value in properties.items():
            bag = ET.SubElement(ET.SubElement(description, q('SIO_000223', SIO_NS)), q('Bag', RDF_NS))
            ET.SubElement(bag, q('li', RDF_NS), {q('SIO_000116', SIO_NS): key, q('value', RDF_NS): str(value)})

    @staticmethod
    def notes(element, items):
        items = [(k, v) for k, v in items if v]
        if not items:
            return
        body = ET.SubElement(ET.SubElement(element, q('notes')), q('body', XHTML_NS))
        for key, value in items:
            ET.SubElement(body, q('p', XHTML_NS)).text = f'{key}: {value}'

    def glyph_element(self, map_element, id_, glyph):
        attributes = {'class': glyph['class'], 'id': id_}
        # entities, and compartments in compartments (nucleus in cytoplasm); processes only for Newt, which
        # otherwise puts them in the first compartment whose box has them
        if glyph.get('compartment') and 'layout_compartment' not in glyph:
            attributes['compartmentRef'] = glyph['compartment']
        elif self.newt and glyph.get('layout_compartment'):
            attributes['compartmentRef'] = glyph['layout_compartment']
        if glyph['class'] in ('process', 'association', 'dissociation'):
            attributes['orientation'] = 'vertical'
        element = ET.SubElement(map_element, q('glyph'), attributes)

        node, reaction = glyph.get('node'), glyph.get('reaction')
        if node:
            self.notes(element, [('description', node.description)])
            self.annotation(element, id_, node.links(self.pss_adapter.species),
                            {'name': node.name, 'location': glyph['location'], 'pathway': node.pathway,
                             'functional_cluster_id': node.functional_cluster_id})
        elif reaction:
            self.notes(element, [('reaction', reaction.reaction_id), ('type', reaction.reaction_type),
                                 ('mechanism', reaction.reaction_mechanism),
                                 ('evidence', reaction.evidence_sentence), ('export notes', reaction.export_notes)])
            self.annotation(element, id_, [f'skm:{reaction.reaction_id}'] + list(reaction.external_links or []),
                            {'name': reaction.reaction_id, 'reaction_type': reaction.reaction_type,
                             'pathway': self.reaction_pathways.get(reaction.reaction_id)})

        if glyph['label']:
            ET.SubElement(element, q('label'), text=glyph['label'])
        x, y, w, h = glyph['bbox']
        self.bbox(element, x, y, w, h)

        # state variable (top centre), then complex components (stacked)
        if glyph.get('state'):
            state = ET.SubElement(element, q('glyph'), {'class': 'state variable', 'id': f'{id_}.state'})
            ET.SubElement(state, q('state'), value=glyph['state'])
            self.bbox(state, x + (w - STATE_SIZE[0]) / 2, y - STATE_SIZE[1] / 2, *STATE_SIZE)
        for i, component in enumerate(glyph.get('components', [])):
            child = ET.SubElement(element, q('glyph'), {'class': component['class'], 'id': f'{id_}.c{i}'})
            ET.SubElement(child, q('label'), text=component['label'])
            self.bbox(child, x + COMPLEX_PADDING,
                      y + STATE_SIZE[1] / 2 + COMPLEX_PADDING + i * (COMPONENT_SIZE[1] + COMPLEX_PADDING),
                      *COMPONENT_SIZE)

        if glyph.get('ports'):
            cx = x + w / 2
            ET.SubElement(element, q('port'), id=f'{id_}.in', x=f'{cx:.1f}', y=f'{y - PORT_DISTANCE:.1f}')
            ET.SubElement(element, q('port'), id=f'{id_}.out', x=f'{cx:.1f}', y=f'{y + h + PORT_DISTANCE:.1f}')

    def point(self, ref, towards):
        ''' Where an arc meets a glyph (or port): the port, or the glyph's border towards the other end '''
        if ref.endswith('.in') or ref.endswith('.out'):
            x, y, w, h = self.glyphs[ref.rsplit('.', 1)[0]]['bbox']
            return x + w / 2, (y - PORT_DISTANCE if ref.endswith('.in') else y + h + PORT_DISTANCE)
        x, y, w, h = self.glyphs[ref]['bbox']
        cx, cy = x + w / 2, y + h / 2
        tx, ty = towards
        dx, dy = tx - cx, ty - cy
        if dx == 0 and dy == 0:
            return cx, cy
        scale = min(w / 2 / abs(dx) if dx else float('inf'), h / 2 / abs(dy) if dy else float('inf'))
        return cx + dx * scale, cy + dy * scale

    def centre(self, ref):
        if ref.endswith('.in') or ref.endswith('.out'):
            return self.point(ref, (0, 0))
        x, y, w, h = self.glyphs[ref]['bbox']
        return x + w / 2, y + h / 2

    def arc_element(self, map_element, arc):
        element = ET.SubElement(map_element, q('arc'), {'class': arc['class'], 'id': arc['id'],
                                                        'source': arc['source'], 'target': arc['target']})
        sx, sy = self.point(arc['source'], self.centre(arc['target']))
        tx, ty = self.point(arc['target'], self.centre(arc['source']))
        ET.SubElement(element, q('start'), x=f'{sx:.1f}', y=f'{sy:.1f}')
        ET.SubElement(element, q('end'), x=f'{tx:.1f}', y=f'{ty:.1f}')

    def render_extension(self, map_element):
        ''' Newt's render information: colours by form, compartment borders by compartment (Newt needs a style
        for every glyph and arc) '''
        colours = {'color_dark': '#555555', 'colour_white': '#ffffff', 'colour_compartment': '#ebebeb',
                   **{f'form_{k}': v for k, v in FORM_COLOUR.items()},
                   **{f'compartment_{k.replace(" ", "_")}': v for k, v in COMPARTMENT_COLOUR.items()},
                   **{f'component_{k.replace(" ", "_")}': v for k, v in COMPONENT_COLOUR.items()}}
        styles = []
        for id_, glyph in self.glyphs.items():
            style = {'stroke': 'color_dark', 'stroke-width': '1.25', 'fill': 'colour_white'}
            if glyph['class'] == 'compartment':
                style.update({'stroke': f'compartment_{glyph["name"].replace(" ", "_")}', 'stroke-width': '3',
                              'fill': 'colour_compartment'})
            elif glyph.get('form') in FORM_COLOUR:
                style['fill'] = f'form_{glyph["form"]}'
            styles.append((id_, style))
            if glyph.get('state'):
                styles.append((f'{id_}.state', {'stroke': 'color_dark', 'stroke-width': '1', 'fill': 'colour_white'}))
            for i, component in enumerate(glyph.get('components', [])):
                styles.append((f'{id_}.c{i}', {'stroke': 'color_dark', 'stroke-width': '1',
                                               'fill': f'component_{component["class"].replace(" ", "_")}'}))
        styles += [(arc['id'], {'stroke': 'color_dark', 'stroke-width': '1.25'}) for arc in self.arcs]

        extension = ET.SubElement(map_element, q('extension'))
        render = ET.SubElement(extension, q('renderInformation', RENDER_NS),
                               {'id': 'renderInformation', 'program-name': 'pss-export',
                                'background-color': '#00000000'})
        colour_list = ET.SubElement(render, q('listOfColorDefinitions', RENDER_NS))
        for id_, value in colours.items():
            ET.SubElement(colour_list, q('colorDefinition', RENDER_NS), id=id_, value=value)
        style_list = ET.SubElement(render, q('listOfStyles', RENDER_NS))
        for id_, style in styles:
            element = ET.SubElement(style_list, q('style', RENDER_NS), id=f'style.{id_}', idList=id_)
            ET.SubElement(element, q('g', RENDER_NS), style)
        # Newt's map properties (all of them: Newt fails on missing ones, e.g. the paddings)
        properties = ET.SubElement(extension, 'mapProperties')
        for key, value in NEWT_MAP_PROPERTIES:
            ET.SubElement(properties, key).text = value

    def write(self, filename):
        adapter = self.pss_adapter
        root = ET.Element(q('sbgn'))
        map_element = ET.SubElement(root, q('map'), id=adapter.model_id or 'pss', language='process description',
                                    version=PD_VERSION)
        self.notes(map_element, [
            ('name', adapter.model_name), ('description', adapter.model_description),
            ('version', adapter.model_version), ('access', adapter.access), ('species', adapter.species_description),
            ('source', 'https://skm.nib.si'), ('export date', adapter.export_datetime),
            ('model fixes', MODEL_FIXES_NOTE.format(n=adapter.model_fixes_applied)
                            if adapter.model_fixes_applied is not None else None)])
        if self.newt:
            self.render_extension(map_element)
        self.bbox(map_element, *self.map_bbox)

        # compartments first (outer before inner), then the other glyphs
        depth = {id_: 0 if name not in PARENT_COMPARTMENT else 1 + (name == 'nucleolus')
                 for name, id_ in self.compartments.items()}
        order = sorted(self.glyphs, key=lambda i: (self.glyphs[i]['class'] != 'compartment', depth.get(i, 0), i))
        for id_ in order:
            self.glyph_element(map_element, id_, self.glyphs[id_])
        for arc in self.arcs:
            self.arc_element(map_element, arc)

        tree = ET.ElementTree(root)
        ET.indent(tree, space='  ')
        tree.write(filename, encoding='utf-8', xml_declaration=True)


def create_sbgn(pss_adapter, filename, newt=False):
    ''' Write the SBGN-ML map of the model reactions (see the module). Returns the number of process nodes. '''
    sbgn = SBGN(pss_adapter, newt=newt)
    sbgn.reaction_pathways = {r: ', '.join(sorted(p for p in pathways if p))
                              for r, pathways in (pss_adapter.reaction_pathways or {}).items()}
    reaction_ids = list(pss_adapter.model_reaction_ids) + list(pss_adapter.additional_reactions)
    for reaction_id in reaction_ids:
        sbgn.add_reaction(pss_adapter.reactions[reaction_id])
    sbgn.layout()
    sbgn.write(filename)
    return len(reaction_ids)
