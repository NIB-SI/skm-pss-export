'''
The export formats pss-export makes: what they are, their files, and how to make them.

The CLI and the SKM web app's downloads page read this. Titles and descriptions are Markdown.

    from pss_export.formats import FORMATS, get_format

    for fmt in FORMATS.values():
        print(fmt.key, fmt.title, [f.extension for f in fmt.files])

    adapter.collect_reactions(access="public", species="stu")      # species: every format
    adapter.export("gene-network", edges_file="edges.tsv", nodes_file="nodes.tsv")

All formats are made for a species (collect_reactions(species=...), one of pss_schema_config.species, default 'ath';
None: no species filter, not for the gene network). The models (model_fixes=True) also come "with model fixes":

    adapter.export("sbml", filename="model.xml")                                    # exactly PSS
    adapter.export("sbml", filename="model-with-model-fixes.xml", model_fixes=True)  # connected

'''

from dataclasses import dataclass
from typing import Tuple

ACCESS_LEVELS = ('public', 'restricted')


@dataclass(frozen=True)
class FormatFile:
    ''' One file of a format '''
    argument: str       # the PSSAdapter method's argument for its path, e.g. "edges_file"
    name: str           # e.g. "edges", "nodes", "model"
    extension: str      # e.g. "tsv", "xml", "xlsx"
    media_type: str     # e.g. "text/tab-separated-values"
    description: str = ''


@dataclass(frozen=True)
class ExportFormat:
    key: str                              # e.g. "interaction-network"
    title: str                            # Markdown
    description: str                      # Markdown
    method: str                           # the PSSAdapter method that writes it
    files: Tuple[FormatFile, ...]
    access: Tuple[str, ...] = ACCESS_LEVELS   # access levels it can be made for
    filterable: bool = True               # can be limited to reactions / pathways
    model_fixes: bool = False             # also has a "with model fixes" variant (export(..., model_fixes=True))


TSV = 'text/tab-separated-values'

EXTENDED_SIF = ('Tab-separated, extended [SIF](https://manual.cytoscape.org/en/stable/Supported_Network_File_Formats.html#sif-format) '
                '(as [Pathway Commons](https://www.pathwaycommons.org/)): a header, the first three columns are '
                'the source, the interaction and the target. Opens in [Cytoscape](https://cytoscape.org) and networkX.')

FORMATS = {f.key: f for f in [

    ExportFormat(
        key='sbml',
        title='[SBML](https://sbml.org)',
        description=('Systems Biology Markup Language (Level 3 Version 2): species per compartment, reactions with '
                     'their participants, SBO terms and annotations. The automated export has no equations or kinetic '
                     'parameters; for work on this, see [Generating dynamic large-scale systems-biology data from a '
                     'knowledge graph](https://github.com/R4d0slav/ThesisRepository).'),
        method='create_sbml',
        files=(FormatFile('filename', 'model', 'xml', 'application/sbml+xml'),),
        model_fixes=True,
    ),

    ExportFormat(
        key='tabularqual',
        title='[TabularQual](https://github.com/sys-bio/TabularQual)',
        description=('A spreadsheet format for logical models (SBML-qual), readable and editable by hand: species, '
                     'Boolean transition rules built from the reactions, and the model information. For curation and '
                     'annotation, for tools such as [BoolDog](https://nib-si.github.io/BoolDog/), or converted to '
                     'SBML-qual with the [TabularQual converter](https://github.com/sys-bio/TabularQual).'),
        method='create_tabularqual',
        files=(FormatFile('filename', 'model', 'xlsx',
                          'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet'),),
        model_fixes=True,
    ),

    ExportFormat(
        key='reaction-graph',
        title='Reaction graph',
        description=('Entities and reactions as nodes, one edge per participant, as in the database and the PSS '
                     'Explorer: nothing is left out (conditions and gene templates included). The second column, `role`, '
                     'is the participant\'s role (substrate, product, catalyst, …). ' + EXTENDED_SIF),
        method='create_reaction_graph',
        files=(FormatFile('edges_file', 'edges', 'tsv', TSV, 'one edge per reaction participant'),
               FormatFile('nodes_file', 'nodes', 'tsv', TSV, 'entities and reactions')),
    ),

    ExportFormat(
        key='interaction-network',
        title='Interaction network',
        description=('Entities as nodes, and their influences on each other through the reactions as edges (an SBGN '
                     'Activity Flow view): `positive-influence`, `negative-influence` or `unknown-influence`, with '
                     'the reaction and its SBO terms. The edges each reaction type gives are documented with '
                     '[pss-export](https://github.com/NIB-SI/skm-pss-export). ' + EXTENDED_SIF),
        method='create_interaction_network',
        files=(FormatFile('edges_file', 'edges', 'tsv', TSV, 'entity -> entity influences'),
               FormatFile('nodes_file', 'nodes', 'tsv', TSV, 'entities')),
    ),

    ExportFormat(
        key='gene-network',
        title='Gene interaction network',
        description=('The interaction network with functional clusters expanded into their genes (from the '
                     'clusters\' homologues) in the species the export is made for; for '
                     '[DiNAR](https://bio.tools/dinar) and CKN. The same edge columns as the '
                     'interaction network. Genes are annotated with their functional cluster(s) only. ' + EXTENDED_SIF),
        method='create_gene_network',
        files=(FormatFile('edges_file', 'edges', 'tsv', TSV, 'gene -> gene (and other node) influences'),
               FormatFile('nodes_file', 'nodes', 'tsv', TSV, 'genes and other nodes')),
    ),
]}


def get_format(key):
    if key not in FORMATS:
        raise ValueError(f"Unknown export format '{key}', one of: {', '.join(FORMATS)}")
    return FORMATS[key]
