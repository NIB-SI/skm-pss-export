"""
SBGN-ML export: glyphs and arcs by form and role, compartments, annotations, layout.
Needs Graphviz (sfdp).
"""

import shutil
from types import SimpleNamespace
from xml.etree import ElementTree as ET

import pytest

import pss_export.entity_classes as ec
from pss_export.sbgn import sbgn_api

pytestmark = pytest.mark.skipif(not shutil.which("sfdp"), reason="needs Graphviz (sfdp)")

NS = {"s": sbgn_api.SBGN_NS, "rdf": sbgn_api.RDF_NS, "bqbiol": sbgn_api.BQ_NS["bqbiol"], "sio": sbgn_api.SIO_NS}
CLUSTER = ["Node", "FunctionalCluster", "Plant", "PlantCoding"]


def species(name, form, location, edge_type):
    s = ec.Species(name, form, location)
    s.edge_type = edge_type
    return s


def reaction(rid, rtype, effect, participants):
    r = ec.Reaction(rid, rtype, {"reaction_effect": effect, "external_links": ["pubmed:1"]})
    for p in participants:
        p = species(*p)
        r.participants.append(p)
        if p.edge_type in ("SUBSTRATE", "TRANSLOCATE_FROM"):
            if not (rtype.startswith("transcriptional") and p.form == "gene"):    # as collected, include_genes=False
                r.add_substrate(p)
        elif p.edge_type in ("PRODUCT", "TRANSLOCATE_TO"):
            r.add_product(p)
        else:
            r.add_modifier(p)
    return r


NODES = {
    "MYC2": ec.Node("MYC2", labels=CLUSTER, functional_cluster_id="fc00001", external_links=["uniprot:Q39204"]),
    "JAZ1": ec.Node("JAZ1", labels=CLUSTER),
    "JA-Ile": ec.Node("JA-Ile", labels=["Node", "Metabolite"], external_links=["CHEBI:51141"]),
    "COI1|JAZ1": ec.Node("COI1|JAZ1", labels=["Node", "Complex"], components=["COI1", "JAZ1"]),
    "COI1": ec.Node("COI1", labels=CLUSTER),
}


def adapter(*reactions):
    return SimpleNamespace(
        nodes=NODES, reactions={r.reaction_id: r for r in reactions},
        model_reaction_ids=[r.reaction_id for r in reactions], additional_reactions=[],
        include_genes=False, nodes_to_ignore=[], species=None, reaction_pathways={"rx1": ["JA"]},
        model_id="pss", model_name="PSS", model_description=None, model_version=None, access="public",
        species_description="all species", export_datetime="2026-10-06", model_fixes_applied=None)


def write(tmp_path, *reactions, newt=False):
    path = tmp_path / "map.sbgn"
    assert sbgn_api.create_sbgn(adapter(*reactions), path, newt=newt) == len(reactions)
    return ET.parse(path).getroot().find("s:map", NS)


TRANSCRIPTION = reaction("rx1", "transcriptional/translational activation", "activation",
                         [("JAZ1", "gene", "nucleus", "SUBSTRATE"), ("JAZ1", "protein", "cytoplasm", "PRODUCT"),
                          ("MYC2", "protein_active", "nucleus", "ACTIVATES")])
BINDING = reaction("rx2", "binding/oligomerisation", "activation",
                   [("JA-Ile", "metabolite", "cytoplasm", "SUBSTRATE"), ("JAZ1", "protein", "cytoplasm", "SUBSTRATE"),
                    ("COI1|JAZ1", "complex_active", "cytoplasm", "PRODUCT")])
DEGRADATION = reaction("rx3", "degradation/secretion", "inhibition",
                       [("JAZ1", "protein", "cytoplasm", "SUBSTRATE"), ("COI1|JAZ1", "complex_active", "cytoplasm", "ACTIVATES")])


def test_glyphs_and_arcs(tmp_path):
    m = write(tmp_path, TRANSCRIPTION, BINDING, DEGRADATION)
    glyphs = {g.get("id"): g for g in m.findall("s:glyph", NS)}
    arcs = {(a.get("source"), a.get("class"), a.get("target")) for a in m.findall("s:arc", NS)}
    # transcription: the gene template as necessary stimulation, made from a source and sink
    assert ("s_JAZ1_nuc_g", "necessary stimulation", "rx1") in arcs
    assert ("rx1.source", "consumption", "rx1.in") in arcs and glyphs["rx1.source"].get("class") == "source and sink"
    assert ("rx1.out", "production", "s_JAZ1_cyt_p") in arcs
    assert ("s_MYC2_nuc_pa", "stimulation", "rx1") in arcs
    # binding: association; degradation: a sink
    assert glyphs["rx2"].get("class") == "association"
    assert ("rx3.out", "production", "rx3.sink") in arcs
    assert ("s_COI1JAZ1_cyt_ca", "stimulation", "rx3") in arcs
    # form: class and state; a complex with its components
    assert glyphs["s_MYC2_nuc_pa"].get("class") == "macromolecule"
    assert glyphs["s_MYC2_nuc_pa"].find("s:glyph/s:state", NS).get("value") == "active"
    complex_ = glyphs["s_COI1JAZ1_cyt_ca"]
    assert [c.find("s:label", NS).get("text") for c in complex_.findall("s:glyph", NS) if c.get("class") == "macromolecule"] == ["COI1", "JAZ1"]
    assert glyphs["s_JA-Ile_cyt_m".replace("-", "")].get("class") == "simple chemical"


def test_compartments(tmp_path):
    m = write(tmp_path, TRANSCRIPTION)
    glyphs = {g.get("id"): g for g in m.findall("s:glyph", NS)}
    assert glyphs["s_MYC2_nuc_pa"].get("compartmentRef") == "c_nuc"
    assert glyphs["c_nuc"].get("compartmentRef") == "c_cyt"
    assert glyphs["rx1"].get("compartmentRef") is None          # processes are in no compartment

    def box(g):
        return [float(g.find("s:bbox", NS).get(k)) for k in "xywh"]

    x, y, w, h = box(glyphs["s_MYC2_nuc_pa"])
    cx, cy, cw, ch = box(glyphs["c_nuc"])
    assert cx <= x and cy <= y and x + w <= cx + cw and y + h <= cy + ch


def test_annotations(tmp_path):
    m = write(tmp_path, TRANSCRIPTION)
    glyphs = {g.get("id"): g for g in m.findall("s:glyph", NS)}
    links = {li.get(f"{{{sbgn_api.RDF_NS}}}resource") for li in glyphs["s_MYC2_nuc_pa"].iterfind(".//bqbiol:is//rdf:li", NS)}
    assert links == {"http://identifiers.org/uniprot:Q39204", "http://identifiers.org/skm:fc00001"}
    reaction_links = {li.get(f"{{{sbgn_api.RDF_NS}}}resource") for li in glyphs["rx1"].iterfind(".//rdf:li", NS)}
    assert reaction_links == {"http://identifiers.org/skm:rx1", "http://identifiers.org/pubmed:1"}
    assert not list(m.iterfind(".//sio:SIO_000223", NS))           # Newt's custom properties: Newt variant only
    assert m.find("s:extension", NS) is None


def test_newt_variant(tmp_path):
    m = write(tmp_path, TRANSCRIPTION, newt=True)
    properties = {li.get(f"{{{sbgn_api.SIO_NS}}}SIO_000116"): li.get(f"{{{sbgn_api.RDF_NS}}}value")
                  for li in m.iterfind(".//s:glyph[@id='s_MYC2_nuc_pa']//sio:SIO_000223//rdf:li", NS)}
    assert properties["name"] == "MYC2" and properties["location"] == "nucleus"
    render = m.find(f"s:extension/{{{sbgn_api.RENDER_NS}}}renderInformation", NS)
    styled = {s.get("idList") for s in render.iter(f"{{{sbgn_api.RENDER_NS}}}style")}
    assert {g.get("id") for g in m.iter(f"{{{sbgn_api.SBGN_NS}}}glyph")} <= styled     # Newt: a style for every glyph


def test_common_compartment():
    assert sbgn_api.common_compartment(["nucleus", "cytoplasm"]) == "cytoplasm"
    assert sbgn_api.common_compartment(["nucleolus", "nucleus"]) == "nucleus"
    assert sbgn_api.common_compartment(["nucleus", "nucleus"]) == "nucleus"
    assert sbgn_api.common_compartment(["cytoplasm", "extracellular"]) is None
