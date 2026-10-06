'''
Network exports (reaction graph, interaction network, gene network): the rules per reaction
type, without a database. See the interaction_rules in pss_export_config.yaml.
'''
import csv
from types import SimpleNamespace

import pytest

import pss_export.entity_classes as ec
from pss_export.networks import network_api as net


def species(name, form, location, edge_type):
    s = ec.Species(name, form, location)
    s.edge_type = edge_type
    return s


def reaction(rid, rtype, effect, participants, mechanism=None):
    r = ec.Reaction(rid, rtype, {"reaction_effect": effect, "reaction_mechanism": mechanism})
    r.participants = [species(*p) for p in participants]
    return r


CLUSTER = ["Node", "FunctionalCluster", "Plant", "PlantCoding"]


def nodes(*specs):
    return {name: ec.Node(name, labels=labels, homologues=homologues) for name, labels, homologues in specs}


NODES = nodes(("ADK", CLUSTER, {"ath": ["AT1", "AT2"], "stu": []}),
              ("PDT", CLUSTER, {"ath": ["AT3", "AT9"]}),
              ("ADT", CLUSTER, {"ath": ["AT9"]}),
              ("cZ", ["Node", "Metabolite"], {}),
              ("cZ-P", ["Node", "Metabolite"], {}))


def edges(r, nodes=NODES):
    return {(e["source"], e["interaction"], e["target"], e["directed"]) for e in net.interaction_edges(r, nodes)}


def test_catalysis():
    r = reaction("rx1", "catalysis", "activation",
                 [("cZ", "metabolite", "cytoplasm", "SUBSTRATE"), ("cZ-P", "metabolite", "cytoplasm", "PRODUCT"),
                  ("ADK", "protein_active", "cytoplasm", "ACTIVATES")])
    assert edges(r) == {("cZ", "positive-influence", "cZ-P", True),
                        ("ADK", "positive-influence", "cZ-P", True),
                        ("ADK", "negative-influence", "cZ", True)}      # no substrate -> catalyst


def test_binding_mutual_both_ways():
    r = reaction("rx2", "binding/oligomerisation", "inhibition",
                 [("A", "protein", "cytoplasm", "SUBSTRATE"), ("B", "protein", "cytoplasm", "SUBSTRATE"),
                  ("A|B", "complex", "cytoplasm", "PRODUCT")])
    assert edges(r) == {("A", "negative-influence", "B", False), ("B", "negative-influence", "A", False)}


def test_binding_activation():
    r = reaction("rx3", "binding/oligomerisation", "activation",
                 [("A", "protein", "cytoplasm", "SUBSTRATE"), ("B", "protein", "cytoplasm", "SUBSTRATE"),
                  ("A|B", "complex", "cytoplasm", "PRODUCT")])
    assert edges(r) == {("A", "unknown-influence", "B", False), ("B", "unknown-influence", "A", False),
                        ("A", "positive-influence", "A|B", True), ("B", "positive-influence", "A|B", True)}


def test_protein_activation_no_self_loop():
    r = reaction("rx4", "protein activation", "activation",
                 [("SNRK2", "protein", "cytoplasm", "SUBSTRATE"), ("SNRK2", "protein_active", "cytoplasm", "PRODUCT"),
                  ("PYL", "protein_active", "cytoplasm", "ACTIVATES")], mechanism="Phosphorylation")
    [edge] = net.interaction_edges(r, {})
    assert (edge["source"], edge["interaction"], edge["target"]) == ("PYL", "positive-influence", "SNRK2")
    assert edge["reaction_sbo"] == "SBO:0000656"                          # the result (activation)
    assert edge["reaction_mechanism_sbo"] == "SBO:0000216"                # the mechanism


def test_autoregulation_kept():
    # WRKY33 activates its own transcription (rx00160): modifier -> its own entity is kept
    r = reaction("rx11", "transcriptional/translational activation", "activation",
                 [("WRKY33", "gene", "nucleus", "SUBSTRATE"), ("WRKY33", "protein_active", "cytoplasm", "PRODUCT"),
                  ("WRKY33", "protein_active", "nucleus", "ACTIVATES")])
    assert edges(r, {}) == {("WRKY33", "positive-influence", "WRKY33", True)}


def test_deactivation_negative_whatever_the_edge_type():
    for edge_type in ["ACTIVATES", "INHIBITS"]:
        r = reaction("rx5", "protein deactivation", "inhibition",
                     [("ARF", "protein_active", "nucleus", "SUBSTRATE"), ("ARF", "protein", "nucleus", "PRODUCT"),
                      ("AUX/IAA", "protein_active", "nucleus", edge_type)])
        assert edges(r, {}) == {("AUX/IAA", "negative-influence", "ARF", True)}


def test_transcription_activation_inhibiting_modifier_and_conditions():
    r = reaction("rx6", "transcriptional/translational activation", "activation",
                 [("JR1", "gene", "nucleus", "SUBSTRATE"), ("JR1", "protein_active", "cytoplasm", "PRODUCT"),
                  ("MYC2", "protein_active", "nucleus", "ACTIVATES"), ("JAZ", "protein_active", "nucleus", "INHIBITS"),
                  ("NPR1 high", "condition", None, "INHIBITS")], mechanism="transcription")
    got = net.interaction_edges(r, {})
    assert {(e["source"], e["interaction"], e["target"]) for e in got} == {
        ("MYC2", "positive-influence", "JR1"), ("JAZ", "negative-influence", "JR1")}   # no condition node
    assert {e["source_role"] for e in got} == {"stimulator", "inhibitor"}
    assert {(e["reaction_sbo"], e["reaction_mechanism_sbo"]) for e in got} == {("SBO:0000589", "SBO:0000183")}


def test_translocation_without_transporter_gives_no_edges():
    r = reaction("rx7", "translocation", "activation",
                 [("EIN2", "protein_active", "endoplasmic reticulum", "TRANSLOCATE_FROM"),
                  ("EIN2", "protein_active", "nucleus", "TRANSLOCATE_TO")])
    assert edges(r, {}) == set()


def test_net_effect_keeps_the_edge_to_the_product():
    # catalysis whose substrate and product are the same entity: catalyst -> product (+) and -> substrate (-)
    r = reaction("rx8", "catalysis", "activation",
                 [("X", "protein", "cytoplasm", "SUBSTRATE"), ("X", "protein_active", "cytoplasm", "PRODUCT"),
                  ("E", "protein_active", "cytoplasm", "ACTIVATES")])
    [edge] = net.interaction_edges(r, {})
    assert (edge["source"], edge["interaction"], edge["target_role"]) == ("E", "positive-influence", "product")


def test_type_without_rules_warns(capsys):
    r = reaction("rx9", "cleavage/auto-cleavage", "activation", [("A", "protein", None, "SUBSTRATE")])
    assert net.interaction_edges(r, {}) == [] and "no interaction rules" in capsys.readouterr().out


def test_location_and_location_putative():
    r = reaction("rx10", "degradation/secretion", "inhibition",
                 [("VPg", "protein", "putative:cytoplasm", "SUBSTRATE"), ("PUB15", "protein_active", None, "ACTIVATES")])
    [edge] = net.interaction_edges(r, {})
    assert (edge["target_location"], edge["target_location_putative"], edge["source_location"]) == ("cytoplasm", True, None)


def adapter(*reactions, nodes=NODES, species="ath"):
    return SimpleNamespace(reactions={r.reaction_id: r for r in reactions},
                           reaction_ids=[r.reaction_id for r in reactions], nodes=nodes, species=species)


def read(path):
    return list(csv.DictReader(open(path), delimiter="\t"))


def test_gene_network(tmp_path):
    kept = reaction("rx1", "catalysis", "activation",
                    [("cZ", "metabolite", "cytoplasm", "SUBSTRATE"), ("cZ-P", "metabolite", "cytoplasm", "PRODUCT"),
                     ("ADK", "protein_active", "cytoplasm", "ACTIVATES")])
    a = adapter(kept)
    net.create_gene_network(a, tmp_path / "ath.tsv", tmp_path / "ath-nodes.tsv")
    rows = read(tmp_path / "ath.tsv")
    assert {(r["source"], r["target"], r["source_entity"], r["source_type"]) for r in rows if r["source_role"] == "catalyst"} == {
        ("AT1", "cZ-P", "ADK", "gene"), ("AT2", "cZ-P", "ADK", "gene"),
        ("AT1", "cZ", "ADK", "gene"), ("AT2", "cZ", "ADK", "gene")}
    assert list(rows[0]) == net.INTERACTION_EDGE_COLUMNS
    nodes = {r["id"]: r for r in read(tmp_path / "ath-nodes.tsv")}
    assert (nodes["AT1"]["node_type"], nodes["AT1"]["display_label"], nodes["cZ"]["node_type"]) == ("gene", "ADK", "Metabolite")
    assert "entity" not in nodes["AT1"]



def test_gene_in_two_clusters(tmp_path):
    r = reaction("rx2", "binding/oligomerisation", "inhibition",
                 [("PDT", "protein", "cytoplasm", "SUBSTRATE"), ("ADT", "protein", "cytoplasm", "SUBSTRATE"),
                  ("PDT|ADT", "complex", "cytoplasm", "PRODUCT")])
    net.create_gene_network(adapter(r), tmp_path / "e.tsv", tmp_path / "n.tsv")
    pairs = {(e["source"], e["target"]) for e in read(tmp_path / "e.tsv")}
    assert ("AT9", "AT9") not in pairs and ("AT3", "AT9") in pairs
    nodes = {n["id"]: n for n in read(tmp_path / "n.tsv")}
    assert nodes["AT9"]["display_label"] == "ADT|PDT"


def test_node_synonyms_all_pathways_components(tmp_path):
    """ lists joined with '|' (chemical and cluster names contain commas, gene symbols ';'); synonyms: the short name first """
    ns = nodes(("PDT", CLUSTER, {"ath": ["AT3", "AT9"]}), ("ADT", CLUSTER, {"ath": ["AT9"]}),
               ("cZ", ["Node", "Metabolite"], {}), ("PDT|ADT", ["Node", "Complex"], {}))
    ns["PDT|ADT"].components = ["ADT", "AHK2,3,4[fc00012]", "cZ"]
    ns["ADT"].functional_cluster_id = "fc00002"
    ns["PDT"].mapman = ["4.1.1_Amino acid metabolism.biosynthesis.prephenate dehydratase"]
    ns["ADT"].mapman = ["4.1.2_Amino acid metabolism.biosynthesis.arogenate dehydratase"]
    ns["PDT"].short_name, ns["PDT"]._synonyms, ns["PDT"].all_pathways = "PDT", ["PDT", "AROGENATE DH"], ["P1", "P2"]
    ns["ADT"]._synonyms, ns["ADT"].all_pathways = ["ADT1"], ["P2"]
    ns["cZ"]._synonyms = ["10,11-EHT"]
    r = reaction("rx2", "binding/oligomerisation", "activation",
                 [("PDT", "protein", "cytoplasm", "SUBSTRATE"), ("ADT", "protein", "cytoplasm", "SUBSTRATE"),
                  ("cZ", "metabolite", "cytoplasm", "SUBSTRATE"), ("PDT|ADT", "complex", "cytoplasm", "PRODUCT")])
    net.create_interaction_network(adapter(r, nodes=ns), tmp_path / "e.tsv", tmp_path / "n.tsv")
    entities = {n["id"]: n for n in read(tmp_path / "n.tsv")}
    assert (entities["PDT"]["synonyms"], entities["PDT"]["all_pathways"]) == ("PDT|AROGENATE DH", "P1|P2")
    assert (entities["cZ"]["synonyms"], entities["PDT"]["genes"]) == ("10,11-EHT", "AT3|AT9")
    assert (entities["PDT|ADT"]["components"], entities["PDT"]["components"]) == ("ADT|AHK2,3,4[fc00012]|cZ", "")
    # the cluster ids of the components that are clusters in the export's nodes (not cZ; AHK2,3,4 not collected here)
    assert entities["PDT|ADT"]["component_cluster_ids"] == "fc00002"
    assert entities["PDT"]["mapman"] == "4.1.1_Amino acid metabolism.biosynthesis.prephenate dehydratase"
    net.create_gene_network(adapter(r, nodes=ns), None, tmp_path / "g.tsv")
    genes = {n["id"]: n for n in read(tmp_path / "g.tsv")}
    assert (genes["AT9"]["synonyms"], genes["AT9"]["all_pathways"]) == ("ADT1|PDT|AROGENATE DH", "P1|P2")
    assert genes["AT9"]["mapman"].split("|") == ["4.1.1_Amino acid metabolism.biosynthesis.prephenate dehydratase",
                                               "4.1.2_Amino acid metabolism.biosynthesis.arogenate dehydratase"]


def test_gene_autoregulation_kept(tmp_path):
    wrky = nodes(("WRKY33", ["FunctionalCluster", "PlantCoding"], {"ath": ["AT2G38470"]}))
    r = reaction("rx11", "transcriptional/translational activation", "activation",
                 [("WRKY33", "gene", "nucleus", "SUBSTRATE"), ("WRKY33", "protein_active", "cytoplasm", "PRODUCT"),
                  ("WRKY33", "protein_active", "nucleus", "ACTIVATES")])
    net.create_gene_network(adapter(r, nodes=wrky), tmp_path / "e.tsv", None)
    assert {(e["source"], e["target"]) for e in read(tmp_path / "e.tsv")} == {("AT2G38470", "AT2G38470")}

    # several genes: each regulates itself only
    wrky["WRKY33"].homologues["ath"] = ["G1", "G2", "G3"]
    net.create_gene_network(adapter(r, nodes=wrky), tmp_path / "e3.tsv", None)
    assert {(e["source"], e["target"]) for e in read(tmp_path / "e3.tsv")} == {("G1", "G1"), ("G2", "G2"), ("G3", "G3")}


def test_gene_network_needs_a_species():
    with pytest.raises(ValueError, match="needs a species"):
        net.create_gene_network(adapter(species=None), None, None)


def test_reaction_has_genes_in():
    r = reaction("rx1", "catalysis", "activation",
                 [("cZ", "metabolite", "cytoplasm", "SUBSTRATE"), ("ADK", "protein_active", "cytoplasm", "ACTIVATES")])
    assert r.has_genes_in(NODES, "ath") and not r.has_genes_in(NODES, "stu")
    assert reaction("rx2", "catalysis", "activation", [("cZ", "metabolite", None, "SUBSTRATE")]).has_genes_in(
        NODES, "stu")                                                        # no functional clusters


def test_reaction_graph(tmp_path):
    r = reaction("rx6", "transcriptional/translational activation", "activation",
                 [("JR1", "gene", "nucleus", "SUBSTRATE"), ("JR1", "protein_active", "cytoplasm", "PRODUCT"),
                  ("MYC2", "protein_active", "nucleus", "ACTIVATES"), ("NPR1 high", "condition", None, "INHIBITS")])
    assert net.create_reaction_graph(adapter(r), tmp_path / "e.tsv", tmp_path / "n.tsv") == 4
    rows = {(e["source"], e["role"], e["target"], e["edge_type"]) for e in read(tmp_path / "e.tsv")}
    assert rows == {("JR1", "template", "rx6", "SUBSTRATE"), ("rx6", "product", "JR1", "PRODUCT"),
                    ("MYC2", "stimulator", "rx6", "ACTIVATES"), ("NPR1 high", "inhibitor", "rx6", "INHIBITS")}
    nodes = {n["id"]: n for n in read(tmp_path / "n.tsv")}
    assert nodes["rx6"]["node_type"] == "reaction" and nodes["rx6"]["reaction_sbo"] == "SBO:0000589"


def test_node_type_most_specific():
    assert ec.Node("x", labels=["Node", "Foreign", "ForeignCoding", "ForeignNonCoding"]).type == "ForeignCoding"
    assert ec.Node("x", labels=["Node", "Metabolite", "MetaboliteFamily"]).type == "MetaboliteFamily"
    assert ec.Node("x", labels=CLUSTER).type == "PlantCoding"


def test_node_label_and_links():
    n = ec.Node("PYL[AT5G46790]", functional_cluster_id="fc00001", external_links=["tair.name:AT5G46790"],
                homologues={"ath": ["AT5G46790"]})
    assert n.display_label == "PYL"
    assert ec.Node("PYL[fc00001]", display_label="PYL").display_label == "PYL"
    assert n.links("ath") == ["tair.name:AT5G46790", "skm:fc00001"]          # not twice
    assert n.links("stu") == ["tair.name:AT5G46790", "skm:fc00001"]        # curated TAIR link kept
    links = n.links("ath"); links.append("x:y")
    assert n.links("ath") == ["tair.name:AT5G46790", "skm:fc00001"]   # a new list


def test_gene_rows_one_entry_per_cluster(tmp_path):
    """ a gene in several clusters: entity, display_label, short_name, pathway, functional_cluster_id have an
    entry per cluster, in the same order, also when empty; labels with commas stay one entry """
    ns = nodes(("TCP8,14,15[fc00332]", CLUSTER, {"ath": ["AT1G58100"]}), ("TCP8[fc00258]", CLUSTER, {"ath": ["AT1G58100"]}),
               ("cZ", ["Node", "Metabolite"], {}))
    ns["TCP8,14,15[fc00332]"].short_name, ns["TCP8,14,15[fc00332]"].functional_cluster_id = "TCP8,14,15", "fc00332"
    ns["TCP8,14,15[fc00332]"].pathway = "P1"
    ns["TCP8[fc00258]"].short_name, ns["TCP8[fc00258]"].functional_cluster_id = "TCP8", "fc00258"
    r = reaction("rx3", "binding/oligomerisation", "inhibition",
                 [("TCP8,14,15[fc00332]", "protein", "nucleus", "SUBSTRATE"),
                  ("cZ", "metabolite", "nucleus", "SUBSTRATE"), ("TCP8|cZ", "complex", "nucleus", "PRODUCT")])
    r2 = reaction("rx4", "binding/oligomerisation", "inhibition",
                  [("TCP8[fc00258]", "protein", "nucleus", "SUBSTRATE"),
                   ("cZ", "metabolite", "nucleus", "SUBSTRATE"), ("TCP8|cZ", "complex", "nucleus", "PRODUCT")])
    net.create_gene_network(adapter(r, r2, nodes=ns), None, tmp_path / "n.tsv")
    gene = {n["id"]: n for n in read(tmp_path / "n.tsv")}["AT1G58100"]
    columns = ["display_label", "short_name", "pathway", "functional_cluster_id"]
    assert {c: gene[c].split("|") for c in columns} == {
        "display_label": ["TCP8,14,15", "TCP8"],
        "short_name": ["TCP8,14,15", "TCP8"], "pathway": ["P1", ""], "functional_cluster_id": ["fc00332", "fc00258"]}


def test_list_separator_in_a_value_warns(caplog):
    from pss_export.utils import clean
    assert clean(["RAPTOR2|TOR", "TORC1"], "TORC1, synonyms") == "RAPTOR2/TOR|TORC1"
    assert "TORC1, synonyms" in caplog.text and "RAPTOR2|TOR" in caplog.text
    assert clean(["PIP1;3", "PIP1C"]) == "PIP1;3|PIP1C"                           # gene symbols have ';'
    assert clean("free text; with a semicolon") == "free text; with a semicolon"     # not a list


def test_rank_0_on_all_edges(tmp_path):
    """ rank as in CKN: PSS edges are rank 0 (best supported) """
    r = reaction("rx1", "catalysis", "activation",
                 [("cZ", "metabolite", "cytoplasm", "SUBSTRATE"), ("cZ-P", "metabolite", "cytoplasm", "PRODUCT"),
                  ("ADK", "protein_active", "cytoplasm", "ACTIVATES")])
    net.create_interaction_network(adapter(r), tmp_path / "i.tsv", None)
    net.create_gene_network(adapter(r), tmp_path / "g.tsv", None)
    net.create_reaction_graph(adapter(r), tmp_path / "r.tsv", None)
    for name in ("i.tsv", "g.tsv", "r.tsv"):
        rows = read(tmp_path / name)
        assert rows and {e["rank"] for e in rows} == {"0"} and list(rows[0])[4] == "rank"


def test_abstract_clusters_are_like_metabolites(tmp_path):
    """ abstract clusters (PlantAbstract, no genes) don't decide a reaction's species, and stay named nodes in
    the gene network """
    ns = nodes(("ADK", CLUSTER, {"ath": ["AT1"]}),
               ("TNL[fc00500]", ["Node", "FunctionalCluster", "Plant", "PlantAbstract"], {}),
               ("cZ", ["Node", "Metabolite"], {}))
    ns["TNL[fc00500]"].functional_cluster_id = "fc00500"
    r = reaction("rx7", "protein activation", "activation",
                 [("ADK", "protein", "cytoplasm", "SUBSTRATE"), ("ADK", "protein_active", "cytoplasm", "PRODUCT"),
                  ("TNL[fc00500]", "protein_active", "cytoplasm", "ACTIVATES")])
    assert r.has_genes_in(ns, "ath") and not r.has_genes_in(ns, "stu")      # decided by ADK only
    assert r.gene_clusters(ns) == {"ADK"}
    net.create_gene_network(adapter(r, nodes=ns), tmp_path / "e.tsv", tmp_path / "n.tsv")
    assert {(e["source"], e["target"], e["source_type"]) for e in read(tmp_path / "e.tsv")} == {
        ("TNL[fc00500]", "AT1", "PlantAbstract")}
    rows = {n["id"]: n for n in read(tmp_path / "n.tsv")}
    assert (rows["TNL[fc00500]"]["node_type"], rows["TNL[fc00500]"]["functional_cluster_id"]) == ("PlantAbstract", "fc00500")


def test_complex_needs_genes_of_its_component_clusters():
    """ a reaction with a complex is in a species only if the complex's component clusters have genes there """
    ns = nodes(("GID", CLUSTER, {"ath": ["AT3G05120"], "osa": ["LOC_OS05G33730"]}),
               ("SLR1", CLUSTER, {"osa": ["LOC_OS06G03710"]}),
               ("GA|GID1|SLR1", ["Node", "Complex"], {}), ("GA", ["Node", "Metabolite"], {}))
    ns["GA|GID1|SLR1"].components = ["GA", "GID", "SLR1"]
    r = reaction("rx00298", "degradation/secretion", "activation",
                 [("GA|GID1|SLR1", "complex", "nucleus", "SUBSTRATE")])
    assert r.has_genes_in(ns, "osa") and not r.has_genes_in(ns, "ath")


def test_ignored_nodes_do_not_decide_the_species():
    """ the nodes to ignore (models) are as if they weren't there: an ignored cluster without genes in the
    species, or one among the components of a complex, doesn't leave the reaction out """
    ns = nodes(("HSP90", CLUSTER, {"stu": ["S1"]}), ("HSP70", CLUSTER, {"ath": ["A1"]}),
               ("HSP70|HSP90", ["Node", "Complex"], {}), ("COI1", CLUSTER, {"stu": ["S2"]}))
    ns["HSP70|HSP90"].components = ["HSP70", "HSP90"]
    forms = reaction("rx00729", "binding/oligomerisation", "activation",
                     [("HSP90", "protein", "cytoplasm", "SUBSTRATE"), ("HSP70", "protein", "cytoplasm", "SUBSTRATE"),
                      ("HSP70|HSP90", "complex", "cytoplasm", "PRODUCT")])
    uses = reaction("rx00730", "protein activation", "activation",
                    [("COI1", "protein", "cytoplasm", "SUBSTRATE"), ("COI1", "protein_active", "cytoplasm", "PRODUCT"),
                     ("HSP70|HSP90", "complex", "cytoplasm", "ACTIVATES")])
    for r in (forms, uses):
        assert not r.has_genes_in(ns, "stu") and r.has_genes_in(ns, "stu", ignore=["HSP70"])


def test_faidare(tmp_path):
    """ an entry per gene and species of the gene clusters in reactions, with the cluster's reactions, MapMan bins
    and a link to the cluster in the PSS Explorer; every species (no species filter) """
    import json
    from pss_export.faidare import create_faidare
    ns = nodes(("ADK", CLUSTER, {"ath": ["AT1", "AT2"], "stu": ["S1"]}), ("cZ", ["Node", "Metabolite"], {}),
               ("cZ-P", ["Node", "Metabolite"], {}))
    ns["ADK"].functional_cluster_id, ns["ADK"].short_name, ns["ADK"]._synonyms = "fc00001", "ADK", ["ADK", "ATADK"]
    ns["ADK"].mapman = ["26.11_External stimuli response.pathogen"]
    r = reaction("rx1", "catalysis", "activation",
                 [("cZ", "metabolite", "cytoplasm", "SUBSTRATE"), ("cZ-P", "metabolite", "cytoplasm", "PRODUCT"),
                  ("ADK", "protein_active", "cytoplasm", "ACTIVATES")])
    with pytest.raises(ValueError, match="every species"):
        create_faidare(adapter(r, nodes=ns), tmp_path / "f.json")
    assert create_faidare(adapter(r, nodes=ns, species=None), tmp_path / "f.json") == 3
    entries = json.load(open(tmp_path / "f.json"))
    assert [(e["name"], e["species"]) for e in entries] == [
        ("AT1", ["Arabidopsis thaliana"]), ("AT2", ["Arabidopsis thaliana"]), ("S1", ["Solanum tuberosum"])]
    e = entries[0]
    assert e["url"] == "https://skm.nib.si/pss/?functional_cluster_id=fc00001"
    assert "ADK takes part in catalysis with cZ, cZ-P. Synonyms are: ATADK. " in e["description"]
    assert (e["annotationId"], e["annotationName"]) == (
        ["MapMan4:26.11"], ["External stimuli response.pathogen (MapMan4:26.11)"])
