"""
Tests for pss-export Python API (no database connection required).

Covers:
  - Config loading
  - Species / Reaction / IDTracker entity classes
  - Reaction subtype and role assignment
  - AnnotationManager reference processing
  - PSSCollector access levels and nodes_to_ignore
  - Connection settings and connection lifecycle
"""

import os
import re
import tempfile
import pytest

import pss_export.pss.config
import pss_export.pss.pss_reaction_definitions as rdef
import pss_export.entity_classes as ec
from pss_export.pss.config import Config, pss_export_config
from pss_export.annotations.annotation_manager import AnnotationManager
from pss_export.pss.collectors import PSSCollector
from pss_export.graph_db import resolve_connection_settings
import pss_export.pss.pss_adapter as pss_adapter_module
from pss_export import PSSAdapter


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

class TestConfig:

    def test_kwargs_override(self):
        """Config accepts arbitrary keyword arguments as attributes."""
        c = Config(filename=None, foo="bar", count=42)
        assert c.foo == "bar"
        assert c.count == 42

    def test_pss_export_config_loaded(self):
        """pss_export_config loads the YAML and exposes expected keys."""
        assert isinstance(pss_export_config.nodes_to_ignore, list)
        assert len(pss_export_config.nodes_to_ignore) > 0
        assert hasattr(pss_export_config, "compartment_to_short")
        assert hasattr(pss_export_config, "node_form_to_SBO")
        assert hasattr(pss_export_config, "reaction_subtype_to_SBO")

    def test_compartment_short_codes(self):
        """Key compartments map to their expected short codes."""
        m = pss_export_config.compartment_to_short
        assert m["cytoplasm"] == "cyt"
        assert m["nucleus"] == "nuc"
        assert m["chloroplast"] == "chl"


# ---------------------------------------------------------------------------
# Species
# ---------------------------------------------------------------------------

class TestSpecies:

    def test_basic_construction(self):
        s = ec.Species("WRKY33", "protein", "nucleus")
        assert s.name == "WRKY33"
        assert s.form == "protein"
        assert s.compartment == "nucleus"

    def test_gene_is_constant(self):
        assert ec.Species("WRKY33", "gene", "nucleus").constant is True

    def test_non_gene_not_constant(self):
        assert ec.Species("WRKY33", "protein", "nucleus").constant is False

    def test_none_compartment_defaults_to_cytoplasm(self):
        assert ec.Species("X", "protein", None).compartment == "cytoplasm"

    def test_unknown_compartment_defaults_to_cytoplasm(self):
        assert ec.Species("X", "protein", "unknown").compartment == "cytoplasm"

    def test_putative_prefix_stripped(self):
        assert ec.Species("X", "protein", "putative:nucleus").compartment == "nucleus"

    def test_sbo_term_set_for_known_form(self):
        s = ec.Species("X", "protein", "cytoplasm")
        assert s.sbo_term == 297  # SBO:0000297 — protein complex

    def test_sbo_term_none_for_unknown_form(self):
        s = ec.Species("X", "unknown_form_xyz", "cytoplasm")
        assert s.sbo_term is None

    def test_label_strips_locus_tag(self):
        s = ec.Species("WRKY33[AT2G38470]", "protein", "cytoplasm")
        assert s.label == "WRKY33"


# ---------------------------------------------------------------------------
# Reaction subtype and role assignment
# ---------------------------------------------------------------------------

class TestReaction:

    def _make(self, rtype, props=None):
        return ec.Reaction("rx", rtype, props or {})

    def test_catalysis_no_participants(self):
        r = self._make(rdef.reaction_types.CATALYSIS)
        assert r.reaction_subtype == rdef.reaction_subtypes.CATALYSIS_WITHOUT_SUBSTRATE_AND_WITHOUT_MODIFIER
        assert r.modifier_role == rdef.participant_roles.CATALYST

    def test_binding_without_modifier(self):
        r = self._make(rdef.reaction_types.BINDING_OLIGOMERISATION)
        assert r.reaction_subtype == rdef.reaction_subtypes.BINDING_WITHOUT_MODIFIER
        assert r.substrate_role == rdef.participant_roles.INTERACTOR

    def test_translocation_without_modifier(self):
        r = self._make(rdef.reaction_types.TRANSLOCATION)
        assert r.reaction_subtype == rdef.reaction_subtypes.TRANSLOCATION_WITHOUT_MODIFIER
        assert r.modifier_role == rdef.participant_roles.TRANSPORTER

    def test_transcriptional_activation(self):
        r = self._make(
            rdef.reaction_types.TRANSCRIPTIONAL_TRANSLATIONAL_ACTIVATION,
            {"reaction_mechanism": "transcription"},
        )
        assert r.reaction_subtype == rdef.reaction_subtypes.TRANSCRIPTIONAL_ACTIVATION_WITH_MODIFIER
        assert r.modifier_role == rdef.participant_roles.STIMULATOR
        assert r.substrate_role == rdef.participant_roles.TEMPLATE

    def test_transcriptional_repression(self):
        r = self._make(
            rdef.reaction_types.TRANSCRIPTIONAL_TRANSLATIONAL_REPRESSION,
            {"reaction_mechanism": "transcription"},
        )
        assert r.reaction_subtype == rdef.reaction_subtypes.TRANSCRIPTIONAL_REPRESSION_WITH_MODIFIER
        assert r.modifier_role == rdef.participant_roles.INHIBITOR

    def test_unknown_activation(self):
        r = self._make(rdef.reaction_types.UNKNOWN, {"reaction_effect": "activation"})
        assert r.reaction_subtype == rdef.reaction_subtypes.UNKNOWN_ACTIVATION_WITHOUT_MODIFIER

    def test_sbo_term_assigned(self):
        r = self._make(rdef.reaction_types.CATALYSIS)
        assert r.reaction_type_sbo == 176  # SBO:0000176 — biochemical reaction
        assert r.kinetic_law_sbo == 29     # SBO:0000029 — Henri-Michaelis-Menten

    def test_add_substrate_product_modifier(self):
        r = self._make(rdef.reaction_types.CATALYSIS)
        r.add_substrate(ec.Species("MPK4", "protein", "cytoplasm"))
        r.add_product(ec.Species("WRKY33", "protein_active", "nucleus"))
        r.add_modifier(ec.Species("MPK3", "protein_active", "cytoplasm"))
        assert r.has_substrates()
        assert r.has_modifiers()
        assert len(r.substrates) == 1
        assert len(r.products) == 1
        assert len(r.modifiers) == 1

    def test_subtype_follows_participants(self):
        """Subtype describes the reaction as it currently is, including
        participants added after construction (e.g. by model fixes)."""
        r = self._make(rdef.reaction_types.CATALYSIS)
        assert r.reaction_subtype == rdef.reaction_subtypes.CATALYSIS_WITHOUT_SUBSTRATE_AND_WITHOUT_MODIFIER
        r.add_substrate(ec.Species("X", "protein", "cytoplasm"))
        assert r.reaction_subtype == rdef.reaction_subtypes.CATALYSIS_WITH_SUBSTRATE_AND_WITHOUT_MODIFIER
        r.add_modifier(ec.Species("E", "protein_active", "cytoplasm"))
        assert r.reaction_subtype == rdef.reaction_subtypes.CATALYSIS_WITH_SUBSTRATE_AND_WITH_MODIFIER
        r.substrates.clear()
        assert r.reaction_subtype == rdef.reaction_subtypes.CATALYSIS_WITHOUT_SUBSTRATE_AND_WITH_MODIFIER

    def test_sbo_follows_participants(self):
        """Translocation with a transporter is active transport, without passive."""
        r = self._make(rdef.reaction_types.TRANSLOCATION)
        assert r.reaction_type_sbo == 658  # SBO:0000658 — passive transport
        r.add_modifier(ec.Species("T", "protein_active", "cytoplasm"))
        assert r.reaction_type_sbo == 657  # SBO:0000657 — active transport

    def test_dissociation_sbo_with_and_without_modifier(self):
        """Both variants are SBO:0000180 dissociation (0000015 is the
        participant role 'substrate', not a process)."""
        r = self._make(rdef.reaction_types.DISSOCIATION)
        assert r.reaction_type_sbo == 180
        r.add_modifier(ec.Species("M", "protein_active", "cytoplasm"))
        assert r.reaction_type_sbo == 180

    def test_subtype_is_read_only(self):
        r = self._make(rdef.reaction_types.CATALYSIS)
        with pytest.raises(AttributeError):
            r.reaction_subtype = "something"

    def test_every_subtype_has_sbo_terms(self):
        """Every subtype has both SBO terms configured (else they'd be None)."""
        configured = pss_export_config.reaction_subtype_to_SBO
        missing = [s for s in rdef.ALL_REACTION_SUBTYPES if s not in configured]
        assert missing == []
        for subtype in rdef.ALL_REACTION_SUBTYPES:
            assert int(configured[subtype]["reaction_type_SBO"]) > 0
            assert int(configured[subtype]["kinetic_law_SBO"]) > 0


# ---------------------------------------------------------------------------
# IDTracker
# ---------------------------------------------------------------------------

class TestIDTracker:

    def test_new_species_gets_unique_id(self):
        tracker = ec.IDTracker()
        s = ec.Species("WRKY33", "protein", "nucleus")
        id_, status = tracker.get_species_id(s)
        assert status == 0  # new
        assert "WRKY33" in id_
        assert "nuc" in id_

    @pytest.mark.parametrize("name", ["WRKY33", "12-OH-JA-Ile", "4CL", "6K2"])
    def test_species_id_is_a_valid_sbml_id(self, name):
        ''' SBML ids start with a letter or _ (names like 4CL don't) '''
        id_, _ = ec.IDTracker().get_species_id(ec.Species(name, "protein", "cytoplasm"))
        assert id_.startswith("s_")
        assert re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", id_)

    def test_same_species_returns_existing_id(self):
        tracker = ec.IDTracker()
        s = ec.Species("WRKY33", "protein", "nucleus")
        id1, _ = tracker.get_species_id(s)
        tracker.set_species_id(s, id1)
        id2, status = tracker.get_species_id(s)
        assert status == 1  # existing
        assert id1 == id2

    def test_different_form_gets_different_id(self):
        tracker = ec.IDTracker()
        s1 = ec.Species("WRKY33", "protein", "nucleus")
        s2 = ec.Species("WRKY33", "protein_active", "nucleus")
        id1, _ = tracker.get_species_id(s1)
        tracker.set_species_id(s1, id1)
        id2, _ = tracker.get_species_id(s2)
        assert id1 != id2

    def test_counters_increment(self):
        tracker = ec.IDTracker()
        s = ec.Species("X", "protein", "cytoplasm")
        id_, _ = tracker.get_species_id(s)
        tracker.set_species_id(s, id_)
        assert tracker.counters["species"] == 1

    def test_write_entities_table(self, tmp_path):
        tracker = ec.IDTracker()

        s1 = ec.Species("MPK4", "protein", "cytoplasm")
        id1, _ = tracker.get_species_id(s1)
        tracker.set_species_id(s1, id1)

        s2 = ec.Species("WRKY33", "protein_active", "nucleus")
        id2, _ = tracker.get_species_id(s2)
        tracker.set_species_id(s2, id2)

        r = ec.Reaction("r1", rdef.reaction_types.CATALYSIS, {})
        tracker.set_reaction_id(r, "r_cat_1")

        outfile = str(tmp_path / "entities.tsv")
        tracker.write_entities_table(outfile)

        with open(outfile) as f:
            lines = f.readlines()

        assert lines[0].strip() == "id\ttype\tname\tform\tcompartment"
        assert len(lines) == 4  # header + 2 species + 1 reaction
        ids = {line.split("\t")[0] for line in lines[1:]}
        assert id1 in ids
        assert id2 in ids
        assert "r_cat_1" in ids


# ---------------------------------------------------------------------------
# AnnotationManager
# ---------------------------------------------------------------------------

class DummyStrategy:
    """Passthrough strategy that returns valid_records for inspection."""
    def format_node(self, valid_records):
        return valid_records


class TestAnnotationManager:

    @pytest.fixture
    def am(self):
        manager = AnnotationManager()
        manager.register_export_strategy("test", DummyStrategy())
        return manager

    def test_registry_loaded(self, am):
        assert "uniprot" in am._registry
        assert "chebi" in am._registry

    def test_known_ref_is_valid(self, am):
        result, invalid = am.process_node("test", ["uniprot:P12345"])
        assert len(result) == 1
        assert invalid == []
        assert result[0]["local_id"] == "P12345"
        assert result[0]["qualifier"] == "bqbiol:is"

    def test_ref_without_colon_is_invalid(self, am):
        result, invalid = am.process_node("test", ["BADREF"])
        assert result == []
        assert "BADREF" in invalid

    def test_unknown_prefix_is_invalid(self, am):
        result, invalid = am.process_node("test", ["notadb:XYZ"])
        assert result == []
        assert "notadb:XYZ" in invalid

    def test_mixed_refs(self, am):
        refs = ["uniprot:P12345", "NOCODON", "chebi:12345", "ghost:000"]
        result, invalid = am.process_node("test", refs)
        assert len(result) == 2
        assert len(invalid) == 2

    def test_unsupported_format_raises(self, am):
        with pytest.raises(ValueError, match="Unsupported format"):
            am.process_node("nonexistent_format", ["uniprot:P12345"])

# ---------------------------------------------------------------------------
# PSSCollector access levels
# ---------------------------------------------------------------------------

class TestCollectorAccess:

    @pytest.mark.parametrize("access, expected", [
        ("public", "public"),
        ("restricted", "restricted"),
        ("all", "restricted"),   # legacy alias
    ])
    def test_access_levels(self, access, expected):
        collector = PSSCollector(None, access=access)
        assert collector.access == expected

    def test_default_is_public(self):
        assert PSSCollector(None).access == "public"

    @pytest.mark.parametrize("access", ["Restricted", "private", "", None])
    def test_invalid_access_raises(self, access):
        with pytest.raises(ValueError, match="Invalid access"):
            PSSCollector(None, access=access)

    def test_public_filters_external_links(self):
        where, args = PSSCollector(None, access="public", nodes_to_ignore=None)._build_where_clause()
        assert "external_links" in where
        assert "invented_reason_allowlist" in args

    def test_pathway_filter_selects_whole_reactions(self):
        """The pathway filter selects reactions (EXISTS on a participant);
        it must not filter the collected edges (n)."""
        collector = PSSCollector(None, pathways=["Hormone - Abscisic acid (ABA)"], nodes_to_ignore=None)
        where, args = collector._build_where_clause()
        assert "EXISTS" in where
        assert "n." not in where
        assert args["pathways"] == ["Hormone - Abscisic acid (ABA)"]

    def test_restricted_has_no_external_links_filter(self):
        where, args = PSSCollector(None, access="restricted", nodes_to_ignore=None)._build_where_clause()
        assert "external_links" not in where
        assert "invented_reason_allowlist" not in args


# ---------------------------------------------------------------------------
# PSSCollector nodes_to_ignore
# ---------------------------------------------------------------------------

class FakeEdge(dict):
    """Stands in for a neo4j Relationship (type, start/end nodes, properties)."""
    def __init__(self, type_, start_node, end_node, **props):
        super().__init__(**props)
        self.type = type_
        self.start_node = start_node
        self.end_node = end_node


class FakePath:
    """Stands in for a neo4j Path p=(r)-[]-(n): end_node is the participant."""
    def __init__(self, edge, participant):
        self.relationships = [edge]
        self.end_node = participant


REACTION_NODE = {"name": "rx"}


def substrate(name):
    node = {"name": name}
    return FakePath(FakeEdge("SUBSTRATE", node, REACTION_NODE,
                             source_location="cytoplasm", source_form="protein"), node)


def product(name):
    node = {"name": name}
    return FakePath(FakeEdge("PRODUCT", REACTION_NODE, node,
                             target_location="cytoplasm", target_form="protein"), node)


def modifier(name):
    node = {"name": name}
    return FakePath(FakeEdge("ACTIVATES", node, REACTION_NODE,
                             source_location="cytoplasm", source_form="protein_active"), node)


class TestCollectorIgnoreNodes:

    @pytest.fixture
    def collector(self):
        return PSSCollector(None, nodes_to_ignore=["SCF"])

    def _build(self, collector, paths, reaction_type=rdef.reaction_types.BINDING_OLIGOMERISATION):
        return collector._build_reaction("rx00001", {"reaction_type": reaction_type}, paths)

    def test_nodes_to_ignore_not_in_cypher(self, collector):
        where, args = collector._build_where_clause()
        assert "nodes_to_ignore" not in where
        assert "nodes_to_ignore" not in args

    def test_other_substrate_left_keeps_reaction(self, collector):
        r = self._build(collector, [substrate("A"), substrate("SCF"), product("A|SCF")])
        assert [s.name for s in r.substrates] == ["A"]
        assert [p.name for p in r.products] == ["A|SCF"]

    def test_only_substrate_ignored_drops_reaction(self, collector):
        assert self._build(collector, [substrate("SCF"), product("B")]) is None

    def test_only_product_ignored_drops_reaction(self, collector):
        assert self._build(collector, [substrate("A"), product("SCF")]) is None

    def test_modifier_ignored_keeps_reaction(self, collector):
        r = self._build(collector, [substrate("X"), product("Y"), modifier("SCF")],
                        rdef.reaction_types.CATALYSIS)
        assert r.modifiers == []
        assert len(r.substrates) == 1 and len(r.products) == 1

    def test_all_edges_ignored_drops_reaction(self, collector):
        assert self._build(collector, [modifier("SCF")], rdef.reaction_types.CATALYSIS) is None

    def test_empty_side_by_design_keeps_reaction(self, collector):
        """Gene substrate is skipped (include_genes=False), so an empty
        substrate side is not caused by ignoring and must not drop it."""
        r = self._build(collector, [substrate("GENE"), product("P"), modifier("SCF")],
                        rdef.reaction_types.TRANSCRIPTIONAL_TRANSLATIONAL_ACTIVATION)
        assert r is not None
        assert r.substrates == []
        assert [p.name for p in r.products] == ["P"]

    def test_nothing_ignored(self):
        collector = PSSCollector(None, nodes_to_ignore=None)
        r = self._build(collector, [substrate("SCF"), product("B")])
        assert [s.name for s in r.substrates] == ["SCF"]


# ---------------------------------------------------------------------------
# Connection settings
# ---------------------------------------------------------------------------

ENV_VARS = ("MY_NEO4J_URI", "MY_NEO4J_USER", "MY_NEO4J_PASSWORD")


@pytest.fixture
def clean_env(monkeypatch, tmp_path):
    """No connection settings in the environment, and no .env file."""
    for v in ENV_VARS:
        monkeypatch.delenv(v, raising=False)
    monkeypatch.chdir(tmp_path)
    return tmp_path


class TestConnectionSettings:

    def test_arguments(self, clean_env):
        s = resolve_connection_settings(uri="bolt://a:7687", user="u", pwd="p")
        assert s == {"uri": "bolt://a:7687", "user": "u", "pwd": "p"}

    def test_empty_strings_are_valid(self, clean_env):
        s = resolve_connection_settings(uri="bolt://a:7687", user="", pwd="")
        assert s["user"] == "" and s["pwd"] == ""

    def test_environment(self, clean_env, monkeypatch):
        monkeypatch.setenv("MY_NEO4J_URI", "bolt://env:7687")
        monkeypatch.setenv("MY_NEO4J_USER", "env_user")
        monkeypatch.setenv("MY_NEO4J_PASSWORD", "env_pwd")
        s = resolve_connection_settings()
        assert s == {"uri": "bolt://env:7687", "user": "env_user", "pwd": "env_pwd"}

    def test_dotenv_file(self, clean_env):
        (clean_env / ".env").write_text(
            "MY_NEO4J_URI=bolt://file:7687\nMY_NEO4J_USER=file_user\nMY_NEO4J_PASSWORD=file_pwd\n")
        s = resolve_connection_settings()
        assert s == {"uri": "bolt://file:7687", "user": "file_user", "pwd": "file_pwd"}

    def test_precedence(self, clean_env, monkeypatch):
        """argument > environment > .env file"""
        (clean_env / ".env").write_text(
            "MY_NEO4J_URI=bolt://file:7687\nMY_NEO4J_USER=file_user\nMY_NEO4J_PASSWORD=file_pwd\n")
        monkeypatch.setenv("MY_NEO4J_USER", "env_user")
        s = resolve_connection_settings(uri="bolt://arg:7687")
        assert s == {"uri": "bolt://arg:7687", "user": "env_user", "pwd": "file_pwd"}

    def test_missing_raises(self, clean_env):
        with pytest.raises(ValueError, match="MY_NEO4J_USER, MY_NEO4J_PASSWORD"):
            resolve_connection_settings(uri="bolt://a:7687")

    def test_adapter_fails_early_on_missing_settings(self, clean_env):
        with pytest.raises(ValueError, match="Missing database connection settings"):
            PSSAdapter()


# ---------------------------------------------------------------------------
# Connection lifecycle
# ---------------------------------------------------------------------------

class FakeGraphDB:
    """Records connections; returns no data, or raises on query if asked to."""
    instances = []

    def __init__(self, uri, user, pwd, fail=False):
        self.closed = False
        self.fail = fail
        FakeGraphDB.instances.append(self)

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.closed = True

    def run_query(self, query_function, *args):
        if self.fail:
            raise RuntimeError("query failed")
        return []


class TestConnectionLifecycle:

    @pytest.fixture(autouse=True)
    def fake_db(self, monkeypatch):
        FakeGraphDB.instances = []
        monkeypatch.setattr(pss_adapter_module, "GraphDB", FakeGraphDB)

    def _adapter(self):
        return PSSAdapter(neo4j_uri="bolt://x:7687", neo4j_user="", neo4j_password="")

    def test_no_connection_at_construction(self):
        self._adapter()
        assert FakeGraphDB.instances == []

    def test_connection_closed_after_collect(self):
        self._adapter().collect_reactions(nodes_to_ignore=None)
        assert len(FakeGraphDB.instances) == 1
        assert FakeGraphDB.instances[0].closed

    def test_connection_closed_on_error(self, monkeypatch):
        monkeypatch.setattr(pss_adapter_module, "GraphDB",
                            lambda **kw: FakeGraphDB(**kw, fail=True))
        with pytest.raises(RuntimeError, match="query failed"):
            self._adapter().collect_reactions(nodes_to_ignore=None)
        assert FakeGraphDB.instances[0].closed


class TestModelFixTransport:
    """Transport reactions added by the model fixes get a Boolean rule."""

    def test_transport_reaction_has_rule(self):
        pytest.importorskip("matplotlib")
        from types import SimpleNamespace
        from pss_export.model_fixes.model_fixes import ModelFixer, TransportReaction
        from pss_export.boolean.boolean import reaction_rule_constructor

        adapter = SimpleNamespace(reactions={}, additional_reactions=[], include_genes=False)
        fix = TransportReaction("AREB/ABF[AT1G45249]", "protein_active", "cytoplasm", "nucleus")
        assert ModelFixer(adapter).apply_model_fixes([fix]) == 1

        reaction = adapter.reactions[adapter.additional_reactions[0]]
        assert reaction.reaction_type == rdef.reaction_types.TRANSLOCATION
        assert reaction.reaction_effect == "activation"
        assert reaction_rule_constructor(reaction) is not None
