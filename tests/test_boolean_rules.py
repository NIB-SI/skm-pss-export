"""The Boolean rules of the models (pss_export/boolean/boolean.py), used by TabularQual and BoolNet: the rule of each
reaction type and effect (its target and inputs), which combinations have none, and how a species' rules from
several reactions are combined (rule_composer). No database needed."""
import re

import pytest

import pss_export.entity_classes as ec
import pss_export.pss.pss_reaction_definitions as rdef
from pss_export.boolean.boolean import reaction_rule_constructor, rule_composer

T = rdef.reaction_types


def species(name):
    s = ec.Species(name, "protein", "cytoplasm")
    s.set_id(f"s_{name}")
    return s


def reaction(reaction_type, effect, substrates=("A", "B"), products=("P",), modifiers=("M", "N")):
    r = ec.Reaction("rx00001", reaction_type, {"reaction_effect": effect})
    for name in substrates:
        r.add_substrate(species(name))
    for name in products:
        r.add_product(species(name))
    for name in modifiers:
        r.add_modifier(species(name))
    return r


def rule_of(r):
    constructor = reaction_rule_constructor(r)
    assert constructor is not None, (r.reaction_type, r.reaction_effect)
    return constructor(r)


SUBSTRATES, PRODUCTS = {"s_A", "s_B"}, {"s_P"}
ALL_INPUTS = "s_A & s_B & s_M & s_N"            # substrates and modifiers
MODIFIERS = "s_M & s_N"


# (type, effect, targets, rule) for a reaction A + B -> P, modified by M and N
RULES = [
    (T.BINDING_OLIGOMERISATION, "activation", PRODUCTS, ALL_INPUTS),
    (T.BINDING_OLIGOMERISATION, "inhibition", SUBSTRATES, f"!({ALL_INPUTS})"),
    (T.DISSOCIATION, "activation", PRODUCTS, ALL_INPUTS),
    (T.CATALYSIS, "activation", PRODUCTS, ALL_INPUTS),
    (T.PROTEIN_ACTIVATION, "activation", PRODUCTS, ALL_INPUTS),
    (T.TRANSCRIPTIONAL_TRANSLATIONAL_ACTIVATION, "activation", PRODUCTS, MODIFIERS),
    (T.TRANSLOCATION, "activation", PRODUCTS, ALL_INPUTS),
    (T.UNKNOWN, "activation", PRODUCTS, ALL_INPUTS),
    (T.DEGRADATION_SECRETION, "inhibition", SUBSTRATES, f"!({MODIFIERS})"),
    (T.PROTEIN_DEACTIVATION, "inhibition", SUBSTRATES, f"!({MODIFIERS})"),
    (T.TRANSCRIPTIONAL_TRANSLATIONAL_REPRESSION, "inhibition", PRODUCTS, f"!({MODIFIERS})"),
    (T.UNKNOWN, "inhibition", SUBSTRATES, f"!({MODIFIERS})"),
]


class TestReactionRules:
    """The rule each reaction type and effect gives: which species it sets (the targets) and from what."""

    @pytest.mark.parametrize("reaction_type, effect, targets, rule", RULES,
                             ids=[f"{t}-{e}" for t, e, _, _ in RULES])
    def test_rule(self, reaction_type, effect, targets, rule):
        assert rule_of(reaction(reaction_type, effect)) == (targets, rule)

    def test_every_reaction_type_has_a_rule(self):
        types = {getattr(T, k) for k in dir(T) if k.isupper()}
        assert {t for t, _, _, _ in RULES} == types

    def test_inputs_sorted(self):
        """the same reaction always gives the same rule"""
        r = reaction(T.CATALYSIS, "activation", substrates=("Z", "A"), modifiers=("M", "B"))
        assert rule_of(r)[1] == "s_A & s_B & s_M & s_Z"


class TestNoRule:
    """A reaction without inputs gives no rule (None), and is left out of the model's transitions."""

    @pytest.mark.parametrize("reaction_type, effect", [
        (T.TRANSCRIPTIONAL_TRANSLATIONAL_ACTIVATION, "activation"),
        (T.TRANSCRIPTIONAL_TRANSLATIONAL_REPRESSION, "inhibition"),
        (T.DEGRADATION_SECRETION, "inhibition"),
        (T.PROTEIN_DEACTIVATION, "inhibition"),
        (T.UNKNOWN, "inhibition"),
    ])
    def test_modifier_rules_need_a_modifier(self, reaction_type, effect):
        assert rule_of(reaction(reaction_type, effect, modifiers=()))[1] is None

    @pytest.mark.parametrize("reaction_type, effect", [
        (T.BINDING_OLIGOMERISATION, "activation"),
        (T.CATALYSIS, "activation"),
        (T.PROTEIN_ACTIVATION, "activation"),
        (T.TRANSLOCATION, "activation"),
    ])
    def test_rules_need_an_input(self, reaction_type, effect):
        assert rule_of(reaction(reaction_type, effect, substrates=(), modifiers=()))[1] is None

    @pytest.mark.parametrize("reaction_type, effect", [
        (T.CATALYSIS, "inhibition"),
        (T.DEGRADATION_SECRETION, "activation"),
        (T.PROTEIN_ACTIVATION, "inhibition"),
        (T.PROTEIN_DEACTIVATION, "activation"),
        (T.TRANSCRIPTIONAL_TRANSLATIONAL_REPRESSION, "activation"),
        (T.DISSOCIATION, "inhibition"),
        (T.TRANSLOCATION, "inhibition"),
        (T.CATALYSIS, None),
    ])
    def test_unsupported_combinations(self, reaction_type, effect, capsys):
        """no rule for these type and effect combinations: the reaction is reported and left out"""
        assert reaction_rule_constructor(reaction(reaction_type, effect)) is None
        assert "boolnet: issue with rx00001" in capsys.readouterr().out


def evaluate(rule, values):
    """a rule's value for the given species values (BoolNet syntax: ! before & before |, as Python's not/and/or)"""
    expression = re.sub(r"!", " not ", rule).replace("&", " and ").replace("|", " or ")
    return eval(expression, {}, values)


class TestRuleComposer:
    """A species' rule from the rules of all reactions that set it: any activation, unless inhibited."""

    def test_no_rules(self):
        assert rule_composer("s_X", [], []) is None

    def test_activations_any(self):
        assert rule_composer("s_X", ["s_A", "s_B & s_C"], []) == "s_A | s_B & s_C"

    def test_one_activation_one_inhibition(self):
        assert rule_composer("s_X", ["s_A"], ["!(s_M)"]) == "s_A & !(s_M)"

    def test_inhibitions_all(self):
        assert rule_composer("s_X", ["s_A"], ["!(s_M)", "!(s_N)"]) == "s_A & ( !(s_M) & !(s_N) )"

    def test_only_inhibited_keeps_its_value(self):
        """without an activating reaction the species keeps its value until inhibited (it is an input)"""
        rule = rule_composer("s_X", [], ["!(s_M)"])
        assert rule == "s_X & !(s_M)"
        assert evaluate(rule, {"s_X": True, "s_M": False}) and not evaluate(rule, {"s_X": True, "s_M": True})

    @pytest.mark.xfail(strict=True, reason="the activations are joined with | and the inhibition appended with & "
                       "without brackets: & binds tighter, so the inhibition only applies to the last activation "
                       "(14 of 156 inhibited rules in the public BoolNet export of 2026-10-20, e.g. EIN3)")
    def test_inhibition_applies_to_every_activation(self):
        rule = rule_composer("s_X", ["s_A", "s_B"], ["!(s_M)"])
        # A is on, but M inhibits: X must be off whichever activation is on
        assert not evaluate(rule, {"s_A": True, "s_B": False, "s_M": True})
