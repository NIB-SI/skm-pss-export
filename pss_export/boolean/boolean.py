# class to generate the boolean rule based in the reaction type
# class BooleanRules:
''' Functions to create Boolean rules from PSS reaction types
All functions take four arguments:
s : set of substrate species ids (strings)
p : set of product species ids (strings)
m : set of modifier species ids (strings)
processes : set of the reaction's process species ids (strings)
Each function returns a dict of the species the reaction sets (target ids) and their rule:
{target_species_id (string): boolean_rule (string)}, empty without a rule (no inputs).
The species in a rule are sorted, so that the same reactions always give the same rule.
'''

def rule_composer(species_id, activation_rules, inhibition_rules):
    ''' A species' rule from the rules of the reactions that set it: any activation, and none of the inhibitions
    (without an activation, the species keeps its value until inhibited). None without rules. '''

    if len(activation_rules) + len(inhibition_rules) > 0:

        # (copies: the callers' lists are kept as they are; the same rule from several reactions only once)
        activation_rules = list(dict.fromkeys(activation_rules)) or [species_id]
        inhibition_rules = list(dict.fromkeys(inhibition_rules))

        update_function = f"{' | '.join(activation_rules)}"
        # several activations before an inhibition: in brackets, as & binds tighter than |
        # (else the inhibition would only apply to the last activation)
        if len(activation_rules) > 1 and inhibition_rules:
            update_function = f"( {update_function} )"

        if len(inhibition_rules) == 0:
            pass
        elif len(inhibition_rules) == 1:
            update_function += f" & {inhibition_rules[0]}"
        else:
            update_function += f" & ( {' & '.join(inhibition_rules)} )"

        return update_function

def reaction_rule_constructor(reaction):

    rule_constructor = {
        'activation': {
            'binding/oligomerisation': binding_oligomerisation,
            'dissociation': dissociation,
            'catalysis': catalysis,
            'protein activation': protein_activation,
            'transcriptional/translational activation':
            transcriptional_translational_activation,
            'translocation': translocation,
            'unknown': unknown_activation,
        },
        'inhibition': {
            'binding/oligomerisation': binding_oligomerisation_inhibition,
            'degradation/secretion': degradation_secretion,
            'protein deactivation': protein_deactivation,
            'transcriptional/translational repression':
            transcriptional_translational_repression,
            'unknown': unknown_inhibition
        }
    }

    reaction_type = reaction.reaction_type
    reaction_effect = reaction.reaction_effect

    if not (reaction_effect in rule_constructor) or\
       not (reaction_type in rule_constructor[reaction_effect]):

        print(f"boolnet: issue with {reaction.reaction_id} "
              f"reaction effect: {reaction_effect} "
              f"reaction_type: {reaction_type}")
        return None

    return rule_constructor[reaction_effect][reaction_type]


def preparerule(function):

    def wrapper(reaction):
        # print(f"Processing reaction {reaction.reaction_id}...")
        # print(f"  substrates: {[n for n in reaction.substrates]}")
        # print(f"  products: {[n for n in reaction.products]}")
        # print(f"  modifiers: {[n for n in reaction.modifiers]}")
        s_set = {n.id for n in reaction.substrates}
        p_set = {n.id for n in reaction.products}
        m_set = {n.id for n in reaction.modifiers}
        processes = {n.id for n in reaction.substrates + reaction.products + reaction.modifiers
                     if n.form in ('process', 'process_active')}
        return function(s_set, p_set, m_set, processes)

    return wrapper


def rules_for(targets, factors, negate=False, processes=frozenset()):
    ''' A rule per target: all its inputs (factors), or none of them (negate); no rule without inputs.
    A process is no input of its own rule: it has no inactive form, so a regulated process is its own substrate and
    product (Thermotolerance + HSFA6b -> Thermotolerance: Thermotolerance = HSFA6b). Other species keep themselves
    (Phytophthora + cholesterol -> Phytophthora: it persists, given cholesterol). '''
    rules = {}
    for target in targets:
        inputs = factors - {target} if target in processes else factors
        if len(inputs) > 0:
            rule = ' & '.join(sorted(inputs))
            rules[target] = f"!({rule})" if negate else rule
    return rules


@preparerule
def binding_oligomerisation(s, p, m, processes):
    return rules_for(p, s | m, processes=processes)

@preparerule
def binding_oligomerisation_inhibition(s, p, m, processes):
    # each binding partner is off while the others (and the modifiers) are on: a rule per partner, without itself
    rules = {}
    for target in s:
        rules.update(rules_for({target}, (s - {target}) | m, negate=True, processes=processes))
    return rules

@preparerule
def dissociation(s, p, m, processes):
    return rules_for(p, s | m, processes=processes)

@preparerule
def catalysis(s, p, m, processes):
    return rules_for(p, s | m, processes=processes)

@preparerule
def protein_activation(s, p, m, processes):
    return rules_for(p, s | m, processes=processes)

@preparerule
def transcriptional_translational_activation(s, p, m, processes):
    return rules_for(p, m, processes=processes)

@preparerule
def translocation(s, p, m, processes):
    return rules_for(p, s | m, processes=processes)

@preparerule
def unknown_activation(s, p, m, processes):
    return rules_for(p, s | m, processes=processes)

@preparerule
def degradation_secretion(s, p, m, processes):
    return rules_for(s, m, negate=True, processes=processes)

@preparerule
def protein_deactivation(s, p, m, processes):
    return rules_for(s, m, negate=True, processes=processes)

@preparerule
def transcriptional_translational_repression(s, p, m, processes):
    return rules_for(p, m, negate=True, processes=processes)

@preparerule
def unknown_inhibition(s, p, m, processes):
    return rules_for(s, m, negate=True, processes=processes)
