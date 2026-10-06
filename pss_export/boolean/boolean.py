# class to generate the boolean rule based in the reaction type
# class BooleanRules:
''' Functions to create Boolean rules from PSS reaction types
All functions take three arguments:
s : set of substrate species names (strings)
p : set of product species names (strings)
m : set of modifier species names (strings)
Each function returns a tuple:
(target_species_name (string), boolean_rule (string))
The species in a rule are sorted, so that the same reactions always give the same rule.
'''

def generate_transition(rules):
    """
    List of rules
    """

    activation_rules = rule_dict["activation"]
    inhibition_rules = rule_dict["inhibition"]

    update_function = rule_composer(species_id, activation_rules, inhibition_rules)



def rule_composer(species_id, activation_rules, inhibition_rules):

    if len(activation_rules) + len(inhibition_rules) > 0:

        if len(activation_rules) == 0:
            activation_rules.append(species_id)
            # activation_reactions.append('self')

        update_function = f"{' | '.join(activation_rules)}"
        # update_reaction = f"( {' | '.join(activation_reactions)} )"

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
        return function(s_set, p_set, m_set)

    return wrapper


@preparerule
def binding_oligomerisation(s, p, m):
    target = p
    factors = s.union(m)
    if len(factors) == 0:
        return target, None
    rule = ' & '.join(sorted(factors))
    return target, rule

@preparerule
def binding_oligomerisation_inhibition(s, p, m):
    target = s
    factors = s.union(m)
    if len(factors) == 0:
        return target, None
    rule = f"!({' & '.join(sorted(factors))})"
    return target, rule

@preparerule
def dissociation(s, p, m):
    target = p
    factors = s.union(m)
    if len(factors) == 0:
        return target, None
    rule = ' & '.join(sorted(factors))
    return target, rule

@preparerule
def catalysis(s, p, m):
    target = p
    factors = s.union(m)
    if len(factors) == 0:
        return target, None
    rule = ' & '.join(sorted(factors))
    return target, rule

@preparerule
def protein_activation(s, p, m):
    target = p
    factors = s.union(m)
    if len(factors) == 0:
        return target, None
    rule = ' & '.join(sorted(factors))
    return target, rule

@preparerule
def transcriptional_translational_activation(s, p, m):
    target = p
    if len(m) == 0:
        return target, None
    rule = ' & '.join(sorted(m))
    return target, rule

@preparerule
def translocation(s, p, m):
    target = p
    factors = s.union(m)
    if len(factors) == 0:
        return target, None
    rule = ' & '.join(sorted(factors))
    return target, rule

@preparerule
def unknown_activation(s, p, m):
    target = p
    factors = s.union(m)
    if len(factors) == 0:
        return target, None
    rule = ' & '.join(sorted(factors))
    return target, rule

@preparerule
def degradation_secretion(s, p, m):
    target = s
    if len(m) == 0:
        return target, None
    rule = f"!({' & '.join(sorted(m))})"
    return target, rule

@preparerule
def protein_deactivation(s, p, m):
    target = s
    if len(m) == 0:
        return target, None
    rule = f"!({' & '.join(sorted(m))})"
    return target, rule

@preparerule
def transcriptional_translational_repression(s, p, m):
    target = p
    if len(m) == 0:
        return target, None
    rule = f"!({' & '.join(sorted(m))})"
    return target, rule

@preparerule
def unknown_inhibition(s, p, m):
    target = s
    if len(m) == 0:
        return target, None
    rule = f"!({' & '.join(sorted(m))})"
    return target, rule
