'''
This module defines constants for different types of reactions in PSS.
'''


class reaction_types:
    """
    Class to hold constants for different types of reactions.
    """

    BINDING_OLIGOMERISATION = 'binding/oligomerisation'
    DISSOCIATION = 'dissociation'
    CATALYSIS = 'catalysis'
    DEGRADATION_SECRETION = 'degradation/secretion'
    PROTEIN_DEACTIVATION = 'protein deactivation'
    PROTEIN_ACTIVATION = 'protein activation'
    TRANSCRIPTIONAL_TRANSLATIONAL_ACTIVATION = 'transcriptional/translational activation'
    TRANSCRIPTIONAL_TRANSLATIONAL_REPRESSION = 'transcriptional/translational repression'
    TRANSLOCATION = 'translocation'
    UNKNOWN = 'unknown'


class reaction_classes:
    """Class to hold constants for different classes of reaction types.
    These classes are used to categorize reactions based on their biological function.
    """

    TRANSCRIPTIONAL_TRANSLATIONAL = [
        reaction_types.TRANSCRIPTIONAL_TRANSLATIONAL_ACTIVATION,
        reaction_types.TRANSCRIPTIONAL_TRANSLATIONAL_REPRESSION
    ]

class participant_roles:
    """
    Class to hold constants for different roles of participants in reactions.
    """
    SUBSTRATE = 'substrate'
    PRODUCT = 'product'
    MODIFIER = 'modifier'
    INTERACTOR = 'interactor'
    TEMPLATE = 'template'
    STIMULATOR = 'stimulator'
    INHIBITOR = 'inhibitor'
    CATALYST = 'catalyst'
    TRANSPORTER = 'transporter'

    @staticmethod
    def assign_roles(reaction):
        """
        Assign roles to reaction participants based on the reaction type.
        This method sets the roles of substrates, products, and modifiers
        according to predefined rules for each reaction type.
        """

        # reactant roles
        if reaction.reaction_type == reaction_types.BINDING_OLIGOMERISATION:
            reaction.substrate_role = participant_roles.INTERACTOR
        elif reaction.reaction_type in reaction_classes.TRANSCRIPTIONAL_TRANSLATIONAL:
            reaction.substrate_role = participant_roles.TEMPLATE
        else:
            reaction.substrate_role = participant_roles.SUBSTRATE

        # modifier roles
        if reaction.reaction_type == reaction_types.TRANSCRIPTIONAL_TRANSLATIONAL_ACTIVATION:
            reaction.modifier_role = participant_roles.STIMULATOR
        elif reaction.reaction_type == reaction_types.TRANSCRIPTIONAL_TRANSLATIONAL_REPRESSION:
            reaction.modifier_role = participant_roles.INHIBITOR
        elif reaction.reaction_type == reaction_types.CATALYSIS:
            reaction.modifier_role = participant_roles.CATALYST
        elif reaction.reaction_type == reaction_types.TRANSLOCATION:
            reaction.modifier_role = participant_roles.TRANSPORTER
        else:
            reaction.modifier_role = participant_roles.MODIFIER

        # product roles
        reaction.product_role = participant_roles.PRODUCT




class reaction_subtypes:
    """
    Class to hold constants for different types of kinetic laws.
    These constants are used to categorize reactions based on their kinetic behavior.
    """
    # Define constants for reaction subtypes
    BINDING_WITH_MODIFIER = 'binding-with-modifier'
    BINDING_WITHOUT_MODIFIER = 'binding-without-modifier'

    DISSOCIATION_WITH_MODIFIER = 'dissociation-with-modifier'
    DISSOCIATION_WITHOUT_MODIFIER = 'dissociation-without-modifier'

    CATALYSIS_WITH_SUBSTRATE_AND_WITH_MODIFIER = 'catalysis-with-substrate-and-with-modifier'
    CATALYSIS_WITHOUT_SUBSTRATE_AND_WITH_MODIFIER = 'catalysis-without-substrate-and-with-modifier'
    CATALYSIS_WITH_SUBSTRATE_AND_WITHOUT_MODIFIER = 'catalysis-with-substrate-and-without-modifier'
    CATALYSIS_WITHOUT_SUBSTRATE_AND_WITHOUT_MODIFIER = 'catalysis-without-substrate-and-without-modifier'

    DEGRADATION_WITH_MODIFIER = 'degradation-with-modifier'
    DEGRADATION_WITHOUT_MODIFIER = 'degradation-without-modifier'

    PROTEIN_DEACTIVATION_WITH_MODIFIER = 'protein-deactivation-with-modifier'

    PROTEIN_ACTIVATION_WITH_MODIFIER = 'protein-activation-with-modifier'

    TRANSCRIPTIONAL_ACTIVATION_WITH_MODIFIER = 'transcriptional-activation-with-modifier'
    TRANSLATIONAL_ACTIVATION_WITH_MODIFIER = 'translational-activation-with-modifier'
    TRANSCRIPTIONAL_OR_TRANSLATIONAL_ACTIVATION_WITH_MODIFIER = 'transcriptional-or-translational-activation-with-modifier'

    TRANSCRIPTIONAL_REPRESSION_WITH_MODIFIER = 'transcriptional-repression-with-modifier'
    TRANSLATIONAL_REPRESSION_WITH_MODIFIER = 'translational-repression-with-modifier'
    TRANSCRIPTIONAL_OR_TRANSLATIONAL_REPRESSION_WITH_MODIFIER = 'transcriptional-or-translational-repression-with-modifier'

    TRANSLOCATION_WITH_MODIFIER = 'translocation-with-modifier'
    TRANSLOCATION_WITHOUT_MODIFIER = 'translocation-without-modifier'

    UNKNOWN_ACTIVATION_WITH_MODIFIER = 'unknown-activation-with-modifier'
    UNKNOWN_ACTIVATION_WITHOUT_MODIFIER = 'unknown-activation-without-modifier'
    UNKNOWN_INHIBITION_WITH_MODIFIER = 'unknown-inhibition-with-modifier'
    UNKNOWN_INHIBITION_WITHOUT_MODIFIER = 'unknown-inhibition-without-modifier'

    UNKNOWN = 'unknown'

    @staticmethod
    def classify(reaction):
        """
        Return the subtype of the reaction based on its type, effect, mechanism,
        and the presence of substrates and modifiers (as currently in the reaction).
        """

        match reaction.reaction_type:
            case reaction_types.BINDING_OLIGOMERISATION:
                if reaction.has_modifiers():
                    return reaction_subtypes.BINDING_WITH_MODIFIER
                else:
                    return reaction_subtypes.BINDING_WITHOUT_MODIFIER

            case reaction_types.DISSOCIATION:
                if reaction.has_modifiers():
                    return reaction_subtypes.DISSOCIATION_WITH_MODIFIER
                else:
                    return reaction_subtypes.DISSOCIATION_WITHOUT_MODIFIER

            case reaction_types.CATALYSIS:
                if reaction.has_substrates() and reaction.has_modifiers():
                    return reaction_subtypes.CATALYSIS_WITH_SUBSTRATE_AND_WITH_MODIFIER
                elif reaction.has_modifiers():
                    return reaction_subtypes.CATALYSIS_WITHOUT_SUBSTRATE_AND_WITH_MODIFIER
                elif reaction.has_substrates():
                    return reaction_subtypes.CATALYSIS_WITH_SUBSTRATE_AND_WITHOUT_MODIFIER
                else:
                    return reaction_subtypes.CATALYSIS_WITHOUT_SUBSTRATE_AND_WITHOUT_MODIFIER

            case reaction_types.DEGRADATION_SECRETION:
                if reaction.has_modifiers():
                    return reaction_subtypes.DEGRADATION_WITH_MODIFIER
                else:
                    return reaction_subtypes.DEGRADATION_WITHOUT_MODIFIER

            case reaction_types.PROTEIN_DEACTIVATION:
                return reaction_subtypes.PROTEIN_DEACTIVATION_WITH_MODIFIER

            case reaction_types.PROTEIN_ACTIVATION:
                return reaction_subtypes.PROTEIN_ACTIVATION_WITH_MODIFIER

            case reaction_types.TRANSCRIPTIONAL_TRANSLATIONAL_ACTIVATION:
                if reaction.reaction_mechanism == 'transcription':
                    return reaction_subtypes.TRANSCRIPTIONAL_ACTIVATION_WITH_MODIFIER
                elif reaction.reaction_mechanism == 'translation':
                    return reaction_subtypes.TRANSLATIONAL_ACTIVATION_WITH_MODIFIER
                else:
                    return reaction_subtypes.TRANSCRIPTIONAL_OR_TRANSLATIONAL_ACTIVATION_WITH_MODIFIER

            case reaction_types.TRANSCRIPTIONAL_TRANSLATIONAL_REPRESSION:
                if reaction.reaction_mechanism == 'transcription':
                    return reaction_subtypes.TRANSCRIPTIONAL_REPRESSION_WITH_MODIFIER
                elif reaction.reaction_mechanism == 'translation':
                    return reaction_subtypes.TRANSLATIONAL_REPRESSION_WITH_MODIFIER
                else:
                    return reaction_subtypes.TRANSCRIPTIONAL_OR_TRANSLATIONAL_REPRESSION_WITH_MODIFIER

            case reaction_types.TRANSLOCATION:
                if reaction.has_modifiers():
                    return reaction_subtypes.TRANSLOCATION_WITH_MODIFIER
                else:
                    return reaction_subtypes.TRANSLOCATION_WITHOUT_MODIFIER

            case reaction_types.UNKNOWN:
                if reaction.reaction_effect == 'activation':
                    if reaction.has_modifiers():
                        return reaction_subtypes.UNKNOWN_ACTIVATION_WITH_MODIFIER
                    else:
                        return reaction_subtypes.UNKNOWN_ACTIVATION_WITHOUT_MODIFIER
                elif reaction.reaction_effect == 'inhibition':
                    if reaction.has_modifiers():
                        return reaction_subtypes.UNKNOWN_INHIBITION_WITH_MODIFIER
                    else:
                        return reaction_subtypes.UNKNOWN_INHIBITION_WITHOUT_MODIFIER
                else:
                   return reaction_subtypes.UNKNOWN

            case _:
                return reaction_subtypes.UNKNOWN


ALL_REACTION_SUBTYPES = [
    'binding-with-modifier',
    'binding-without-modifier',
    'dissociation-with-modifier',
    'dissociation-without-modifier',
    'catalysis-with-substrate-and-with-modifier',
    'catalysis-without-substrate-and-with-modifier',
    'catalysis-with-substrate-and-without-modifier',
    'catalysis-without-substrate-and-without-modifier',
    'degradation-with-modifier',
    'degradation-without-modifier',
    'protein-deactivation-with-modifier',
    'protein-activation-with-modifier',
    'transcriptional-activation-with-modifier',
    'translational-activation-with-modifier',
    'transcriptional-repression-with-modifier',
    'translational-repression-with-modifier',
    'translocation-with-modifier',
    'translocation-without-modifier',
    'unknown-activation-with-modifier',
    'unknown-activation-without-modifier',
    'unknown-inhibition-with-modifier',
    'unknown-inhibition-without-modifier',
    'unknown'
]