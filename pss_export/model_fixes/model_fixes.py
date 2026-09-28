'''
The models auto-generated from PSS have two major issues that leave parts of it disconnected. This occurs in the SBML, SBGN, and BoolNet models.

The issues are:
1) Genes and gene products occur in multiple forms, e.g.
      - "protein" form (e.g. unphosphorylated) as output from a translation reaction
      - "protein_active" form (e.g. phosphorylated) acting as a modifier in a downstream reaction.
   However, some are missing conversion reactions to connect the different forms
   (e.g. a protein activation reaction).
2) A number of nodes occur in multiple compartments as different species,
   but are missing transport reactions to connect them.


The following solutions are implemented:
1) We assume the unmodified protein form is the active form (i.e. no protein activation is required for the protein to be in the "active" form). To connect the model, we therefore change the output of translation reactions from "protein" to "protein_active", thus directly connecting them to downstream reactions where they act as the modifier. This is done by:
    a) Identify all nodes that occur as both "protein" and "protein_active" forms, and are not connected.
    b) For each of these nodes, identify all translation reactions that produce the "protein" form.
    c) Change the output of these translation reactions to the "protein_active" form.

2) For each species that occurs in multiple compartments, add transport reactions, from where the node is produced to where it is consumed. This is done by:
    a) Identify all nodes that occur in multiple compartments as different species, and are not connected.
    b) For each of these nodes, identify all compartments where it is produced (output from a reaction) and where it is consumed (input to a reaction).
    c) For each pair of compartments (from production to consumption), add a transport reaction if it does not already exist.

'''

from collections import defaultdict
from rich.console import Console
import matplotlib.pyplot as plt
import networkx as nx

from ..pss.pss_reaction_definitions import reaction_types, reaction_classes
from ..pss.config import pss_export_config

from ..entity_classes import IDTracker, Reaction, Species
from .graph import Graph, GraphVisualizer

console = Console()

class ModelFixTypes:
    """Class to hold different types of fixes and their explanations."""

    # Define fix types as class attributes
    FORM_PROTEIN_PRODUCT_ACTIVE = "fix-form:protein-product-active"
    FORM_COMPLEX_SUBSTRATE_ACTIVE = "fix-form:complex-substrate-active"
    FORM_COMPLEX_PRODUCT_ACTIVE = "fix-form:complex-product-active"

    LOC_TRANSPORT_TO_CONSUMED = "fix-loc:transport-to-consumed-compartment"
    LOC_TRANSPORT_FROM_PRODUCED = "fix-loc:transport-from-produced-compartment"
    LOC_TRANSPORT_PROTEIN_ER_TO_CYT = "fix-loc:protein-er-to-cyt"
    LOC_TRANSPORT_FROM_CYT = "fix-loc:from-cyt"

    # Explanations for each fix type
    _explanations = {
        FORM_PROTEIN_PRODUCT_ACTIVE: "Change the output of these transcription/translation reactions from 'protein' to 'protein_active'.",
        FORM_COMPLEX_SUBSTRATE_ACTIVE: "Change the substrate of these binding/oligomerisation reactions to 'protein_active' form. ",
        FORM_COMPLEX_PRODUCT_ACTIVE: "Change the output of these binding/oligomerisation reactions to 'complex_active' form. ",

        LOC_TRANSPORT_TO_CONSUMED: "Based on where node is consumed, add new transport reaction(s), transporting node from where it is produced to where it is consumed.",
        LOC_TRANSPORT_FROM_PRODUCED: "Based on where node is produced, add new transport reaction(s), transporting node from where it is produced to where it is consumed.",
        LOC_TRANSPORT_PROTEIN_ER_TO_CYT: "For proteins, if it exists in the 'endoplasmic reticulum' and is consumed in the 'cytoplasm', transport it from the 'endoplasmic reticulum' to the 'cytoplasm'.",
        LOC_TRANSPORT_FROM_CYT: "As last resort, if node exists in the 'cytoplasm' and is consumed elsewhere, transport if from the 'cytoplasm' to all consuming compartments, where it is not already present."
    }

    @classmethod
    def get_explanation(cls, fix_type):
        """Get the explanation for a given fix type."""
        return cls._explanations.get(fix_type, "Unknown fix type.") + f" ({fix_type})"


class ReactionFix:
    '''Hold the info to be able to apply a reaction fix later.'''

    def __init__(self,
                 reaction_id,
                 species_role,
                 name,
                 new_form=None,
                 new_location=None,
                 fix_type='model-fix'):
        '''
        Parameters
        ----------
        reaction_id: str
        species_role: 'substrate', 'product', or 'modifier'
        name: str
        new_form: str or None
        new_location: str or None
        fix_type: str
        '''

        if new_location and new_form:
            raise ValueError("Cannot specify both new_form and new_location.")

        self.reaction_id = reaction_id
        self.species_role = species_role
        self.name = name
        self.new_form = new_form
        self.new_location = new_location
        self.fix_type = fix_type

    def __repr__(self):
        if self.new_form:
            change_str = f"new_form='{self.new_form}'"
        elif self.new_location:
            change_str = f"new_location='{self.new_location}'"
        else:
            change_str = "no change specified"
        return f"""ReactionFix(reaction_id='{self.reaction_id}', species_role='{self.species_role}', name='{self.name.split("[")[0]}', {change_str})"""


class TransportReaction:
    '''Hold the information to create a new transport reaction.'''

    def __init__(self,
                 name,
                 form,
                 source_compartment,
                 target_compartment,
                 fix_type='model-fix'):
        '''
        Parameters
        ----------
        species_name: str
            The name of the node to be transported.
        form: str
            The form of the species (e.g., 'protein', 'protein_active').
        source_compartment: str
            The compartment where the species originates.
        target_compartment: str
            The compartment where the species is transported to.
        '''
        self.name = name
        self.form = form
        self.source_compartment = source_compartment
        self.target_compartment = target_compartment
        self.fix_type = fix_type

    def __repr__(self):
        return (
            f"""TransportReaction(name='{self.name.split("[")[0]}', form='{self.form}', """
            f"""source_compartment='{self.source_compartment}', target_compartment='{self.target_compartment}')"""
        )

class ModelFixer:

    def __init__(self, pss_adapter, interactive=False, apply_fixes=True):
        self.pss_adapter = pss_adapter
        self.interactive = interactive
        self.apply_fixes = apply_fixes
        self.console = Console()

    def identify_model_fixes(self, max_iterations=5):
        ''' Identify and optionally apply model fixes to the collected reactions.
        Parameters
        ----------
        pss_adapter: PSSAdapter
        interactive: bool
            If True, show plots and ask user to confirm fixes.
        apply_fixes: bool
            If True, apply the fixes to the pss_adapter (only relevant if interactive is False).
        '''
        self.console.rule("[bold red]Starting model fixing")

        # 1) Fix protein activation issues, ignoring location for now
        self.console.print(
            "Step 1: Fixing node 'form' issues, ignoring location for now.")

        # rerun the form fixing until no more fixes are found or we reach a max number of iterations
        self.console.print(
            f"Running form fixing until no more fixes are found (max {max_iterations} iterations)."
        )

        if self.interactive:
            plt.ion()

        num_form_fixes_applied = 0
        for iteration in range(max_iterations):
            if iteration > 0:
                console.print("Rerunning to check for more fixes...\n")
            self.console.print(f"Iteration {iteration + 1}:")
            num_form_fixes_applied_here = self._identify_model_fixes(
                part='form')
            if num_form_fixes_applied_here == 0:
                break

            num_form_fixes_applied += num_form_fixes_applied_here
            self.console.print(f"Applied {num_form_fixes_applied} fixes\n")

        if num_form_fixes_applied_here > 0:
            self.console.print("Max iterations reached, stopping.")

        self.console.print(
            f"Step 1 complete, applied {num_form_fixes_applied} fixes to species forms.\n"
        )

        # 2) Fix location issues
        self.console.print("Step 2: Fixing node 'location' issues.")

        # rerun the location fixing until no more fixes are found or we reach a max number of iterations
        self.console.print(
            f"Running location fixing until no more fixes are found (max {max_iterations} iterations)."
        )

        num_location_fixes_applied = 0
        for iteration in range(max_iterations):
            if iteration > 0:
                console.print("Rerunning to check for more fixes...\n")
            self.console.print(f"Iteration {iteration + 1}:")
            num_location_fixes_applied_here = self._identify_model_fixes(
                part='location')
            if num_location_fixes_applied_here == 0:
                break

            num_location_fixes_applied += num_location_fixes_applied_here
            self.console.print(f"Applied {num_location_fixes_applied} fixes\n")

        if num_location_fixes_applied_here > 0:
            self.console.print("Max iterations reached, stopping.")

        self.console.print(
            f"Step 2 complete, applied {num_location_fixes_applied} fixes to species locations.\n"
        )

        # Summary
        self.console.rule(
            f"[bold red]Model fixing complete, applied {num_form_fixes_applied + num_location_fixes_applied} fixes."
        )

        if self.interactive:
            plt.ioff()
            plt.close('all')

    def _identify_model_fixes(self, part='location'):

        num_fixes_applied = 0

        graph = Graph(self.pss_adapter, location=part == 'location')
        problematic_nodes = graph.find_problematic_nodes()
        self.console.print(
            f"Model fixing (considering {part}): Found {len(problematic_nodes)} problematic nodes."
        )

        if len(problematic_nodes) == 0:
            self.console.print("Nothing to fix.")
            return num_fixes_applied

        self.console.print()

        if part == 'location':
            fixing_func = self._suggest_fixes_location
        else:
            fixing_func = self._suggest_fixes_form

        for node_name, data in problematic_nodes.items():
            self.console.rule(f"[bold red]Node {node_name.split('[')[0]}")
            species = data['species']
            subgraph = data['subgraph']
            fixes = fixing_func(node_name, species, subgraph)
            if len(fixes) == 0:
                console.print("\nNo automatic fix identified for this node.")
            else:
                console.print("\nSuggested fixes:")
                # print the suggested fixes
                for fix in fixes:
                    console.print(f" - {fix}")
                console.print()

            if self.interactive:
                num_fixes_applied_here = self._interactive_fixing(
                    node_name, species, subgraph, fixes)

                # if returns -1, user chose to quit
                if num_fixes_applied_here == -1:
                    break
                num_fixes_applied += num_fixes_applied_here

            else:
                if fixes and self.apply_fixes:
                    num_fixes_applied += self.apply_model_fixes(fixes)
                else:
                    self.console.print("No fixes applied.")

        self.console.rule()
        return num_fixes_applied

    def _suggest_fixes_form(self, node_name, species, subgraph):

        # suggest fixes by case...
        fixes = []
        fix_identified = False

        forms = {
            s: subgraph.nodes[s].get('form')
            for s in species if subgraph.nodes[s].get('form')
        }
        self.console.print(f"All node forms: {forms}\n")

        # what are the separate components of the subgraph?
        components = list(nx.weakly_connected_components(subgraph))
        self.console.print(f"{len(components)} components found:")

        # which forms are in each component?
        comp_forms = []
        for i, comp in enumerate(components):
            forms_in_comp = {
                subgraph.nodes[n].get('form')
                for n in comp if n in species
            }
            self.console.print(f" - Component {i} has forms: {forms_in_comp}")
            comp_forms.append(forms_in_comp)
        console.print()

        # if there are more than two components, too complex to fix automatically
        if len(components) > 2:
            console.print(
                "More than two components, too complex to fix automatically.")

        # Multiple components, some with "protein" and some with "protein_active"
        # and the components with "protein" has one or more translation reactions,
        # suggest changing the output of those translation reactions to "protein_active"
        elif ( (protein_alone := {'protein'} in comp_forms) or ({'gene', 'protein'} in comp_forms) ) and \
                       ( ({'protein_active'} in comp_forms) or ({'gene', 'protein_active'} in comp_forms) ):
            # find the component with "protein"
            if protein_alone:
                protein_comp = components[comp_forms.index({'protein'})]
            else:
                protein_comp = components[comp_forms.index({'gene', 'protein'})]

            # option 1 - we can change translation/transcription to output protein active
            # find all translation reactions in this component
            console.print(
                "Looking for translation reactions in the 'protein' component..."
            )
            translation_reactions = [
                r for r in protein_comp if subgraph.nodes[r].get('type') == 'reaction' \
                # reaction is translation/transcription
                and subgraph.nodes[r].get('reaction_type') in reaction_classes.TRANSCRIPTIONAL_TRANSLATIONAL \
                # node (node_name) has a 'product' edge from the reaction
                and any(data.get('type') == 'product' and v in species for u, v, data in subgraph.out_edges(r, data=True))
            ]
            if translation_reactions:
                fix_type = ModelFixTypes.FORM_PROTEIN_PRODUCT_ACTIVE
                console.print(
                    f"[green]Suggestion: {ModelFixTypes.get_explanation(fix_type)}"
                )
                for reaction_id in translation_reactions:
                    console.print(
                        f" - transcription/translation reaction {reaction_id} producing 'protein' form."
                    )
                    fix = ReactionFix(reaction_id=reaction_id,
                                      species_role='product',
                                      name=node_name,
                                      new_form='protein_active',
                                      fix_type=fix_type)
                    fixes.append(fix)
                fix_identified = True
            else:
                console.print(
                    "[orange]No translation reactions found in the 'protein' component."
                )

            if not fix_identified:

                # option 2 - we can change binding/oligomerisation to use protein_active as substrate
                # find all binding/oligomerisation reactions in these components
                console.print(
                    "Looking for binding/oligomerisation reactions in the 'protein' component..."
                )
                binding_reactions = [
                    r for r in protein_comp if subgraph.nodes[r].get('type') == 'reaction' \
                    # reaction is binding/oligomerisation
                    and subgraph.nodes[r].get('reaction_type') == reaction_types.BINDING_OLIGOMERISATION \
                    # node (node_name) has a 'substrate' edge to the reaction
                    and any(data.get('type') == 'substrate' and u in species for u, v, data in subgraph.in_edges(r, data=True))
                ]
                if binding_reactions:
                    fix_type = ModelFixTypes.FORM_COMPLEX_SUBSTRATE_ACTIVE
                    console.print(
                        f"[green]Suggestion: {ModelFixTypes.get_explanation(fix_type)}"
                    )
                    for reaction_id in binding_reactions:
                        console.print(
                            f" - binding/oligomerisation reaction {reaction_id} using 'protein' form as substrate."
                        )
                        fix = ReactionFix(reaction_id=reaction_id,
                                          species_role='substrate',
                                          name=node_name,
                                          new_form='protein_active',
                                          fix_type=fix_type)
                        fixes.append(fix)
                    fix_identified = True
                else:
                    console.print(
                        "[dark_orange]No binding/oligomerisation reactions found in the 'protein' component, cannot suggest fix."
                    )

        # similarly if we have complex and complex_active (with no link between), then suggest to make all 'complex' to 'complex_active'
        elif {'complex'} in comp_forms and {'complex_active'} in comp_forms:
            console.print(
                "Found components with expected forms 'complex' and 'complex_active'."
            )
            # find the components with "complex"
            complex_comps = [
                components[i] for i, forms in enumerate(comp_forms)
                if forms == {'complex'}
            ]
            for complex_comp in complex_comps:
                # find binding/oligomerisation reactions in this component
                binding_reactions = [
                    n for n in complex_comp
                    if subgraph.nodes[n].get('type') == 'reaction'
                    and subgraph.nodes[n].get('reaction_type') ==
                    reaction_types.BINDING_OLIGOMERISATION
                ]
                if binding_reactions:
                    fix_type = ModelFixTypes.FORM_COMPLEX_PRODUCT_ACTIVE
                    console.print(
                        f"[green]Suggestion: {ModelFixTypes.get_explanation(fix_type)}"
                    )
                    for reaction_id in binding_reactions:
                        console.print(
                            f" - binding/oligomerisation reaction {reaction_id} producing 'complex' form."
                        )
                        fix = ReactionFix(reaction_id=reaction_id,
                                          species_role='product',
                                          name=node_name,
                                          new_form='complex_active',
                                          fix_type=fix_type)
                        fixes.append(fix)
                    fix_identified = True

            if not fix_identified:
                console.print(
                    "[dark_orange]No binding/oligomerisation reactions found in the 'complex' component, cannot suggest fix."
                )

        return fixes

    def _suggest_fixes_location(self, node_name, species, subgraph):

        # suggest fixes by case...
        fixes = []
        fix_identified = False

        # in the case of fixing the locations, we have to treat each form of the node separately...
        # collect all the forms of this node {dict form --> list of species ids}
        forms_dict = defaultdict(list)
        for s in species:
            form = subgraph.nodes[s].get('form')
            if form:
                forms_dict[form].append(s)

        # what are the separate components of the subgraph?
        components = list(nx.weakly_connected_components(subgraph))
        self.console.print(f"{len(components)} components found.")

        for form, form_species in forms_dict.items():

            console.print(f"[bold blue]\nConsidering form '{form}'")

            locations = {
                s: subgraph.nodes[s].get('compartment')
                for s in form_species if subgraph.nodes[s].get('compartment')
            }

            if len(locations) <= 1:
                self.console.print(
                    f"[dark_orange]Form '{form}' only occurs in one location, no fix suggested."
                )
                continue

            self.console.print(f"All node ('{form}') locations: {locations}")

            # which locations are in each component?
            comp_locations = []
            for i, comp in enumerate(components):
                locations_in_comp = {
                    subgraph.nodes[n].get('compartment')
                    for n in comp if n in form_species
                }
                if locations_in_comp:
                    self.console.print(
                        f" - Component {i} has locations: {locations_in_comp}")
                    comp_locations.append(locations_in_comp)
            console.print()

            # find all compartments where the form_species is produced (product of a reaction)
            # (excluding produced by transport reactions)
            producing_compartments = {
                subgraph.nodes[v].get('compartment')
                for r in subgraph.nodes
                if subgraph.nodes[r].get('type') == 'reaction' and subgraph.
                nodes[r].get('reaction_type') != reaction_types.TRANSLOCATION
                for u, v, data in subgraph.out_edges(r, data=True)
                if data.get('type') == 'product' and v in form_species
            }

            # find all places the form_species is produced by transport reactions
            transport_producing_compartments = {
                subgraph.nodes[v].get('compartment')
                for r in subgraph.nodes
                if subgraph.nodes[r].get('type') == 'reaction' and subgraph.
                nodes[r].get('reaction_type') == reaction_types.TRANSLOCATION
                for u, v, data in subgraph.out_edges(r, data=True)
                if data.get('type') == 'product' and v in form_species
            }

            # find all compartments where the form_species is consumed (substrate of a reaction)
            # (including consumed by transport reactions)
            consuming_compartments = {
                subgraph.nodes[u].get('compartment')
                for r in subgraph.nodes
                if subgraph.nodes[r].get('type') == 'reaction'
                for u, v, data in subgraph.in_edges(r, data=True)
                if data.get('type') in ['substrate', 'modifier']
                and u in form_species
            }

            console.print(
                f"Producing compartments: {','.join(list(producing_compartments)) if producing_compartments else 'None'}"
            )
            console.print(
                f"Produced by transport into compartments: {','.join(list(transport_producing_compartments)) if transport_producing_compartments else 'None'}"
            )
            console.print(
                f"Consuming compartments: {','.join(list(consuming_compartments)) if consuming_compartments else 'None'}"
            )
            console.print()

            # check if there are compartments where it is consumed but not
            # produced or transported to. Then we need to figure out how to get
            # it to these compartments from where it is produced (but at the
            # moment don't transport it from where it is transported to).
            need_to_transport_to = []
            for comp in consuming_compartments:
                if (comp in producing_compartments) or (
                        comp in transport_producing_compartments):
                    continue
                console.print(
                    f"Node is consumed in '{comp}' but not produced or transported to there."
                )
                need_to_transport_to.append(comp)

            if need_to_transport_to:
                console.print()
                # create a transport reaction to transport to these comps
                # first figure out where it comes from...
                if len(producing_compartments) == 0:
                    console.print(
                        "[dark_orange]Node is not produced in any compartment, no obvious compartment to transport from.\n"
                    )
                elif len(producing_compartments) > 1:
                    console.print(
                        "[dark_orange]Node is produced in multiple compartments, no obvious compartment to transport from.\n"
                    )
                else:
                    source_compartment = producing_compartments.pop()

                    fix_type = ModelFixTypes.LOC_TRANSPORT_TO_CONSUMED
                    console.print(
                        f"[green]Suggestion: {ModelFixTypes.get_explanation(fix_type)}"
                    )
                    for comp in need_to_transport_to:
                        console.print(
                            f" - transport '{form}' from '{source_compartment}' to '{comp}'"
                        )
                        fix = TransportReaction(name=node_name,
                                                form=form,
                                                source_compartment=source_compartment,
                                                target_compartment=comp,
                                                fix_type=fix_type)
                        fixes.append(fix)
                        fix_identified = True

            # check if there are compartments where it is produced but not
            # consumed. Then we need to figure out how to transport it from
            # these compartments to where it is consumed.
            need_to_transport_from = []
            for comp in producing_compartments:
                if comp in consuming_compartments:
                    continue
                console.print(
                    f"Node is produced in '{comp}' but not consumed there.")
                need_to_transport_from.append(comp)

            if need_to_transport_from:
                console.print()
                # create a transport reaction to transport from these comps
                # first figure out where it goes to...
                if len(consuming_compartments) == 0:
                    console.print(
                        "[dark_orange]Node is not consumed in any compartment, no obvious compartments to transport to.\n"
                    )
                else:
                    target_compartments = [c for c in consuming_compartments
                        # don't bother moving it to places it is already being produced
                        if not ((c in producing_compartments) or (c in transport_producing_compartments))]
                    if target_compartments:
                        fix_type = ModelFixTypes.LOC_TRANSPORT_FROM_PRODUCED
                        console.print(
                            f"[green]Suggestion: {ModelFixTypes.get_explanation(fix_type)}"
                        )
                        for target_compartment in target_compartments:
                            for comp in need_to_transport_from:
                                console.print(
                                    f" - transport '{form}' from '{comp}' to '{target_compartment}'"
                                )
                                fix = TransportReaction(name=node_name,
                                                        form=form,
                                                        source_compartment=comp,
                                                        target_compartment=target_compartment,
                                                        fix_type=fix_type)
                                fixes.append(fix)
                                fix_identified = True

            if not fix_identified:
                # if it's a protein or protein_active, and it exists in the
                # endoplasmic reticulum suggest transporting it from the er to the cytoplasm
                if form in [
                        'protein', 'protein_active'
                ] and 'endoplasmic reticulum' in locations.values():
                    if 'cytoplasm' in consuming_compartments:
                        fix_type = ModelFixTypes.LOC_TRANSPORT_PROTEIN_ER_TO_CYT
                        console.print(
                            f"[green]Suggestion: {ModelFixTypes.get_explanation(fix_type)}"
                        )
                        console.print(
                            f" - transport '{form}' from 'endoplasmic reticulum' to 'cytoplasm'."
                        )
                        fix = TransportReaction(name=node_name,
                                                form=form,
                                                source_compartment='endoplasmic reticulum',
                                                target_compartment='cytoplasm',
                                                fix_type=fix_type)
                        fixes.append(fix)
                        fix_identified = True

                # if it exists in the cytoplasm and is consumed in other compartments,
                # suggest transporting from cytoplasm to those compartments
                # if it's not already being transported there
                elif 'cytoplasm' in locations.values():
                    target_compartments = [c for c in consuming_compartments
                        # don't bother moving it to places it is already being produced
                        if not ((c == 'cytoplasm') or (c in transport_producing_compartments))]
                    if target_compartments:
                        fix_type = ModelFixTypes.LOC_TRANSPORT_FROM_CYT
                        console.print(
                            f"[green]Suggestion: {ModelFixTypes.get_explanation(fix_type)}"
                        )
                        for comp in target_compartments:
                            console.print(
                                f" - transport '{form}' from 'cytoplasm' to '{comp}'."
                            )
                            fix = TransportReaction(name=node_name,
                                                    form=form,
                                                    source_compartment='cytoplasm',
                                                    target_compartment=comp,
                                                    fix_type=fix_type)
                            fixes.append(fix)
                        fix_identified = True
        return fixes

    def _interactive_fixing(self, node_name, species, subgraph, fixes):
        # plot the subgraph
        GraphVisualizer.plot_subgraph(subgraph, species)

        num_applied = 0

        # if no fixes, wait for the user to skip to the next one or quit
        if not fixes:
            action = console.input(
                "Enter 's' to skip to the next node, or 'q' to quit: ").strip(
                ).lower()
            plt.clf()  # close plot once user has made a choice
            if action == 's':
                console.print(f"   - Skipping node {node_name.split('[')[0]}")
            elif action == 'q':
                console.print("Exiting model fixing.")
                return -1
            else:
                console.print("Invalid input, skipping.")

        # ask the user what to do if we found a fix
        else:
            action = console.input(
                "Enter 'a' to apply the fix(es), 's' to skip, or 'q' to quit: "
            ).strip().lower()
            plt.clf()  # close plot once user has made a choice
            if action == 'a':
                num_applied = self.apply_model_fixes(fixes)
            elif action == 's':
                console.print(f"   - Skipping node {node_name}")
            elif action == 'q':
                console.print("Exiting model fixing.")
                return -1
            else:
                console.print("Invalid input, skipping.")

        return num_applied

    def apply_model_fixes(self, fixes):
        num_fixes_applied = 0
        for fix in fixes:
            if isinstance(fix, ReactionFix):
                reaction = self.pss_adapter.reactions.get(fix.reaction_id)
                if not reaction:
                    self.console.print(
                        f"Warning: Reaction {fix.reaction_id} not found, cannot apply fix."
                    )
                    continue
                species_list = getattr(reaction, fix.species_role + 's', [])
                for species in species_list:
                    if species.name == fix.name:
                        if fix.new_form:
                            self.console.print(
                                f"Applying fix: Changing form of {species.name} in reaction {reaction.id} from {species.form} to {fix.new_form}."
                            )
                            species.form = fix.new_form
                        elif fix.new_location:
                            self.console.print(
                                f"Applying fix: Changing location of {species.name} in reaction {reaction.id} from {species.compartment} to {fix.new_location}."
                            )
                            species.compartment = fix.new_location
                        num_fixes_applied += 1
                reaction.export_notes += fix.fix_type

            elif isinstance(fix, TransportReaction):
                # Create a new transport reaction

                node_symbol = IDTracker.remove_nonalphanum(IDTracker.get_display_label(fix.name))
                reaction_id = f"transport_{node_symbol}_{pss_export_config.node_form_to_short[fix.form]}_{pss_export_config.compartment_to_short[fix.source_compartment]}_to_{pss_export_config.compartment_to_short[fix.target_compartment]}"
                if reaction_id in self.pss_adapter.reactions:
                    self.console.print(
                        f"Warning: Transport reaction {reaction_id} already exists, skipping."
                    )
                    continue
                self.console.print(
                    f"Applying fix: Adding transport reaction for {fix.name} from {fix.source_compartment} to {fix.target_compartment}."
                )
                new_reaction = Reaction(
                    reaction_id=reaction_id,
                    reaction_type=reaction_types.TRANSLOCATION,
                    export_notes=fix.fix_type,
                    reaction_properties={},
                    include_genes=self.pss_adapter.include_genes
                )
                new_reaction.add_substrate(Species(fix.name, fix.form, fix.source_compartment))
                new_reaction.add_product(Species(fix.name, fix.form, fix.target_compartment))

                self.pss_adapter.reactions[reaction_id] = new_reaction
                self.pss_adapter.additional_reactions.append(reaction_id)
                num_fixes_applied += 1

        self.console.print("[bold green]Model fixes applied!")
        return num_fixes_applied
