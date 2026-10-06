'''
BoolNet: the Boolean model in the BoolNet format ("targets, factors", one rule per line), readable by BoolNet
(R), pyboolnet and BoolDog. The rules are the TabularQual model's (pss_export.boolean.TabularQual: the same
machinery, species ids and rules); species without a rule of their own (inputs) keep their value ("x, x"),
as BoolNet needs every variable defined. A comment above each rule lists the reactions it is built from.
Optionally a node file: per species id its PSS node, display label, form, location and type.
'''
from ..entity_classes import Node
from ..utils import write_tsv

NODE_COLUMNS = ['id', 'display_label', 'name', 'node_type', 'form', 'location']


def create_boolnet(tabqual, filename, nodes_file=None):
    ''' Write the BoolNet file (and the node file) of a TabularQual model (create_transitions() done); returns
    the number of rules (inputs included) '''
    transitions = {t.target: t for t in tabqual.transitions}
    targets = sorted(tabqual.species_dict)
    lines = ['targets, factors']
    for target in targets:
        if target in transitions:
            reactions = tabqual.rules_rx[target]
            lines += [f"# {', '.join(reactions)}", f'{target}, {transitions[target].rule}']
    inputs = [t for t in targets if t not in transitions]
    if inputs:
        lines.append('# inputs: no rule from PSS, they keep their value')
        lines += [f'{target}, {target}' for target in inputs]
    with open(filename, 'w', encoding='utf-8') as out:
        out.write('\n'.join(lines) + '\n')

    if nodes_file:
        nodes = tabqual.pss_adapter.nodes
        rows = []
        for (name, form, compartment), species_id in tabqual.species_ids.items():
            node = nodes.get(name) or Node(name)
            rows.append({'id': species_id, 'display_label': node.display_label, 'name': name,
                         'node_type': node.type, 'form': form, 'location': compartment})
        write_tsv(nodes_file, NODE_COLUMNS, sorted(rows, key=lambda r: r['id']))
    return len(targets)
