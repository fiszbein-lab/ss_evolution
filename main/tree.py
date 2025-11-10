import ete3

from itertools import combinations


def get_node_age(tree: ete3.Tree):
    node_age = dict()

    for node in tree.traverse('postorder'):
        node_name = node.name

        if node.is_leaf():
            node_age[node_name] = node.name
        else:
            _, age = node.get_farthest_leaf()
            node_age[node_name] = age

    return node_age


def fitch(tree: ete3.Tree):
    tree_cp = tree.copy('deepcopy')
    return _fitch(tree_cp)


def viz(tree: ete3.Tree, save_as: str) -> None:
    tree_style = ete3.TreeStyle()
    tree_style.show_leaf_name = False

    pos_node_style = ete3.NodeStyle()
    pos_node_style['fgcolor'] = "red"

    ori_node_style = ete3.NodeStyle()
    ori_node_style['fgcolor'] = "yellow"

    for node in tree.traverse('preorder'):
        if node.is_leaf():
            node.add_face(ete3.TextFace(node.name), column=0)
            node.add_face(ete3.TextFace(node.state), column=0)

        if node.state == {1}:
            node.set_style(pos_node_style)

        if node.gain:
            node.set_style(ori_node_style)

    # tree.show(tree_style=tree_style)
    tree.render(save_as, w=6, units="in", tree_style=tree_style)


def _fitch(tree: ete3.Tree):
    annot_trees = set()

    min_changes = _fitch_pass1(tree)
    ambiguous_nodes = _fitch_pass2(tree)

    if ambiguous_nodes:
        r = len({node for node, _ in ambiguous_nodes})

        if r > 10:
            # Ignore — accounts for a minority of cases and becomes too time-
            # consuming.
            return min_changes, annot_trees

        for tree_state in combinations(ambiguous_nodes, r=r):
            if not _unique_nodes(tree_state):
                continue

            node_state = {node.name: state for node, state in tree_state}

            changes = {'gain': set(), 'loss': set()}
            tot_changes = 0

            for node in tree.traverse('preorder'):
                node_name = node.name

                if node.is_leaf():
                    continue

                if node_name in node_state:
                    state = node_state[node_name]
                else:
                    state, = node.state

                if node.is_root() and state == 1:
                    # Here, the assignment of state 1 to the root is a de
                    # novo gain we've introduced, but it won't be counted
                    # below because there is no more ancestral 0 from which
                    # to track the change as a gain, so we increment the
                    # introduce change manually.
                    tot_changes += 1
                    changes['gain'].add(node.name)

                for desc in node.children:
                    desc_name = desc.name

                    if desc_name in node_state:
                        desc_state = node_state[desc_name]
                    else:
                        desc_state, = desc.state

                    match [state, desc_state]:
                        case [0, 1]:
                            changes['gain'].add(desc_name)
                            tot_changes += 1
                        case [1, 0]:
                            changes['loss'].add(desc_name)
                            tot_changes += 1
                        case _:
                            continue

            if tot_changes == min_changes:
                tree_cp = tree.copy('deepcopy')

                for node in tree_cp.traverse('postorder'):
                    name = node.name

                    if name in changes['gain']:
                        node.gain = True

                    if name in changes['loss']:
                        node.loss = True

                    # `node` was ambiguous and assigned a state; annotate the
                    #  tree accordingly.
                    if name in node_state:
                        state = node_state[name]
                        node.state = {state}

                annot_trees.add(tree_cp)
    else:
        tree_cp = tree.copy('deepcopy')

        for node in tree_cp.traverse(strategy='preorder'):
            if node.state == {0}:
                for desc in node.children:
                    if desc.state == {1}:
                        desc.gain = True

            if node.state == {1}:
                for desc in node.children:
                    if desc.state == {0}:
                        desc.loss = True

        annot_trees.add(tree_cp)

    pseudo_root = _get_pseudo_root(tree)
    if pseudo_root.state == {1}:
        # If the origin was assigned to the root, the true origin is unknown;
        # it could be more ancestral than the common ancestor.
        annot_trees = set()

    return min_changes, annot_trees


def _fitch_pass1(tree: ete3.Tree) -> int:
    min_changes = 0

    for node in tree.traverse(strategy='postorder'):
        node.gain = False
        node.loss = False

        if node.is_leaf():
            continue

        desc_states = [d.state for d in node.children]
        shared_state = set.intersection(*desc_states)

        if shared_state:
            operation = '∩'
            state = shared_state
        else:
            operation = '∪'
            state = set.union(*desc_states)
            min_changes += 1

        node.state = state
        node.operation = operation

    return min_changes


def _fitch_pass2(tree: ete3.Tree):
    ambiguous_nodes = list()

    for node in tree.traverse(strategy='preorder'):
        if node.is_leaf():
            continue

        if node.is_root():
            if len(node.state) > 1:
                ambiguous_nodes += [(node, state) for state in (0, 1)]
        else:
            mra, *_ = node.get_ancestors()

            if mra.state.issubset(node.state):
                updated_state = node.state & mra.state
            else:
                if node.operation == '∪':
                    updated_state = mra.state | node.state
                else:
                    desc_states = set.union(*[d.state for d in node.children])
                    updated_state = node.state | (mra.state & desc_states)

            node.state = updated_state
            if len(updated_state) > 1:
                ambiguous_nodes += [(node, state) for state in (0, 1)]

    return ambiguous_nodes


def _unique_nodes(node_states):
    return len({node for node, _ in node_states}) == len(node_states)


def _get_pseudo_root(tree: ete3.Tree):
    leaves = tree.get_leaf_names()
    return tree.get_common_ancestor(leaves)
