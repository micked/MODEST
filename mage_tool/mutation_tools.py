"""
General Mutation tools
"""

from oligo_design import Mutation


def find_mutation_box(parent, child):
    """Find a mutation box based on parent and mutation

    parent and mutation must be same length
    """
    if len(parent) != len(child):
        raise ValueError("parent {} and child {} must be same length".format(len(parent), len(child)))

    parent = parent.upper()
    child = child.upper()

    offset = 0
    for p,m in zip(parent,child):
        if p != m:
            break
        offset += 1

    length = len(parent)
    for p,m in zip(parent[::-1],child[::-1]):
        if p != m:
            break
        length -= 1

    mut = "{}={}".format(parent[offset:length], child[offset:length])
    mut = Mutation("eq", mut, offset)
    return mut

def compare_seqs(parent, child):
    """Compare two sequences and return number of mutations

    return value is a tuple containing number of differences and a bool
    describing number of groups.
    """
    if len(parent) != len(child):
        raise ValueError("parent {} and child {} must be same length".format(len(parent), len(child)))

    parent = parent.upper()
    child = child.upper()

    muts = 0
    groups = 0
    on_mut = False

    for p,m in zip(parent,child):
        if p == m:
            on_mut = False
        else:
            muts += 1
            if not on_mut:
                groups += 1
            on_mut = True

    return muts, groups