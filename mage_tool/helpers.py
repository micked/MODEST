from string import maketrans

"""
DNA string tools
"""
    
def complement(x):
    """Complement a DNA string"""
    try:
        return x.complement()
    except AttributeError:
        return x.translate(maketrans("ATGCatgc", "TACGtacg"))
    
def reverse_complement(x):
    """Reverse complement"""
    try:
        return x.reverse_complement()
    except AttributeError:
        return x.translate(maketrans("ATGCatgc", "TACGtacg"))[::-1]

"""
Mutation tools
"""

def find_mutation_box(parent, mutation):
    """Find a mutation box based on parent and mutation

    parent and mutation must be same length
    """

    if len(parent) != len(mutation):
        raise ValueError("parent {} and mutation {} must be same length".format(len(parent), len(mutation)))

    #
    offset = 0
    for p,m in zip(parent,mutation):
        if p != m:
            break
        offset += 1

    length = len(parent)
    for p,m in zip(parent[::-1],mutation[::-1]):
        if p != m:
            break
        length -= 1

    mut = "{}={}".format(parent[offset:length], mutation[offset:length])
    return offset, mut


if __name__ == '__main__':
    print find_mutation_box("ATG", "AGC")