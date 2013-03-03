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