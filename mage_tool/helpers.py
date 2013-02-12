from string import maketrans

"""
DNA string tools
"""

def reverse(x):
    """Reverse string"""
    return x[::-1]
    
def complement(x):
    """Complement a DNA string"""
    return x.translate(maketrans("ATGCatgc", "TACGtacg"))
    
def reversecomplement(x):
    """Reverse complement"""
    return x.translate(maketrans("ATGCatgc", "TACGtacg"))[::-1]
