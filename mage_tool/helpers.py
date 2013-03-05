from string import maketrans
import re

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

def valid_na(seq):
    """Check for valid DNA/RNA"""
    return bool(re.match(r"^[ATGCU]+$", seq, re.IGNORECASE))

def valid_rna(seq):
    """Check for valid RNA"""
    return bool(re.match(r"^[AGCU]+$", seq, re.IGNORECASE))

def valid_dna(seq):
    """Check for valid DNA"""
    return bool(re.match(r"^[ATGC]+$", seq, re.IGNORECASE))