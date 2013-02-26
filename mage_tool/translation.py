#!/usr/bin/env python

from IO import Mutation
from helpers import find_mutation_box

def replace_start_codon(seq, pos, start_codon="ATG"):
    """Replace start codon"""

    if str(seq[0:3]).upper() == str(start_codon).upper():
        return False
    if len(start_codon) != 3:
        raise ValueError("start_codon must be 3 nucleotides long. Recieved {}".format(len(start_codon)))
    
    offset, mut = find_mutation_box(seq[0:3], start_codon)
    return Mutation("eq", mut, pos+offset)