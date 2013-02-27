#!/usr/bin/env python

"""
Mutation generation routines related to translation initiation
"""


from IO import Mutation
from helpers import find_mutation_box

def replace_start_codon(gene, start_codon="ATG"):
    """Replace start codon"""

    if str(gene.cds[0:3]).upper() == str(start_codon).upper():
        return False
    if len(start_codon) != 3:
        raise ValueError("start_codon must be 3 nucleotides long. Recieved {}".format(len(start_codon)))
    
    offset, mut = find_mutation_box(gene.cds[0:3], start_codon)
    mutation = Mutation("eq", mut, offset)

    return gene.do_mutation(mutation)