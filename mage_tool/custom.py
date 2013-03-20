#!/usr/bin/env python

"""
Custom mutations
"""

from __future__ import print_function

import random
import math
import re

from mutation_tools import find_mutation_box
from mutation_tools import compare_seqs
from oligo_design import Mutation
from helpers import degenerate_nucleotides


def custom_mutation(gene, mutation):
    """Do specified custom mutations"""
    #Find specified mutation. Specified as "upstream nucleotides"["original nucleotides" = "new nucleotides"]"downstream nucleotides".
    #Any part can be omitted, ie. new nucleotides can be omitted for a deletion.
    m = re.match("(\w*)\[(\w*)=(\w*)\](\w*)", mutation.upper())
    m = m.groups()
    before = m[0]+m[1]+m[3]
    after = m[0]+m[2]+m[3]
    #Find mutation offset on gene.
    offset = 0
    while str(gene.cds[offset:len(before)+offset]) != before:
        offset += 1
        if len(before)+offset > len(gene.cds):
            return False
    #Find mutationbox offset from supplied sequence.
    mutation = find_mutation_box(gene.cds[offset:len(before)+offset], after)
    #Add gene offset.
    mutation.pos += offset
    #Return mutation.after as uppercase?! WHY?!
    mutation.after = mutation.after.upper() 
    #Return gene mutation.
    return gene.do_mutation(mutation)
