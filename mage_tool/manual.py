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


def gene_mutation(gene, mutation):
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
    #Return gene mutation.
    return gene.do_mutation(mutation)


def genome_mutation(genome, mutation):
    """Validate a mutation (string in eq format) and return Mutation object."""
    pass


def residue_mutation(gene, mutations, codon_table, dgn_table):
    """Substitue residue.

    Substitutions are given as a list of amino acid mutations:

        X5A #Mutate X5 to A

    Deletions are denoted with a star:

        X5* #Delete X5

    Insertions are denoted with a star and a lower letter index:

        *5aA #Insert A between residue 5 and 6

    Stop-codons are denoted by $:

        X10$ #Mutate X10 to a STOP

    """
    #Find desired residue substitution. 
    m = re.match("([A-Z*$]+)(\d+)([A-Z*$]+)", mutation.upper())
    m = m.groups()
    
    old_AA = m[0]
    new_AA = m[2]
    pos = m[1]

