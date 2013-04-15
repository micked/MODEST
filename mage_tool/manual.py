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
from helpers import dgn_to_nts


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
    mutation = find_mutation_box(before, after)
    #Add gene offset.
    mutation.pos += offset
    #Return gene mutation.
    return gene.do_mutation(mutation)


def genome_mutation(genome, mutation):
    """Validate a mutation (string in eq format) and return Mutation object."""
    pass


def residue_mutation(gene, mutations, codon_table, dgn_table, usage_table):
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
    new_dna_string = str(gene.cds)
    
    mutations = sorted(mutations, key=lambda x: x[1:-1])

    
    for mut in mutations:     
        #Find desired residue substitution. 
        m = re.match("([A-Z*$])(\d+)([a-z]?)([A-Z*$])", mut)
        m = m.groups()
        old_AA = m[0]
        new_AA = m[3]
        pos = m[1]
        pos_letter = m[2]
        
        dna_pos = int(pos)*3-3
        if old_AA != "*":
            dna_AA = str(gene.cds[dna_pos:dna_pos+3])
        else:
            dna_AA = "NNN"
            new_dna_string = "{}{}{}".format(new_dna_string[:dna_pos],dna_AA,new_dna_string[dna_pos:])
            
        new_codon = ""
        
        if new_AA != "*":
            
            for dnt, nt in zip(dgn_table[new_AA], dna_AA):
                if nt in dgn_to_nts[dnt]:
                    new_codon = "{}{}".format(new_codon, nt)
                elif len(dgn_to_nts[dnt]) == 1:
                    new_codon = "{}{}".format(new_codon, list(dgn_to_nts[dnt])[0])
                else:
                    new_codon = "{}{}".format(new_codon, dnt)
        
            pos_cds = [new_codon]
            for i, dnt in enumerate(new_codon):
                if len(dgn_to_nts[dnt]) > 1:
                    tmplist = list()
                    for dgn in dgn_to_nts[dnt]:
                        cdn = "{}{}{}".format(new_codon[:i], dgn, new_codon[i+1:])
                        tmplist.append(cdn)
                    pos_cds = tmplist
                    
            usage = -9
            for codon in pos_cds:
                if usage_table[codon][1] > usage:
                    new_codon = codon
                    usage = usage_table[codon][1]
        
        new_dna_string = "{}{}{}".format(new_dna_string[:dna_pos],new_codon,new_dna_string[dna_pos+3:])
    
    mutation = find_mutation_box(gene.cds, new_dna_string)
    return gene.do_mutation(mutation)
    pos = m[1]

