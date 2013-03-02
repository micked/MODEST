#!/usr/bin/env python

"""
Mutation generation routines related to translation initiation
"""


from oligo_design import Mutation
from DNA_tools import find_mutation_box
from DNA_tools import compare_seqs

def replace_start_codon(gene, start_codon="ATG"):
    """Replace start codon"""

    if str(gene.cds[0:3]).upper() == str(start_codon).upper():
        return False
    if len(start_codon) != 3:
        raise ValueError("start_codon must be 3 nucleotides long. Recieved {} ({})".format(start_codon, len(start_codon)))
    
    offset, mut = find_mutation_box(gene.cds[0:3], start_codon)
    mutation = Mutation("eq", mut, offset)

    return gene.do_mutation(mutation)


def translational_KO(gene, stop_codons=["TAG", "TAA", "TGA"], KO_frame=10):
    """Knock out a gene with a stop-codon with the least possible mutations

    KO_frame is how many codons the method will look for stop codon mutations.

    This method will prioritise single mutations, double mutations located
    beside each other, double mutations with a match inbetween, and finally
    triple mutations.
    """
    KO_frame = min(KO_frame, (len(gene.cds)-3)/3)

    KO = gene.cds[3:KO_frame*3+3]
    #m = target mutations, g = target groups
    for m,g in [(1,1), (2,1), (2,2), (3,1)]:
        #Loop codons
        for i in range(0, len(KO), 3):
            #Parent codon
            parent = KO[i:i+3]
            for child in stop_codons:
                #test mutitions and test groups
                tm,tg = compare_seqs(parent, child)
                if tm <= m and tg <= g:
                    #Do mutation
                    offset, mut = find_mutation_box(parent, child)
                    mutation = Mutation("eq", mut, offset+i+3)
                    return gene.do_mutation(mutation)
