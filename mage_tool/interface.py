#!/usr/bin/env python

"""
General interface to <>
"""

from translation import replace_start_codon

def interface(adjustments, genes, genome):
    """General interface to <>

    adjustments is a parsed list
    genes is a dictionary with Gene objects
    genome is a Biopython Seq object or a string
    """

    mutations = list()

    for a in adjustments:
        op = operations[a["operation"]]
        gene = genes[a["gene"]]
        muts = op(gene, *a["options"])
        mutations.extend(muts)

    for m in mutations:
        print m

def start_codon_optimal(gene):
    mut = replace_start_codon(gene.cds, gene.pos, "ATG")
    if not mut:
        return []
    return [mut]

operations = {
    "start_codon_optimal": start_codon_optimal
    }