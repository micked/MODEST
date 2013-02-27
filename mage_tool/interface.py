#!/usr/bin/env python

"""
General interface to <>
"""

import logging

import oligo_design
from translation import replace_start_codon

#Define a log
log = logging.getLogger("MODEST")
log.addHandler(logging.NullHandler())

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
        for mut in muts:
            oligo = oligo_design.mut_to_oligo(mut, genome, 90)
            mutations.append({"mutation": mut, "gene": gene, "oligo": oligo})

    for m in mutations:
        print m["mutation"], m["gene"].cds[0:10]

def start_codon_optimal(gene):
    mut = replace_start_codon(gene, "ATG")
    if not mut:
        log.debug("start_codon_optimal({}): Not mutating, optimal start codon found.".format(gene))
        return []
    log.info("start_codon_optimal({}): {}".format(gene, str(mut)))
    return [mut]

operations = {
    "start_codon_optimal": start_codon_optimal
    }