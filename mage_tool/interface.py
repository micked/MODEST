#!/usr/bin/env python

"""
General interface to <>
"""

import logging

from oligo_design import Oligo
from translation import replace_start_codon

#Define a log
log = logging.getLogger("MODEST")
log.addHandler(logging.NullHandler())

def interface(adjustments, genes, genome, config, project=None):
    """General interface to <>

    adjustments is a parsed list
    genes is a dictionary with Gene objects
    genome is a Biopython Seq object or a string
    """

    oligos = list()

    i = 0
    for a in adjustments:
        op = operations[a["operation"]]
        gene = genes[a["gene"]]
        muts = op(gene, *a["options"])
        for mut in muts:
            oligo = Oligo(mut, gene, project, i, oligo_len=90)
            oligo.make_oligo(genome)
            oligo.target_lagging_strand(config["replication"]["ori"], config["replication"]["ter"])
            oligos.append(oligo)
            i += 1

    return oligos

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