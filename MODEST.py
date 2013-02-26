#!/usr/bin/env python

"""
Commandline interface to <>

Yadda yadda yadda
"""

import argparse

from Bio import SeqIO

from mage_tool.interface import interface
from mage_tool.IO import seqIO_to_genelist


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("adjustments", help="Adjustment list")
    parser.add_argument("genome", help="Annotated genome")
    args = parser.parse_args()

    include_genes = set()
    print("Parsing adjustments..")
    adjustments = list()
    with open(args.adjustments) as f:
        for line in f:
            line = line.split()
            adjustments.append({"gene": line[0], "operation": line[1], "options": ()})
            include_genes.add(line[0])

    print("Parsing genome..")
    genome = SeqIO.read(args.genome, "genbank")
    print("Collecting gene list..")
    genes = seqIO_to_genelist(genome, include_genes)

    interface(adjustments, genes, genome.seq)