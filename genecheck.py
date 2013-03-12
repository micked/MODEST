#!/usr/bin/env python

"""
Search genome and find gene properties
"""

from __future__ import print_function
import sys

from Bio import SeqIO

from mage_tool.IO import seqIO_to_genelist
from mage_tool.translation import RBSPredict


if __name__ == '__main__':
    #Print HELP
    if len(sys.argv) == 1 or "-h" in sys.argv:
        print("usage:", sys.argv[0], "genome gene1 [gene2] [gene3] ..")
        exit()

    print("Parsing genome..")
    genome = SeqIO.read(sys.argv[1], "genbank")

    include_genes = sys.argv[2:]
    genes = seqIO_to_genelist(genome, include_genes)

    for gene in include_genes:
        if gene not in genes:
            print(gene, "not found..")

    for gene in genes.values():
        expr_lvl = RBSPredict(str(gene.leader), str(gene.cds))["expr_lvl"]

        print(">>{}".format(gene))
        print("Strand:    ", gene.strand)
        print("Pos:       ", gene.pos)
        print("promoter:  ", gene.promoter)
        print("Leader:    ", gene.leader)
        print("Expression: {:.2f}".format(expr_lvl))
        print("cds:")
        for i in range(0,len(gene.cds), 60):
            print(gene.cds[i:i+60])
        print("\n")
