#!/usr/bin/env python2.7

"""
Search genome and find gene properties
"""

from __future__ import print_function
import sys
import yaml

from Bio import SeqIO

from mage_tool.IO import seqIO_to_genelist
from mage_tool.IO import create_config_tables
from mage_tool.translation import RBS_predict
from mage_tool.translation import dG_to_AU


if __name__ == '__main__':
    #Print HELP
    if len(sys.argv) <= 2 or "-h" in sys.argv or "--help" in sys.argv:
        print("usage:", sys.argv[0], "genome genome.config gene1 [gene2] [gene3] ..")
        exit()

    print("Parsing genome..")
    genome = SeqIO.read(sys.argv[1], "genbank")

    print("Loading config file..")
    with open(sys.argv[2]) as cfg:
        config = yaml.safe_load(cfg)
        config = create_config_tables(config)

    include_genes = sys.argv[3:]
    genes = seqIO_to_genelist(genome, config, include_genes)

    for gene in include_genes:
        if gene not in genes:
            print(gene, "not found..")

    for gene in genes.values():
        expr_lvl = RBS_predict(str(gene.leader), str(gene.cds))
        expr_lvl = dG_to_AU(expr_lvl)
        print(">>{}".format(gene))
        print("Strand:    ", gene.strand)
        print("Pos:       ", gene.pos)
        print("promoter:  ", gene.promoter)
        print("Leader:    ", gene.leader)
        print("Leader wb :", gene.leader_wobble)
        print("Expression: {:.2f}".format(expr_lvl))
        print("cds:")
        for i in range(0,len(gene.cds), 60):
            print(gene.cds[i:i+60])
        print("\n")
