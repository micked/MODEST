#!/usr/bin/env python2.7

"""
Search genome and find gene properties
"""

from __future__ import print_function
import os
import sys
import yaml
import argparse

from Bio import SeqIO

from mage_tool.IO import seqIO_to_genelist
from mage_tool.IO import create_config_tables
from mage_tool.operations.translation import RBS_predict
from mage_tool.operations.translation import dG_to_AU


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("config", help="Genome configuration", type=argparse.FileType("r"))
    parser.add_argument("--genome", help="Annotated genome. Leave empty to locate automatically.", default=None)
    parser.add_argument("genes", help="A number of genes to probe", nargs="*")
    args = parser.parse_args()

    print("Loading config file..")
    config = yaml.safe_load(args.config)
    cfg_basedir = os.path.abspath(os.path.dirname(args.config.name))
    config = create_config_tables(config, cfg_basedir)

    if not args.genome:
        cfg_basename = os.path.splitext(args.config.name)[0]
        genome_locations = [cfg_basename + ".gb",
                            cfg_basename + ".genbank",
                            config["Locus"],
                            config["Locus"] + ".gb",
                            config["Locus"] + ".genbank",
                            os.path.join(cfg_basedir, config["Locus"]),
                            os.path.join(cfg_basedir, config["Locus"] + ".gb"),
                            os.path.join(cfg_basedir, config["Locus"] + ".genbank")]

        for loc in genome_locations:
            if os.path.isfile(loc):
                genome_loc = loc
                break
        else:
            print("Genome not found, searched the following locations:")
            for loc in genome_locations:
                print(" *", loc)
            print("Use --genome to manually supply a genome")
            exit(1)
    else:
        genome_loc = args.genome

    print("Parsing genome..")
    genome = SeqIO.read(genome_loc, "genbank")

    include_genes = args.genes
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
        print("Leader:    ", gene.leader)
        print("Leader wb :", gene.leader.get_wobble_str())
        print("Expression: {:.2f}".format(expr_lvl))
        print()
        print("In operon: ", gene.in_operon)
        print("promoter:")
        for i in range(0,len(gene.promoter), 60):
            print(gene.promoter[i:i+60])
        print("@:", gene.promoter_pos)
        print()
        print("cds:")
        for i in range(0,len(gene.cds), 60):
            print(gene.cds[i:i+60])
        print("\n")
