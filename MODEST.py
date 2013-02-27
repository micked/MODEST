#!/usr/bin/env python

"""
Commandline interface to <>

Yadda yadda yadda
"""

import logging
import argparse

from Bio import SeqIO

from mage_tool.interface import interface
from mage_tool.IO import seqIO_to_genelist


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("adjustments", help="Adjustment list")
    parser.add_argument("genome", help="Annotated genome")
    parser.add_argument("--log", help="Logfile, default MODEST.log. Use STDOUT to log all messages to screen, use - to disable", default="MODEST.log")
    args = parser.parse_args()

    #Set up logging
    format = "%(asctime)s %(name)-12s: %(levelname)-8s %(message)s"
    if args.log == "-":
        logging.basicConfig(level=logging.ERROR)
        console = logging.StreamHandler()
        console.setLevel(logging.ERROR)
        formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
        logging.getLogger('').addHandler(console)
    elif args.log.lower() == "stdout":
        logging.basicConfig(level=logging.DEBUG, format=format, datefmt='%Y-%m-%d %H:%M')
    else:
        logging.basicConfig(level=logging.DEBUG, format=format, datefmt='%Y-%m-%d %H:%M', filename='MODEST.log', filemode='w')
        console = logging.StreamHandler()
        console.setLevel(logging.WARNING)
        formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

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