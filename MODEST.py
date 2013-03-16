#!/usr/bin/env python

"""
Commandline interface to <>

Yadda yadda yadda
"""

import logging
import argparse

import yaml
from Bio import SeqIO

from mage_tool.interface import interface
from mage_tool.IO import seqIO_to_genelist
from mage_tool.IO import oligolist_to_tabfile
from mage_tool.IO import parse_barcode_library
from mage_tool.IO import create_config_tables


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("adjustments", help="Adjustment list")
    parser.add_argument("genome", help="Annotated genome")   
    parser.add_argument("barcodes", help="Barcode library")    
    parser.add_argument("config", help="Genome configuration")
    parser.add_argument("--log", help="Logfile, default MODEST.log. Use STDOUT to log all messages to screen, use - to disable", default="MODEST.log")
    parser.add_argument("-p", "--project", help="Project name", default="Untitled")
    parser.add_argument("-o", "--output", help="Output file. Default <project>.out", default=False)
    args = parser.parse_args()

    #Set up logging
    format = "%(asctime)s %(name)-12s: %(levelname)-8s %(message)s"
    if args.log == "-":
        #Only errors are logged
        logging.basicConfig(level=logging.ERROR, format='%(name)-12s: %(levelname)-8s %(message)s')
    elif args.log.lower() == "stdout":
        #Everything to stdout
        logging.basicConfig(level=logging.DEBUG, format=format, datefmt='%Y-%m-%d %H:%M')
    else:
        #Log to file and warnings to screen
        logging.basicConfig(level=logging.DEBUG, format=format, datefmt='%Y-%m-%d %H:%M', filename=args.log, filemode='w')
        console = logging.StreamHandler()
        console.setLevel(logging.WARNING)
        formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

    include_genes = set()
    print("Parsing adjustments..")
    adjustments = list()
    with open(args.adjustments) as f:
        for i,line in enumerate(f,1):
            if line[0] == "#":
                continue
            line = line.split()
            if len(line) < 3:
                pass #TODO
            elif len(line) == 3:
                options = ""
            else:
                options = line[3]

            adjustments.append({"gene": line[0], "operation": line[1], "barcode_id": line[2], "options": options, "line": i})
            include_genes.add(line[0])


    print("Loading barcode file..")
    with open(args.barcodes) as bcs:
        barcoding_lib = parse_barcode_library(bcs)
    
    print("Loading config file..")
    with open(args.config) as cfg:
        config = yaml.safe_load(cfg)
        config = create_config_tables(config)
        #TODO: Validate config

    print("Parsing genome..")
    genome = SeqIO.read(args.genome, "genbank")

    print("Collecting gene list..")
    genes = seqIO_to_genelist(genome, config, include_genes)

    print("Making oligos..")
    oligos = interface(adjustments, genes, genome.seq, config, barcoding_lib, args.project)

    if args.output:
        output = args.output
    else:
        output = args.project + ".out"

    print("Writing to {}..".format(output))
    oligolist_to_tabfile(oligos, output)
