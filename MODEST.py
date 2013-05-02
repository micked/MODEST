#!/usr/bin/env python

"""
Commandline interface to <>

Yadda yadda yadda
"""
from __future__ import print_function
import logging
import argparse
import os.path

import yaml
from Bio import SeqIO

from mage_tool.interface import run_adjustments
from mage_tool.interface import run_adjustments_unthreaded
from mage_tool.IO import seqIO_to_genelist
from mage_tool.IO import oligolist_to_tabfile
from mage_tool.IO import parse_barcode_library
from mage_tool.IO import create_config_tables
from mage_tool.IO import oligolist_to_csv
from mage_tool.IO import OligoLibraryReport


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("adjustments", help="Adjustment list")
    parser.add_argument("barcodes", help="Barcode library")
    parser.add_argument("config", help="Genome configuration")
    parser.add_argument("--genome", help="Annotated genome. Leave empty to"
                        " locate automatically.", default=None)
    parser.add_argument("--log", help="Logfile, default <project>.log. "
                        "Use STDOUT to log all messages to screen, "
                        "use - to disable", default="--")
    parser.add_argument("-p", "--project", help="Project name",
                        default="Untitled")
    parser.add_argument("-o", "--output", help="Output file. "
                        "Default <project>.out", default=False)
    parser.add_argument("-T", help="Run unthreaded", action="store_true")
    args = parser.parse_args()

    if args.log == "--":
        args.log = args.project + ".log"

    # Set up logging
    format = "%(asctime)s %(name)-12s: %(levelname)-8s %(message)s"
    if args.log == "-":
        # Only errors are logged
        logging.basicConfig(level=logging.ERROR, format="%(name)-12s: "
                            "%(levelname)-8s %(message)s")
    elif args.log.lower() == "stdout":
        # Everything to stdout
        logging.basicConfig(
            level=logging.DEBUG, format=format, datefmt='%Y-%m-%d %H:%M')
    else:
        # Log to file and warnings to screen
        logging.basicConfig(level=logging.DEBUG, format=format,
                            datefmt='%Y-%m-%d %H:%M',
                            filename=args.log, filemode='w')

        console = logging.StreamHandler()
        console.setLevel(logging.WARNING)
        formatter = logging.Formatter(
            '%(name)-12s: %(levelname)-8s %(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

    print("Reading adjustments..")
    with open(args.adjustments) as f:
        adjustlist = f.readlines()

    include_genes = set([line.split()[0] for line in adjustlist
                         if line.strip() and line[0] != "#"])

    print("Loading config file..")
    with open(args.config) as cfg:
        config = yaml.safe_load(cfg)
        config = create_config_tables(config)

    if not args.genome:
        genome_locations = [os.path.splitext(args.config)[0] + ".gb",
                            os.path.splitext(args.config)[0] + ".genbank",
                            config["Locus"],
                            config["Locus"] + ".gb",
                            config["Locus"] + ".genbank",
                            os.path.join(os.path.dirname(
                                args.config), config["Locus"]),
                            os.path.join(os.path.dirname(
                                args.config), config["Locus"] + ".gb"),
                            os.path.join(os.path.dirname(args.config),
                                         config["Locus"] + ".genbank")]

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

    print("Parsing genome ({})..".format(genome_loc))
    genome = SeqIO.read(genome_loc, "genbank")

    print("Collecting gene list..")
    genes = seqIO_to_genelist(genome, config, include_genes)

    print("Loading barcode file..")
    with open(args.barcodes) as bcs:
        barcoding_lib = parse_barcode_library(bcs)

    print("Making oligos..")
    adj_args = (adjustlist, genes, genome.seq, config, args.project, barcoding_lib)
    if not args.T:
        oligos, errors = run_adjustments(*adj_args)
    else:
        oligos, errors = run_adjustments_unthreaded(*adj_args)
    if errors:
        print("Computation not started, errors in adjustments file:")
        print("\n".join(errors))
        exit(1)

    if args.output:
        output = args.output
    else:
        output = args.project + ".out"

    print("Writing to {}..".format(output))
    oligolist_to_tabfile(oligos, output)

    output_csv = args.project + ".csv"
    print("Writing report CSV to {}..".format(output_csv))
    csvlist = oligolist_to_csv(oligos, output_csv)

    output_pdf = args.project + ".pdf"
    print("Writing report PDF to {}..".format(output_pdf))
    report = OligoLibraryReport(args.project)
    report.parse_and_generate(csvlist, csv_file=False)
    report.write_pdf(output_pdf)
