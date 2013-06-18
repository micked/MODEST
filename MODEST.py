#!/usr/bin/env python2

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

from mage_tool.IO import ParserError
from mage_tool.IO import seqIO_to_genelist
from mage_tool.IO import oligolist_to_tabfile
from mage_tool.IO import parse_barcode_library
from mage_tool.IO import create_config_tables
from mage_tool.IO import oligolist_to_csv
from mage_tool.IO import OligoLibraryReport
from mage_tool.IO import oligolist_to_mascfile
from mage_tool.IO import raw_adjlist_to_adjlist
from mage_tool.interface import parse_adjustments
from mage_tool.interface import run_adjustments


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("adjustments", help="Adjustment list", type=argparse.FileType("U"))
    parser.add_argument("barcodes", help="Barcode library", type=argparse.FileType("U"))
    parser.add_argument("config", help="Genome configuration", type=argparse.FileType("U"))
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
    parser.add_argument("--MASC", help="Design MASC PCR primers to file", default=False, nargs="?", metavar="mascfile")
    parser.add_argument("--PDF", help="Output report PDF", default=False, nargs="?", metavar="pdffile")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-Y", action="store_true", help="Continue despite errors")
    group.add_argument("-N", action="store_true", help="Exit on errors")
    parser.add_argument("--operations", help="Display registered operations and exit", action="store_true")
    args = parser.parse_args()

    if args.operations:
        from mage_tool.operations import OPERATIONS
        for op in OPERATIONS:
            print(op, OPERATIONS[op].__doc__)
        exit(0)

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
        logging.basicConfig(level=logging.DEBUG, format=format, datefmt='%Y-%m-%d %H:%M')
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
    adjustlist = args.adjustments.readlines()

    include_genes = set([line.split()[0] for line in adjustlist
                         if line.strip() and line[0] != "#"])

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

    print("Parsing genome ({})..".format(genome_loc))
    genome = SeqIO.read(genome_loc, "genbank")

    print("Collecting gene list..")
    incl_gnm = "genome" in include_genes
    #include_genes = None if "all" in include_genes
    try:
        genes = seqIO_to_genelist(genome, config, include_genes,
                                  include_genome=incl_gnm)
    except ParserError as e:
        print("Fatal error:", e)
        exit(1)

    print("Loading barcode file..")
    #with open(args.barcodes) as bcs:
    barcoding_lib = parse_barcode_library(args.barcodes)

    print("Making oligos..")
    threaded = not args.T
    adjlist, error_list = raw_adjlist_to_adjlist(adjustlist)
    oplist, op_error_list = parse_adjustments(adjlist, genes, config, barcoding_lib)
    error_list.extend(op_error_list)
    if error_list:
        print("Following errors found:")
        for e in error_list:
            print("[E]:", e)

        if args.Y:
            cont = "y"
        elif args.N:
            cont = "n"
        else:
            cont = raw_input("\nContinue anyways? y/n [n]: ")

        if not cont or cont.lower()[0] is "n":
            print("Computation stoppet prematurely.")
            exit(1)

    adj_args = (adjustlist, genes, genome.seq, config, args.project,
                barcoding_lib, threaded)
    oligos = run_adjustments(oplist, genome.seq, args.project, barcoding_lib, threaded)

    #if errors:
    #    print("Computation not started, errors in adjustments file:")
    #    print("\n".join(errors))
    #    exit(1)

    if args.output:
        output = args.output
    else:
        output = args.project + ".out"

    print("Writing to {}..".format(output))
    oligolist_to_tabfile(oligos, output)

    output_csv = args.project + ".csv"
    print("Writing report CSV to {}..".format(output_csv))
    csvlist = oligolist_to_csv(oligos, output_csv)

    if args.PDF or args.PDF is None:
        output_pdf = args.project + ".pdf" if not args.PDF else args.pdf
        print("Writing report PDF to {}..".format(output_pdf))
        report = OligoLibraryReport(args.project)
        report.parse_and_generate(csvlist, csv_file=False)
        report.write_pdf(output_pdf)

    if args.MASC or args.MASC is None:
        mascfile = args.MASC if args.MASC else args.project + "_masc_primers.csv"
        print("Designing MASC PCR primers to {}..".format(mascfile))
        masc_kwargs = {"lengths": [100, 150, 200, 250, 300, 400, 500, 600, 700, 850]}
        masc_kwargs["ref_genome"] = genome.seq
        masc_kwargs["temp"] = 62.0
        masc_primers = oligolist_to_mascfile(oligos, masc_kwargs, mascfile)
