#!/usr/bin/env python2

"""
Commandline interface to <>

Yadda yadda yadda
"""

from __future__ import print_function
import logging
import argparse
import os.path
import sys
import re

from Bio import SeqIO

from mage_tool.IO import ParserError
from mage_tool.IO import seqIO_to_genelist
from mage_tool.IO import oligolist_to_tabfile
from mage_tool.IO import parse_barcode_library
from mage_tool.IO import load_config_file
from mage_tool.IO import oligolist_to_csv
from mage_tool.IO import OligoLibraryReport
from mage_tool.IO import oligolist_to_mascfile
from mage_tool.IO import raw_adjlist_to_adjlist
from mage_tool.IO import generate_codon_usage
from mage_tool.IO import create_config_cdntables
from mage_tool.interface import parse_adjustments
from mage_tool.interface import run_adjustments
import mage_tool.run_control as rc

#Define a log
log = logging.getLogger("MODEST.py")


class ListOperations(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if namespace.rc:
            rc.load_rcfile(namespace.rc)

        from mage_tool.operations import load_operations, OPERATIONS
        load_operations()
        for op in OPERATIONS:
            pass
            print(op, OPERATIONS[op].__doc__)

        if not namespace.rc:
            print('', 'To load an rc file before displaying operations, set --rc before --operations.')

        parser.exit()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("adjustments", help="Adjustment list", type=argparse.FileType("U"))
    parser.add_argument("config", help="Genome configuration", type=argparse.FileType("U"))
    parser.add_argument("barcodes", help="Barcode library", nargs="?", default=None, type=argparse.FileType("U"))
    parser.add_argument("--genome", help="Annotated genome. Leave empty to locate automatically.", default=None)
    loghelp = "Logfile, default <project>.log. Use STDOUT to log all messages to screen, use - to disable"
    parser.add_argument("--log", help=loghelp, default="--")
    parser.add_argument("--rc", help='MODEST conf file', default=None)
    parser.add_argument("-p", "--project", help="Project name", default="Untitled")
    parser.add_argument("-o", "--output", help="Output file. Default <project>.out", default=False)
    parser.add_argument("-T", help="Run unthreaded", action="store_true")
    parser.add_argument("--MASC", help="Design MASC PCR primers to file", default=False, nargs="?", metavar="mascfile")
    parser.add_argument("--PDF", help="Output report PDF", default=False, nargs="?", metavar="pdffile")
    parser.add_argument("--CSV", help="Output oligos in csv format", default=False, nargs="?", metavar="csvfile")
    parser.add_argument("--HTML", help="Output report HTML", default=False, nargs="?", metavar="htmlfile")
    parser.add_argument("--operations", help="Display registered operations and exit", action=ListOperations, nargs=0)
    args = parser.parse_args()

    #Load MODEST configuration
    if args.rc:
        rc.load_rcfile(args.rc)

    #Custom printed logging level
    LOG_OUTPUT=55
    logging.addLevelName(LOG_OUTPUT, 'OUTPUT')

    def modestprint(self, message, *args, **kws):
        self.log(LOG_OUTPUT, message, *args, **kws)

    logging.Logger.modestprint = modestprint

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
        logging.basicConfig(level=logging.DEBUG, format=format, datefmt='%Y-%m-%d %H:%M', filename=args.log)

        console = logging.StreamHandler()
        console.setLevel(logging.WARNING)
        formatter = logging.Formatter(
            '%(name)-12s: %(levelname)-8s %(message)s')
        console.setFormatter(formatter)
        logging.getLogger('').addHandler(console)

    log.modestprint("Reading adjustments..")
    try:
        adjlist = raw_adjlist_to_adjlist(args.adjustments.readlines())
    except ParserError as e:
        print(str(e))
        for error in e.error_list:
            print ("[E] - " + error)
        exit(1)

    include_genes = set([l["gene"] for l in adjlist])

    log.modestprint("Loading config file..")
    cfg_basedir = os.path.abspath(os.path.dirname(args.config.name))

    try:
        config = load_config_file(args.config, cfg_basedir)
    except ParserError as e:
        print(str(e))
        for error in e.error_list:
            print("[E] -", error)
        exit(1)

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

    log.modestprint("Parsing genome ({})..".format(genome_loc))
    genome = SeqIO.read(genome_loc, "genbank")

    #Make ekstra config tables
    if 'codon_usage' not in config:
        log.modestprint('Generating codon usage..')
        config = generate_codon_usage(genome, config)
    config = create_config_cdntables(config)

    log.modestprint("Collecting gene list..")
    incl_gnm = "genome" in include_genes
    #include_genes = None if "all" in include_genes
    try:
        genes = seqIO_to_genelist(genome, config, include_genes, include_genome=incl_gnm)
    except ParserError as e:
        print("Fatal error:", e)
        exit(1)

    barcoding_lib = dict()
    if args.barcodes:
        log.modestprint("Loading barcode file..")
        barcoding_lib = parse_barcode_library(args.barcodes)

    log.modestprint("Making oligos..")
    threaded = not args.T
    try:
        oplist = parse_adjustments(adjlist, genes, config, barcoding_lib)
    except ParserError as e:
        print(str(e))
        for error in e.error_list:
            print("[E] -", error)
        exit(1)

    oligos = run_adjustments(oplist, genome.seq, args.project, barcoding_lib, threaded)

    if args.output:
        output = args.output
    else:
        output = args.project + ".out"

    log.modestprint("Writing to {}..".format(output))
    oligolist_to_tabfile(oligos, output)

    if args.CSV or args.CSV is None:
        output_csv = args.project + ".csv" if not args.CSV else args.CSV
        log.modestprint("Writing report CSV to {}..".format(output_csv))
        csvlist = oligolist_to_csv(oligos, output_csv)
    else:
        csvlist = oligolist_to_csv(oligos)

    if args.PDF or args.PDF is None or args.HTML or args.HTML is None:
        report = OligoLibraryReport(args.project)
        report.parse_and_generate(csvlist, csv_file=False)

    if args.PDF or args.PDF is None:
        output_pdf = args.project + ".pdf" if not args.PDF else args.PDF
        log.modestprint("Writing report PDF to {}..".format(output_pdf))
        report.write_pdf(output_pdf)

    if args.HTML or args.HTML is None:
        output_html = args.project if not args.HTML else args.HTML
        output_html = re.sub(r'\.html$', '', output_html, flags=re.IGNORECASE)
        log.modestprint("Writing report PDF to {}.html..".format(output_html))
        report.write_html(output_html)

    if args.MASC or args.MASC is None:
        mascfile = args.MASC if args.MASC else args.project + "_masc_primers.csv"
        log.modestprint("Designing MASC PCR primers to {}..".format(mascfile))
        masc_kwargs = {"lengths": [100, 150, 200, 250, 300, 400, 500, 600, 700, 850]}
        masc_kwargs["ref_genome"] = genome.seq
        masc_kwargs["temp"] = 62.0
        masc_primers = oligolist_to_mascfile(oligos, masc_kwargs, mascfile)

    log.info('Job done')
