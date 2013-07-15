#!/usr/bin/env python2

"""
Make a genome config file from parameters.
"""

from __future__ import print_function

import os
import argparse
from urllib2 import HTTPError

import yaml
from Bio import Entrez
from Bio import SeqIO
#Entrez.email = 'sch@ntz.nu'

from mage_tool.IO import make_genomeconfig
from mage_tool.IO import load_config_file


def int_or_region(val):
    try:
        val = int(val)
        return [val-100, val+100]
    except (TypeError, ValueError):
        return [int(r) for r in val.split(',')[0:2]]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    gb_group = parser.add_mutually_exclusive_group(required=True)
    gb_group.add_argument('--gb', help='GenBank genome file', type=argparse.FileType('U'))
    gb_group.add_argument('--nc', help='NCBI Genbank identifier to automatically fetch genome online.')
    parser.add_argument('-o', help='Origin of replication, region eg. 100,200 or integer. Can be left out if \'rep_origin\' is annotated.', default=None, type=int_or_region)
    parser.add_argument('-t', help='Termination of replication. Region or integer.', default=None, type=int_or_region)
    parser.add_argument('--output', help='Output dir or filename. Defaults to location of genome, current dir if genome is downloaded.', default=None)
    args = parser.parse_args()

    #In current working directory
    if args.output is None:
        args.output = os.getcwd()

    #If output is specifed as a directory
    if os.path.isdir(args.output):
        output_dir = args.output
        gbcfg_filename = None
    #Or a file
    else:
        output_dir = os.path.dirname(args.output)
        gbcfg_filename = args.output

    #Fetch it weth efetch
    if args.nc:
        print("Fetching:", args.nc)
        try:
            handle = Entrez.efetch(db='nuccore', id=args.nc, rettype='gbwithparts', retmode='text')
        except HTTPError:
            print('Could not find', args.nc)
            exit(1)

        #Read genome
        genome = SeqIO.read(handle, 'genbank')

        #Lets use an expressive name
        gb_basename = '{}-{}'.format(genome.id, genome.annotations['source'].replace(' ', '_'))
        gb_filename = os.path.join(output_dir, gb_basename) + '.gb'

        #Write out genome
        SeqIO.write(genome, gb_filename, 'genbank')
    #Load local genome
    else:
        genome = SeqIO.read(args.gb, 'genbank')
        gb_filename = args.gb.name
        gb_basename = os.path.splitext(args.gb.name)[0]

    #Set gbcfg filename if the user hasn't
    if not gbcfg_filename:
        gbcfg_filename = os.path.join(output_dir, gb_basename) + '.gbcfg'

    print('Genome file:       ', gb_filename)
    print('Genome config file:', gbcfg_filename)

    gbcfg = make_genomeconfig(genome, args.o, args.t)

    with open(gbcfg_filename, 'w') as f:
        f.write(yaml.dump(gbcfg))
        print(yaml.dump(gbcfg))
