
"""
Classes dealing with input/output and format changes
"""

import logging
import csv
import re
import codecs

from helpers import reverse_complement
from helpers import reverse_complement_dgn
from helpers import seqs_to_degenerate
from helpers import cds_to_wobble
from oligo_design import Gene

log = logging.getLogger("MODEST.IO")
log.addHandler(logging.NullHandler())


def seqIO_to_genelist(genome, options, include_genes=None, leader_len=35):
    """TODO"""

    genes = dict()

    for g in genome.features:
        if g.type == "CDS":
            name = g.qualifiers["gene"][0]

            if include_genes and name not in include_genes:
                continue

            #Little bit of warning:
            #Bio converts positions to 0-index internally
            #This is mostly a good thing
            strand = g.location.strand
            cds = g.extract(genome).seq
            if strand == 1:
                pos = g.location.start
                #leader = genome.seq[pos-leader_len:pos]
            else:
                pos = g.location.end
                #leader = reverse_complement(genome.seq[pos:pos+leader_len])

            l_start = min(pos, pos-leader_len*strand)
            l_end = max(pos, pos-leader_len*strand)

            #Leader in genome context (+1)
            leader = genome.seq[l_start:l_end]

            leader_wobble = None
            #Look for CDS in leader
            for c in genome.features:
                    if c.type == "CDS":
                        if c.location.start < l_start < c.location.end or \
                           c.location.start < l_end < c.location.end or \
                           (c.location.start > l_start and c.location.end < l_end):
                            w_cds = c.extract(genome).seq
                            # w_seq = w_cds
                            w_seq = cds_to_wobble(w_cds, options)

                            #Wobble sequence in genome context (+1)
                            if c.location.strand == -1:
                                w_seq = reverse_complement_dgn(w_seq)

                            start_offset = l_start - c.location.start
                            end_offset = l_end - c.location.start

                            st = ""
                            ed = ""
                            if start_offset < 0:
                                # st = str(genome.seq[l_start:l_start-start_offset])
                                st = "N"*-start_offset
                                start_offset = 0
                            if end_offset > len(w_seq):
                                #DOUBLE CHECK!!!
                                # ed = str(genome.seq[l_start+start_offset+len(w_seq):l_end])
                                ed = "N"*(end_offset-len(w_seq))
                                end_offset = len(w_seq)

                            leader_wobble = st + str(w_seq[start_offset:end_offset]) + ed

            if strand == -1:
                leader = reverse_complement(leader)
                if leader_wobble:
                    leader_wobble = reverse_complement_dgn(leader_wobble)

            #TODO
            promoter = None
            promoter_pos = None

            if name in genes:
                log.warn("Gene {} found more than once.".format(name))
                #raise Exception("Gene {} found twice!".format(name))
            else:
                genes[name] = Gene(name, pos, strand, cds, leader, leader_wobble, promoter, promoter_pos)

    return genes


def parse_barcode_library(barcode_filehandle):
    """Parse an open barcoding library file.

    Barcoding library file has the following syntax:
    ..

    """

    primers = dict()
    barcodes = dict()

    primer_flag = False
    barcodes_flag = False

    for line in barcode_filehandle:
        if line[0] == "#" or not line.strip():
            continue
        elif line.strip() == ">PRIMERS":
            primer_flag = True
            continue
        elif line.strip() == ">LIBRARY":
            barcodes_flag = True
            continue
        
        if barcodes_flag:
            line = line.split()
            forward = primers[line[1]]
            reverse = reverse_complement(primers[line[2]])
            barcodes[line[0]] = {"forward": forward, "reverse": reverse}
        elif primer_flag:
            line = line.split()
            primers[line[0]] = line[1]
        else:
            raise Exception("Primer header not found")

    return barcodes


def create_config_tables(config):
    """Create additional lookup tables from the parsed config yaml file."""
    #Translation table
    AAtable = dict() #{AA: [] for AA in amino_acids}
    config["codon_table"] = dict()
    config["stop_codons"] = list()
    for triplet, (AA, usage) in config["codon_usage"].items():
        config["codon_table"][triplet] = AA
        if AA in AAtable:
            AAtable[AA].append(triplet)
        else:
            AAtable[AA] = [triplet]
        if AA == "*":
            config["stop_codons"].append(triplet)

    #Reverse translation table
    config["dgn_table"] = dict()
    for AA in AAtable:
        config["dgn_table"][AA] = seqs_to_degenerate(AAtable[AA])

    return config

"""
Oligolist functions
"""

def oligolist_to_tabfile(oligolist, output):
    """Writeout oligolist to a tab separeted file.

    Supply string filename or open file
    
    """
    cls = False
    if not hasattr("write", output):
        output = open(output, "w")
        cls = True

    for o in sorted(oligolist, key=lambda x:x.number):
        output.write(o.id() + "\t" + o.output() + "\n")

    if cls:
        output.close()


def oligolist_to_report(oligolist, output):
    """Print report from oligolist in CSV format."""
    csvoutlist = list()

    op_reg = re.compile(r"\[(\w+)/(\w+)\] (line \d+)\s?(\S*)")

    for o in sorted(oligolist, key=lambda x:x.number):
        m = re.match(op_reg, o.operation)
        op_cmd, gene, line, op_options = m.groups()
        csvoutlist.append([o.short_id(), op_cmd, gene, line, op_options, str(o.mut), "+".join(o.barcode_ids), o.output()] + o.operation_values)

    cls = False
    if not hasattr("write", output):
        output = open(output, "w")
        cls = True

    #For Libreoffice Calc
    output.write(codecs.BOM_UTF8)
    csv_w = csv.writer(output)

    for n in csvoutlist:
        csv_w.writerow(n)

    if cls:
        output.close()
