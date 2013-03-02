
"""
Classes dealing with input/output and format changes
"""

import logging

from DNA_tools import reverse_complement
from oligo_design import Gene

log = logging.getLogger("MODEST.IO")
log.addHandler(logging.NullHandler())


def seqIO_to_genelist(genome, include_genes=None):
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
                leader = genome.seq[pos-35:pos]
            else:
                pos = g.location.end
                leader = reverse_complement(genome.seq[pos:pos+35])
            #TODO
            promoter = None
            promoter_pos = None

            if name in genes:
                log.warn("Gene {} found more than once.".format(name))
                #raise Exception("Gene {} found twice!".format(name))

            genes[name] = Gene(name, pos, strand, cds, leader, promoter, promoter_pos)

    return genes

def oligolist_to_tabfile(oligolist, output):
    """Writeout oligolist to a tab separeted file

    Supply string filename or open file
    """
    cls = False
    if not hasattr("write", output):
        output = open(output, "w")
        cls = True

    for o in oligolist:
        output.write(o.id() + "\t" + o.output() + "\n")

    if cls:
        output.close()