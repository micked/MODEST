#!/usr/bin/env python

from Bio import SeqIO

from mage_tool.IO import Mutation
from mage_tool import oligo_design

if __name__ == "__main__":
    genome = SeqIO.read("data/E_coli_K12_MG1655.gb", "genbank")
    
    for f in genome.features:
        if f.type == "CDS" and f.location.strand == -1:
            o = f
            break

    print o.location.start
    print o.extract(genome).seq
    print
    print genome[o.location.start-6:o.location.end+6].seq.reverse_complement()
    print
    leader = genome.seq[o.location.end:o.location.end+35].reverse_complement()
    print leader