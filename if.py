#!/usr/bin/env python

from Bio import SeqIO

from mage_tool.IO import Mutation
from mage_tool import oligo_design

if __name__ == "__main__":
    #ref_genome = SeqIO.read("data/E_coli_K12_MG1655.gb", "genbank")
    rg2 = open("data/GenomeAncestor_no_header_full2.txt").read()
    rg2 = rg2.strip("_")
    
    mutation_list = [l.strip().split() for l in list(open("data/some_mutations.txt", "U"))]
    
    project = "HeatAdaptation"

    #"""
    with open("master.fasta", "w") as f:
        for mut in mutation_list:
            mutation = Mutation("arrow", mut[1], mut[0], mut[2], rg2)
            oligo = oligo_design.mut_to_oligo(mutation, rg2)
            oligo = oligo_design.target_lagging_strand(oligo, mutation.position, ori=3886229, ter=1493223)
            barcoded = oligo_design.barcoding(oligo, "a", "a")
            header = ">" + project + "_" + "barcode:master_" +  mut[2]  + ":" + mut[1] + "_pos:" + str(mutation.position)
            f.write(header + "\t" + str(barcoded) + "\n")
    #"""