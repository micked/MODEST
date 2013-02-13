#!/usr/bin/env python

from Bio import SeqIO

from mage_tool.IO import Mutation
from mage_tool import oligo_design 

if __name__ == "__main__":
    ref_genome = SeqIO.read("data/E_coli_K12_MG1655.gb", "genbank")
    
    mutation_list = [l.strip().split() for l in list(open("data/all_mutations.txt", "U"))]
    
    project = "HeatAdaptation"
    
    with open("master.fasta", "w") as f:
        for mut in mutation_list:
            mutation = Mutation("arrow", mut[1], mut[0], mut[2], ref_genome)
            oligo = oligo_design.mut_to_oligo(mutation, ref_genome)
            barcoded = oligo_design.barcoding(oligo, "a", "a")
            header = ">" + project + "_" + "barcode:master_" +  mut[2]  + ":" + mut[1] + "_pos:" + str(mutation.position)
            f.write(header + "\t" + str(barcoded) + "\n")
        
