#!/usr/bin/env python

"""
Mutation generation routines related to promoters
"""

from __future__ import print_function

import math
import logging

from oligo_design import Gene
from oligo_design import Sequence

def biggest_difference(promoter, promoterdict, downregulate = False):
    """Calculate single mutation that creates the biggest increase or decrease of promoter activity. Returns None if no mutation is found"""
    biggest_difference = 0
    biggest_pos = -9
    biggest_nt = ""
    
    for i, nt in enumerate(promoter):
        for nts in promoterdict[i]:
            if (not downregulate and promoterdict[i][nt]-promoterdict[i][nts] > biggest_difference) or (downregulate and (promoterdict[i][nt]-promoterdict[i][nts] < biggest_difference)):
                biggest_difference = promoterdict[i][nt]-promoterdict[i][nts]
                biggest_pos = i
                biggest_nt = nts

    newseq = "{}{}{}".format(promoter[:biggest_pos], biggest_nt, promoter[biggest_pos+1:])
    if newseq == promoter:
        return None        
            
    return promoter.mutate(biggest_nt, biggest_pos, in_place=True)
    

def closest_to_target(promoter, promoterdict, target_energy):
    """Calculate single mutation that is closest to target energy. Returns None if no mutation is found."""
    closest_difference = None
    closest_pos = -9
    closest_nt = ""
    
    for i, nt in enumerate(promoter):
        for nts in promoterdict[i]:
            newseq = "{}{}{}".format(promoter[:i], nts, promoter[i+1:])
            newseq_energy = -promoter_energy(newseq, promoterdict)
            if not closest_difference or math.fabs(target_energy-newseq_energy) < closest_difference:
                closest_difference = math.fabs(target_energy-newseq_energy)
                closest_pos = i
                closest_nt = nts
                
    newseq = "{}{}{}".format(promoter[:closest_pos], closest_nt, promoter[closest_pos+1:])
    if newseq == promoter:
        return None    
          
    return promoter.mutate(closest_nt, closest_pos, in_place=True)
    

def promoter_energy(promoter, matrixdict):
    """Calculate promoter energy from matrix"""
    energy = 0
    if not promoter:
        print("wha")
    for i, nt in enumerate(promoter):
        energy += matrixdict[i][nt]
    
    return energy
    
def promoter_library(gene, targets, max_mutations, matrix):
    """Create a library of promoters.
    
    targets is a list of targets. They can be specified as a fold, ie. 2, 5 or 10 fold, specified as 2, 5 or 10 respectively. Or as max or min for finding the strongest/weakest promoter. One promoter will be returned per target.
    
    max_mutations is how many mutations are allowed.
    
    matrix is the energy matrix for the specified promoter.
    
    ."""
    if str(gene) == "genome":
        log.error("Cannot use Promoter_library on genome")
        return None
        
    #list of promoters to return
    promoterlist = list()
    
    
    for target in targets:
        newseq = gene.promoter
        newseq = newseq.copy()
        muts = 0
        if target == "max":
            for i in range(max_mutations):
                tmpseq = biggest_difference(newseq, matrix)
                if tmpseq:
                    newseq = tmpseq
                else:
                    break  
        elif target == "min":
            for i in range(max_mutations):
                tmpseq = biggest_difference(newseq, matrix, True)
                if tmpseq:
                    newseq = tmpseq
                else:
                    break
        else:
            target_energy = -promoter_energy(newseq, matrix)+math.log(target)
            for i in range(max_mutations):
                tmpseq = closest_to_target(newseq, matrix, target_energy)
                if math.exp(target_energy)*0.95 <= math.exp(-promoter_energy(newseq, matrix)) <= math.exp(target_energy)*1.05 or not tmpseq:
                    break
                newseq = tmpseq
        fold_regulation = math.exp(-promoter_energy(newseq, matrix))/math.exp(-promoter_energy(gene.promoter, matrix))
        print(fold_regulation)
        promoterlist.append([newseq, fold_regulation])
        
    muts = list()
    for promoter in promoterlist:
        mut = promoter[0].get_mutation()
        mut.pos -= gene.promoter_pos
        mut = gene.do_mutation(mut)
        mut._fold = promoter[1]
        mut._n = promoter[0].n_muts()
        muts.append(mut)

    return muts
    
    
if __name__ == "__main__":
    
    import yaml
    
    with open("ymldump") as f:
        test = yaml.load(f)

    
    for prom in test["promoters"]:
        for i, l in enumerate(test["promoters"][prom]):
            ntdict = dict()
            ntdict["A"],ntdict["C"],ntdict["G"],ntdict["T"] = float(l[0]), float(l[1]), float(l[2]), float(l[3])
            test["promoters"][prom][i] = ntdict

    gene = Gene("test", 5, 1, "TCGAGTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGG", leader=None, promoter="TCGAGTTTACACTTTATGCTTCCGGCTCGTATAATGTGTGG",
                 promoter_pos=0)
                 
    print(" {}".format(gene.promoter))
    lib = promoter_library(gene, list([2, 4]), 5, test["promoters"]["sigma70"])
    
    print(lib)
    
