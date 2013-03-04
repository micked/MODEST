#!/usr/bin/env python

"""
Mutation generation routines related to translation initiation
"""


from mutation_tools import find_mutation_box
from mutation_tools import compare_seqs
from RBS.RBS_Calculator import RBS_Calculator
from oligo_design import Mutation


def replace_start_codon(gene, start_codon="ATG"):
    """Replace start codon"""

    if str(gene.cds[0:3]).upper() == str(start_codon).upper():
        return False
    if len(start_codon) != 3:
        raise ValueError("start_codon must be 3 nucleotides long. Recieved {} ({})".format(start_codon, len(start_codon)))
    
    mutation = find_mutation_box(gene.cds[0:3], start_codon)

    return gene.do_mutation(mutation)


def translational_KO(gene, stop_codons=["TAG", "TAA", "TGA"], KO_frame=10):
    """Knock out a gene with a stop-codon with the least possible mutations

    KO_frame is how many codons the method will look for stop codon mutations.

    This method will prioritise single mutations, double mutations located
    beside each other, double mutations with a match inbetween, and finally
    triple mutations.
    """
    start_offset = 3
    KO_frame = min(KO_frame, (len(gene.cds)-start_offset)/3)

    KO = gene.cds[start_offset:KO_frame*3+start_offset]
    #m = target mutations, g = target groups
    for m,g in [(1,1), (2,1), (2,2), (3,1)]:
        #Loop codons
        for i in range(0, len(KO), 3):
            #Parent codon
            parent = KO[i:i+start_offset]
            for child in stop_codons:
                #test mutitions and test groups
                tm,tg = compare_seqs(parent, child)
                if tm <= m and tg <= g:
                    #Do mutation
                    mutation = find_mutation_box(parent, child)
                    mutation.pos += i+start_offset
                    #Apply mutation
                    new_mut = gene.do_mutation(mutation)
                    new_mut._codon_offset = (i+start_offset)/3
                    return new_mut


def RBS_single_mutation(gene, maximise=True, insert=False, delete=False, top=3):
    """Maxi/minimise RBS binding with single mutations (Bruteforce approach)"""
    mutation_list = list()
    
    seq = gene.leader + gene.cds
    RBS_calc = RBS_Calculator(str(seq), [35, 35], "")
    RBS_calc.calc_dG()
    old_value = RBS_calc.calc_expression_level(RBS_calc.dG_total_list[0])
    mutation_list.append((old_value, 0, ""))

    mutations = {"A": ["T", "G", "C"],
                 "T": ["A", "G", "C"],
                 "G": ["A", "T", "C"],
                 "C": ["A", "T", "G"]
    }

    if delete:
        for N in mutations:
            mutations[N].append("")

    if insert:
        for N in mutations:
            for M in ["A", "T", "G", "C"]:
                mutations[N].append(N+M)

    leader = str(gene.leader)

    for i,n in enumerate(leader.upper()):
        for m in mutations[n]:
            new_leader = leader[:i] + m + leader[i+1:]
            new_seq = new_leader + str(gene.cds)
            RBS_calc = RBS_Calculator(new_seq, [34+len(m), 34+len(m)], "")
            RBS_calc.calc_dG()
            tmp_mut = "{}={}".format(n, m)
            if not len(RBS_calc.dG_total_list):
                continue
            new_value = RBS_calc.calc_expression_level(RBS_calc.dG_total_list[0])
            mutation_list.append((new_value, i, tmp_mut))

    muts = list()
    for m in sorted(mutation_list, key=lambda x: x[0], reverse=maximise)[0:top]:
        if not m[2]:
            break
        mut = Mutation("eq", m[2], -35+m[1])
        mut._adjustment = m[0] / old_value
        muts.append(gene.do_mutation(mut))

    return muts