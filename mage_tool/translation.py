#!/usr/bin/env python

"""
Mutation generation routines related to translation initiation
"""

from __future__ import print_function

import math
import random
import logging

from helpers import dgn_to_nts
from helpers import valid_rna
from oligo_design import Mutation
from mutation_tools import compare_seqs
from mutation_tools import find_mutation_box
from ViennaRNA import ViennaRNA
from ViennaRNA import brackets_to_basepairing
from ViennaRNA import basepairing_to_brackets

#Define a log
log = logging.getLogger("MODEST")
log.addHandler(logging.NullHandler())

def replace_start_codon(gene, start_codon="ATG"):
    """Replace start codon"""

    if str(gene.cds[0:3]).upper() == str(start_codon).upper():
        return False
    if len(start_codon) != 3:
        raise ValueError("start_codon must be 3 nucleotides long. Recieved: "
                         "{} ({})".format(start_codon, len(start_codon)))

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
                    #print(start_offset, i, parent, child, mutation, gene.cds[0:50])
                    mutation.pos += i + start_offset
                    #Apply mutation
                    new_mut = gene.do_mutation(mutation)
                    new_mut._codon_offset = (i+start_offset)/3
                    return new_mut


def RBS_library_fuzzy(gene, target, n, max_mutations, low_count=0):
    """Generate a fuzzy RBS library.

    'fuzzy' means a single Monte Carlo simulation is run, and a library is
    collected along the way. This method has a tendency to be biased towards a
    higher mutation count or higher expression level.

    target is the target translation initiation rate in AU, n is number of RBS
    sequences to output, and max_mutations is the maximum number of mutations
    to perform.

    low_count can be set to prioritise sequences with number of mutations that
    are below this number. Ie. at low_count=6, the best sequences with 1 to 6
    (included) mutations are added before the rest of the library is selected.

    """
    #Convert target to dG
    target_dG = AU_to_dG(target)

    seq_lib = RBS_Monte_Carlo(gene, target_dG, max_mutations, True)
    seq_lib.sort(key=lambda seq: seq.delta)

    selected_lib = list()

    #Select best sequences with low_count
    if low_count:
        for m in range(1, int(low_count)+1):
            m_muts = [seq for seq in seq_lib if seq.n_muts() == m]
            if m_muts:
                m_muts[0].low_count = True
                selected_lib.append(m_muts[0])

        selected_lib.sort(key=lambda seq: seq.delta)

    if len(selected_lib) >= n:
        seq_lib = selected_lib
        selected_lib = list()

    #Fill in the remaining spots
    n = n - len(selected_lib)
    stepsize = (seq_lib[0].delta - seq_lib[-1].delta) / n

    #Go through each leader in library and select one closest to target
    for i in range(0, n):
        step = seq_lib[0].delta - (stepsize * i)
        #print("Step:", step)
        closest = sorted(seq_lib, key=lambda l: abs(l.delta - step))
        for close in closest:
            if str(close) not in [str(l) for l in selected_lib]:
                selected_lib.append(close)
                break

    selected_lib.sort(key=lambda seq: seq.delta)

    orgAU = dG_to_AU(RBS_predict(gene.leader, gene.cds))

    muts = list()
    for seq in selected_lib:
        mut = seq.get_mutation()
        mut.pos -= len(gene.leader)
        mut = gene.do_mutation(mut)
        mut._AU = dG_to_AU(seq.dG)
        mut._orgAU = orgAU
        mut._n = seq.n_muts()
        muts.append(mut)

    return muts


def RBS_library(gene, target, n, max_mutations, m=None):
    """Create an RBS library.

    The function will attempt to reach as close to target as possible within
    the number of mutations defined in max_mutations.

    n is the number of new RBS sequences generated for the library.

    The library will reach target in m-fold steps, where m is either given as a
    parameter or calculated automatically (m=None). If m given, the function
    will return a maximum of n RBS sequences, but if target is reached before n
    sequences are generated or is unreachable, the number of RBS sequences will
    be less than n.

    """
    targetAU = target
    targetdG = AU_to_dG(target)

    pre_lib = RBS_Monte_Carlo(gene, targetdG, max_mutations, True)
    pre_lib.sort(key=lambda seq: seq.delta)

    orgdG = RBS_predict(gene.leader, gene.cds)
    orgAU = dG_to_AU(orgdG)
    bestAU = dG_to_AU(pre_lib[0].dG)

    if not m:
        m = (bestAU / orgAU)**(1./n)

    #Automatically adjust m
    if bestAU > orgAU and m < 1:
        log.debug("[{}]: m automatically changed from {} to {} because target "
                  "({}) or reached ({}) is higher than wt ({})."
                  "".format(str(gene), m, 1./m, targetAU, bestAU, orgAU))
        m = 1./m
    elif bestAU < orgAU and m > 1:
        log.debug("[{}]: m automatically changed from {} to {} because target "
                  "({}) or reached ({}) is lower than wt ({})."
                  "".format(str(gene), m, 1./m, targetAU, bestAU, orgAU))
        m = 1./m

    lib = list()

    current_target = orgAU

    for j in range(n):
        current_target *= m
        curr_dG = AU_to_dG(current_target)
        curr_delta = abs(targetdG - curr_dG)

        if curr_delta < pre_lib[0].delta:
            break

        initial = pre_lib[-1]
        for i in range(len(pre_lib), 0, -1):
            if pre_lib[i-1].delta < curr_delta: break
            initial = pre_lib[i-1]

        mutgene = gene.copy()
        initial = initial.copy()
        del(initial.dG)
        mutgene.leader = initial

        newseq = RBS_Monte_Carlo(mutgene, curr_dG, max_mutations, False)
        lib.append(newseq[0])

    if len(lib) < n:
        lib.append(pre_lib[0])

    muts = list()
    for seq in lib:
        mut = seq.get_mutation()
        mut.pos -= len(gene.leader)
        mut = gene.do_mutation(mut)
        mut._AU = dG_to_AU(seq.dG)
        mut._orgAU = orgAU
        mut._n = seq.n_muts()
        muts.append(mut)

    return muts


"""
RBS Monte Carlo functions
"""

#Constants
RBS_RT_EFF = 2.222
RBS_LOG_K =  7.824
RBS_K = 2500.0

"""
wt expression = a0
alt expression = ai
max expression = am
number of outputs = o

    Sortings: fuzzy  : quick n' dirty
              expn   : n-fold improvement
                   ai = a0*n^i (until ai > am)
                   (ignores o)
          expauto: auto-fold improvement
                   ai = a0*n^i (until i = o)
                   ao = am
                   n = (am/a0)^(1/o)

Do it as a function?
"""

def RBS_Monte_Carlo(gene, target, maxmuts=10, collect_library=False, **kwargs):
    """TODO"""

    #Start temp is a function of muts
    calc_start_temp = kwargs.get("calc_start_temp", lambda n: (0.6 / n) * 3)
    end_temp        = kwargs.get("end_temp", 0.01)

    #"Main" mutation period
    min_mut_period = kwargs.get("min_mut_period", 1000)
    #Cool mutation period
    max_mut_period = kwargs.get("max_mut_period", 4000)

    #Number of rejects before bumping to an extra mutation
    stationary = kwargs.get("stationary", 500)

    #Tolerance to target (dG)
    tolerance = kwargs.get("tolerance", 0.25)

    #Verbose mode. Handy for plots
    verbose = kwargs.get("verbose", False)

    #Number of rounds before each optimisation
    opti_rounds = kwargs.get("opti_rounds", 25)

    #RNA folder
    RNAfold = kwargs.get("RNAfold", ViennaRNA())

    #Verify target is a float
    target = float(target)

    #Proximity to target
    delta = lambda dG: abs(target - dG)

    #Initial delta
    best_dG = RBS_predict(gene.leader, gene.cds, RNAfold=RNAfold)
    best_delta = optimal_delta = delta(best_dG)

    #Initial leader
    best_leader = gene.leader.copy()

    #List leader sequences
    output = list()

    start_mutations = max(best_leader.n_muts(), 1)
    if start_mutations > maxmuts:
        raise Exception("start_mutations ({}) is larger than maxmuts ({})."
                        "".format(start_mutations, maxmuts))

    new_dG = P = total_its = 0
    def status(msg):
        if verbose:
            fmts = (msg, total_its, best_dG, new_dG, n, P)
            print("{:<18} {:>2} {:>5.2f} {:>5.2f} {:>2} {:4.2f}".format(*fmts))
            # print(msg, total_its, new_dG, new_delta)

    #Number of rejects
    rejects = 0

    #Main loop; n = number of mutations
    for n in range(start_mutations, maxmuts + 1):
        #Temperature control
        start_temp = calc_start_temp(n)
        temp_step = (start_temp - end_temp) / min_mut_period
        current_temp = start_temp

        #TODO: do ALL mutations
        for i in xrange(max_mut_period):
            #Create candidate
            candidate = best_leader.random_mutation(max_mut=n)

            #Calculate new total free energy
            new_dG = RBS_predict(candidate, gene.cds, RNAfold=RNAfold)
            if new_dG > 1e4:
                continue

            new_delta = delta(new_dG)

            accept_move = False
            if new_delta < best_delta:
                accept_move = True
                rejects = 0
                status("Accepted")
            elif new_delta == best_delta:
                rejects += 1
                if random.random() > 0.8:
                    #Randomly accept
                    accept_move = True
                    status("Randomly_accepted")
                else:
                    status("Randomly_rejected")
            else:
                try:
                    P = math.exp((best_delta - new_delta)/current_temp)
                except OverflowError:
                    status("OverflowError")
                else:
                    if P > random.random():
                        #Conditionally accept
                        accept_move = True
                        rejects = 0
                        status("Cond_accepted")
                    else:
                        rejects += 1
                        status("Cond_rejected")

            #Move is accepted
            if accept_move:
                best_leader = candidate
                best_dG = new_dG
                best_delta = new_delta

            #Library generation
            if collect_library:
                if best_delta < optimal_delta:
                    best_leader.dG = best_dG
                    best_leader.delta = best_delta
                    output.append(best_leader.copy())
                    optimal_delta = best_delta
                    #print()
                    #best_leader.pprint()

            #Target reached, STOP
            if best_delta < tolerance:
                if not hasattr(best_leader, "dG"):
                    best_leader.dG = best_dG
                    best_leader.delta = best_delta
                    output.append(best_leader.copy())
                return output

            #Optimise co-insertions/deletions
            if opti_rounds and i % opti_rounds == 0:
                best_leader.optimise_mutations()

            #Decrease temp while in "hot period"
            if i < min_mut_period:
                current_temp -= temp_step
            #Look for stationary phases in "cold period"
            elif rejects >= stationary:
                #Bump number of mutations if we are stationary
                break

            total_its += 1

    #Runout, target not reached
    if not hasattr(best_leader, "dG"):
        best_leader.dG = best_dG
        best_leader.delta = best_delta
        output.append(best_leader.copy())

    return output


def AU_to_dG(AU):
    """Convert a dG value to expression level (AU)"""
    if AU == 0:
        return 1e3 #Close enough
    return RBS_RT_EFF * (RBS_LOG_K - math.log(AU))


def dG_to_AU(dG):
    """Convert an expression level to dG"""
    return RBS_K * math.exp(-dG / RBS_RT_EFF)

"""
RBS Calculator
"""

def RBS_predict(pre_seq, cds, nt_cutoff=35, RNAfold=ViennaRNA(), verbose=False):
    """TODO: Create a more solid RBS calculator"""
    #TODO: Calculate these energies based on start tRNA
    start_codon_energies = {"AUG": -1.194,
                            "GUG": -0.0748,
                            "UUG": -0.0435,
                            "CUG": -0.03406} #hybridization to CAT
    dG_push = (12.2, 2.5, 2.0, 3.0)
    dG_pull = (0.048, 0.24, 0.0)
    # auto_dangles = True
    # dangles_default = "all"
    # temp = 37.0
    optimal_spacing = 5
    energy_cutoff = 3.0
    len_standby_site = 4

    #Footprint of the 30S complex that prevents formation of secondary
    #structures downstream of the start codon. Here, we assume that the entire
    #post-start RNA sequence does not form secondary structures once the 30S
    #complex has bound.
    footprint = 1000

    rRNA = "ACCUCCUUA"

    #TODO:
    #Figure out dangles

    #Convert to Uppercase RNA
    leader = str(pre_seq).upper().replace("T", "U")
    cds = str(cds).upper().replace("T", "U")

    if nt_cutoff:
        cds    = cds[0:nt_cutoff]
        leader = leader[-nt_cutoff:]

    if not leader or not cds:
        raise ValueError("No leader or cds supplied.")

    if not valid_rna(leader) or not valid_rna(cds):
        raise ValueError("Invalid leader or cds input.")

    if cds[0:3] not in start_codon_energies:
        raise ValueError("Invalid start codon in cds: {}".format(cds[0:3]))

    #RNAfold = ViennaRNA()

    """
    Start codon energy
    """
    dG_start_codon = start_codon_energies[cds[0:3]]

    """
    Energy of mRNA folding
    """
    mRNA = leader + cds
    structure_mRNA, dG_mRNA = RNAfold.fold(mRNA, d=2)

    """
    Energy of mRNA:rRNA hybridization & folding
    """
    mRNA_rRNA_list = RNAfold.subopt([leader, rRNA], e=energy_cutoff, d=2)
    calc_dG_spacing = (mRNA_rRNA_list, optimal_spacing, dG_push, dG_pull)
    dG_spacing, bps_mRNA_rRNA = RBS_predict_calc_dG_spacing(*calc_dG_spacing)
    if not bps_mRNA_rRNA:
        return 1e12


    # The rRNA-binding site is between the nucleotides at positions most_5p_mRNA
    # and most_3p_mRNA

    # Now, fold the pre-sequence, rRNA-binding-sequence and post-sequence
    # separately. Take their base pairings and combine them together. Calculate
    # the total energy. For secondary structures, this splitting operation is
    # allowed.

    # We postulate that not all of the post-sequence can form secondary
    # structures. Once the 30S complex binds to the mRNA, it prevents the
    # formation of secondary structures that are mutually exclusive with
    # ribosome binding. We define self.footprint to be the length of the 30S
    # complex footprint. Here, we assume that the entire mRNA sequence
    # downstream of the 16S rRNA binding site can not form secondary structures.

    #Add rRNA:mRNA basepairings
    bp_mRNA, bp_rRNA = zip(*bps_mRNA_rRNA)
    total_bp_x = list(bp_mRNA[:])
    total_bp_y = [y+len(leader)+len(cds) for y in bp_rRNA]

    #Calculate pre_mRNA sequence folding
    pre_mRNA = mRNA[0:min(bp_mRNA)]
    if pre_mRNA:
        pre_mRNA_structure, pre_mRNA_dG = RNAfold.fold(pre_mRNA, d=2)
        strands, bp_x, bp_y = brackets_to_basepairing(pre_mRNA_structure)
        total_bp_x.extend(bp_x)
        total_bp_y.extend(bp_y)

    post_mRNA = cds[footprint:]
    if post_mRNA:
        post_mRNA_structure, post_mRNA_dG = RNAfold.fold(post_mRNA, d=2)
        strands, bp_x, bp_y = brackets_to_basepairing(post_mRNA_structure)
        #TODO: TEST
        offset = len(leader) + footprint
        total_bp_x.extend([x + offset for x in bp_x])
        total_bp_y.extend([y + offset for y in bp_y])

    bpseq = "&".join([mRNA, rRNA])
    brackets = basepairing_to_brackets(total_bp_x, total_bp_y, bpseq)
    dG_mRNA_rRNA = RNAfold.energy(bpseq, brackets, d=2)
    # print("me", brackets)

    #dG with spacing
    dG_mRNA_rRNA_w_spacing = dG_spacing + dG_mRNA_rRNA

    """
    Calculate standby site energy
    """

    # To calculate the mfe structure while disallowing base pairing at the
    # standby site, we split the folded mRNA sequence into three parts:
    # (i) a pre-sequence (before the standby site) that can fold;
    # (ii) the standby site, which can not fold;
    # (iii) the 16S rRNA binding site and downstream sequence, which has been
    # previously folded.
    pre_standby = leader[0:max(0, min(bp_mRNA) - len_standby_site)]
    if not pre_standby:
        strands, standby_bp_x, standby_bp_y = [], [], []
    else:
        #Fold it
        brackets, dG = RNAfold.fold(pre_standby, d=2)
        strands, standby_bp_x, standby_bp_y = brackets_to_basepairing(brackets)

    #Append all basepairing where x > the farthest 5p NT in mRNA:rRNA.
    for x, y in zip(total_bp_x, total_bp_y):
        if x >= min(bp_mRNA):
            standby_bp_x.append(x)
            standby_bp_y.append(y)

    brackets = basepairing_to_brackets(standby_bp_x, standby_bp_y, bpseq)
    dG_standby_after = RNAfold.energy(bpseq, brackets, d=2)
    #print(brackets)

    dG_standby = dG_mRNA_rRNA - dG_standby_after

    #Total energy is mRNA:rRNA + start - rRNA - mRNA - standby_site
    dG_total = dG_mRNA_rRNA_w_spacing + dG_start_codon - dG_mRNA - dG_standby

    if verbose:
        print("{:>15}{:>15}{:>15}{:>15}{:>15}".format("dG(total)", "dG(mRNA:rRNA)", "dG(mRNA)", "dG(spacing)", "dG(standby)"))
        print("{:>15.4f}{:>15.4f}{:>15.4f}{:>15.4f}{:>15.4f}".format(dG_total, dG_mRNA_rRNA, dG_mRNA, dG_spacing, dG_standby))

    return dG_total


def RBS_predict_calc_dG_spacing(mRNA_rRNA_list, optimal_spacing, dG_push, dG_pull):
    """Calculate dG spacing.

    Returns dG_spacing and a list of pairs binding mRNA:rRNA:

        # List from RNAsubopt
        >>> l = [('..................(((.((((.........&.)))).)))', -5.9),
        ... ('.....................(.((((((......&..)))))))', -4.8),
        ... ('.....................(((.((((......&..)))))))', -5.2),
        ... ('.........................((((......&..))))...', -6.3),
        ... ('.........................((((.(....&.)))))...', -6.2),
        ... ('......................((((.........&.))))....', -6.1)]

        >>> RBS_predict_calc_dG_spacing(l, 4, (12.2, 2.5, 2.0, 3.0), (0.048, 0.24, 0.0))
        (0.0, [(25, 5), (26, 4), (27, 3), (28, 2)])

    """
    dG_mRNA_rRNA_spc = 1e12
    dG_mRNA_rRNA     = 1e12
    bps_mRNA_rRNA    = []
    for brackets, dG in mRNA_rRNA_list:
        #Calculate spacing penalty
        spacing, bps = RBS_predict_calc_spacing(brackets)

        #See article for help with these equations
        if spacing < optimal_spacing:
            ds = spacing - optimal_spacing
            c1, c2, c3, c4 = dG_push
            dG_spacing_penalty = c1 / (1.0 + math.exp(c2 * (ds + c3)))**c4
        else:
            ds = spacing - optimal_spacing
            c1, c2, c3 = dG_pull
            dG_spacing_penalty = c1 * ds**2 + c2 * ds + c3

        #Proceed with smallest penalty
        if (dG_spacing_penalty + dG) < (dG_mRNA_rRNA_spc + dG_mRNA_rRNA):
            dG_mRNA_rRNA_spc = dG_spacing_penalty
            dG_mRNA_rRNA = dG
            bps_mRNA_rRNA = bps
            # print("mu", brackets, dG_mRNA_rRNA_spc, dG_mRNA_rRNA, spacing)

    return dG_mRNA_rRNA_spc, bps_mRNA_rRNA


def RBS_predict_calc_spacing(brackets):
    """Calculate spacing (NTs from rRNA to start codon).

    Returns the spacing plus a list of basepairings involved in mRNA:rRNA:

        >>> RBS_predict_calc_spacing("....(((............)))((((((((.....&)).))))))")
        (5, [(22, 8), (23, 7), (24, 6), (25, 5), (26, 4), (27, 3), (28, 1), (29, 0)])

    Consider the following sequence with folded rRNA:

                             AAU CCUCCA
                             ||| ||||
        GATACACCAAAGAGAAAATAATGAGGGAGCGAAGG

    This folding has the bracket notation:

        .....................(((.((((......&..)))))))

    Aligning the bracket notation gives:

                                         oo<<<<
                                 ((( ((((..
            .....................(((.((((......

    Thus spacing is 4 (where the <'s are). To calculate this spacing, the
    nt pairs that are mRNA:rRNA are isolated and the "overhang" (two o's) are
    calculated as the rRNA NT closest to the 5' end, while the spacing minus
    overhang is calculated as the mRNA NT closest to the 3' end.

    """
    #Convert to basepairing
    strands, bp_x, bp_y = brackets_to_basepairing(brackets)
    len_leader = strands[0]

    #Pick only pairs binding mRNA:rRNA
    bps = [(x, y-len_leader) for x,y in zip(bp_x,bp_y) if y >= len_leader]
    if not bps: return 1e12, None

    #Unzip basepairing
    bp_mRNA, bp_rRNA = zip(*bps)
    spacing = len_leader - max(bp_mRNA) - min(bp_rRNA) - 1

    return spacing, bps


def RBS_predict_calc_spacing_old(brackets):
    """Old version from RBScalc, possibly broken.

    Returns the spacing plus a list of basepairings involved in mRNA:rRNA:

        >>> RBS_predict_calc_spacing_old("....(((............)))((((((((.....&)).))))))")
        (4, [(22, 8), (23, 7), (24, 6), (25, 5), (26, 4), (27, 3), (28, 1), (29, 0)])

    """
    #Convert to basepairing
    strands, bp_x, bp_y = brackets_to_basepairing(brackets)
    len_leader = strands[0]

    bps = [(x, y-len_leader) for x,y in zip(bp_x,bp_y) if y >= len_leader]
    if not bp_y or not bps: return 1e12, None

    rRNA_nt = max(bp_y)
    rRNA_pos = bp_y.index(rRNA_nt)
    if bp_x[rRNA_pos] < len_leader:
        farthest_3_prime_rRNA = rRNA_nt - len_leader

        mRNA_nt = bp_x[rRNA_pos]
        distance_to_start = len_leader - mRNA_nt

    spacing = distance_to_start - farthest_3_prime_rRNA -1


    return spacing, bps


if __name__ == "__main__":
    import doctest
    doctest.testmod()
