#!/usr/bin/env python

"""
Mutation generation routines related to translation initiation
"""

from __future__ import print_function
from __future__ import absolute_import

import math
import random
import logging
import itertools

from mage_tool.helpers import valid_rna
from mage_tool.ViennaRNA import ViennaRNA
from mage_tool.ViennaRNA import brackets_to_basepairing
from mage_tool.ViennaRNA import basepairing_to_brackets
from mage_tool.operations import BaseOperation


#Define a log
log = logging.getLogger("Translation")
log.addHandler(logging.NullHandler())


def replace_start_codon(gene, start_codon="ATG"):
    """Replace start codon.

    Start codon will be replaced by start_codon, which must be a string with
    length 3. (Default 'ATG')

        >>> from mage_tool.oligo_design import Gene
        >>> gene = Gene("ficX", 100, 1, "ATTATATAGACT")
        >>> replace_start_codon(gene)
        Mutation: [T->g] at pos 102

        >>> gene.strand = -1
        >>> replace_start_codon(gene, 'GGG')
        Mutation: [AAT->ccc] at pos 97

    If start_codon is equal to gene start codon, this function returns False:

        >>> replace_start_codon(gene, "ATT")
        False

    """
    if str(gene) == "genome":
        log.error("Cannot use replace_start_codon on genome")
        return None

    if str(gene.cds[0:3]).upper() == str(start_codon).upper():
        return False
    if len(start_codon) != 3:
        raise ValueError("start_codon must be 3 nucleotides long. Recieved: "
                         "{} ({})".format(start_codon, len(start_codon)))

    seq = gene.cds.copy()
    for i, nt in enumerate(start_codon):
        seq.mutate(nt, i, in_place=True)

    return gene.do_mutation(seq.get_mutation())


class StartCodonOptimal(BaseOperation):

    """
    ``start_codon_optimal``: Mutates a start codon to the optimal start codon.

    Tries to mutate the start codon to the optimal start codon defined in the
    strain config file (Usually ``ATG``\ ).

    Options and default values:

    - None

    """

    default_options = {}
    required = ()
    genome_allowed = False
    op_str = "start_codon_optimal"

    def run(self):
        mut = replace_start_codon(self.gene, self.config["start_codons"][0])
        if not mut:
            log.debug(str(self) + " Not mutating, optimal start codon found.")
            return []

        code = "OptStart{}".format(len(mut.before))
        return [(mut, code, str(self), [])]


def translational_KO(gene, stop_codons=["TAG", "TAA", "TGA"],
                     KO_mutations=None, KO_frame=None):
    """Knock out a gene with stop-codons with the least possible mutations.

    At least one of each stop_codon will be used.

    KO_mutations is how many subsequent codons should be mutated. Default is
    one per stop codon. If a lower amount of codons is supplied, one per stop
    codon will still be used.  KO_frame is how many codons the method will look
    for stop codon mutations. Default is half the length of the gene.

    """
    if not KO_mutations or KO_mutations < len(stop_codons):
        KO_mutations = len(stop_codons)

    #Do not mutate the start codon
    start_offset = 3
    #Calculate a KO frame if not provided
    if not KO_frame:
        KO_frame = int((len(gene.cds)/3)/2)

    #Sequence to be mutated
    KO = gene.cds[start_offset:KO_frame*3]
    #List of stop codons and how many mutations each codon require to be
    #mutated to this stop codon
    codon_muts = list()
    #Iterate each codon
    for i in range(0, len(KO), 3):
        stop_muts = list()
        parent = KO[i:i+3]
        #Calculate how many mutations each stop codon require
        for child in stop_codons:
            needed_muts = 0
            for p,m in zip(parent,child):
                if not p == m:
                    needed_muts += 1
            stop_muts.append([needed_muts, child])
        codon_muts.append(stop_muts)

    stop_codons = frozenset(stop_codons)
    totalmuts_codon = list()

    #Iterate each codon-position
    for i in range(len(codon_muts)-KO_mutations+1):
        #Find all combinations of stop codons
        combinations = list()
        for lists in itertools.product(*codon_muts[i:i+KO_mutations]):
            combinations.append(list(lists))
        #Sort combinations by total number of required mutations
        combinations.sort(key=lambda x: sum(e[0] for e in x))
        #Iterate each combination
        for comb in combinations:
            stops = frozenset([e[1] for e in comb])
            #Find first combination with every stop codon included
            #TODO: Get smaller window lengths by checking length of stops set
            if stops == stop_codons:
                #Hom many mutations this combination require
                muts = sum([e[0] for e in comb])
                #New sequence
                after = "".join([e[1] for e in comb])
                pos = i+1
                totalmuts_codon.append([muts, after, pos])
                break

    #Sort all best combinations by number of required mutations
    totalmuts_codon.sort(key=lambda x:x[0])
    #Choose best candidate
    n_muts, muts, offset = totalmuts_codon[0]

    #Mutate sequence to stop codons
    seq = gene.cds.copy()
    for i, nt in enumerate(muts):
        seq.mutate(nt, i + offset*3, in_place=True)

    mutation = gene.do_mutation(seq.get_mutation())
    mutation._codon_offset = offset + 1

    return mutation


class TranslationalKnockout(BaseOperation):
    """
    ``translational_knockout``: Gene knock-out by premature stop-codon.

    Tries to knock out a gene by introducing a number of early stop-codons in
    the CDS with the least amount of mutations possible.

    Options and default values:

    - ``ko_frame=0``: Number of codons that are applicable to be mutated.  E.g.
      a value of 10 means the operation will try to mutate stop codons into the
      CDS within 10 codons of the start codon. The default of 0 is within one
      half of the length of the CDS.
    - ``ko_mutations=3``: number of stop codons to introduce. Default (and
      minimum) is the number of different stop codons available in the genome
      configuration file (normally 3).

    """

    default_options = {"ko_frame": (int, 0),
                       "ko_mutations": (int, 3)}
    required = ()
    genome_allowed = False
    op_str = "translational_knockout"

    def run(self):
        stop_codons = self.config["stop_codons"]
        ko_mutations = self.options["ko_mutations"]
        ko_frame = self.options["ko_frame"]
        mut = translational_KO(self.gene, stop_codons, ko_mutations, ko_frame)
        code = "TransKO{}".format(mut._codon_offset)
        return [(mut, code, str(self), [0, mut._codon_offset])]


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
    if n < 1:
        raise ValueError("'n' must be positive integer.")

    targetAU = target
    targetdG = AU_to_dG(target)

    pre_lib = RBS_Monte_Carlo(gene, targetdG, max_mutations, True)
    pre_lib.sort(key=lambda seq: seq.delta)

    orgdG = RBS_predict(gene.leader, gene.cds)
    orgAU = dG_to_AU(orgdG)
    orgAU = 1e-3 if not orgAU else orgAU
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


def RBS_Monte_Carlo(gene, target, maxmuts=10, collect_library=False, **kwargs):
    """Run a Monte Carlo simulation on a Gene to reach target (dG, kcal/mol).

    target is supplied in dG [kcal/mol]. Use maxmuts to set an upper limit on
    how many mutations are allowed. The simulation will attempt to reach target
    with as few as possible mutations.

    kwargs::
        verbose=False
        calc_start_temp=lambda n:(0.6/n)*3
        end_temp=0.01
        mut_period=1500
        tolerance=0.25
        opti_rounds=25
        RNA_fold=ViennaRNA()

    Returns a list of Sequence objects.

    """

    #Start temp is a function of muts
    calc_start_temp = kwargs.get("calc_start_temp", lambda n: (0.6 / n) * 3)
    end_temp        = kwargs.get("end_temp", 0.01)

    #Rounds per simulation
    mut_period = kwargs.get("mut_period", 1500)

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

    new_dG = P = total_its = current_temp = 0
    def status(msg):
        if verbose:
            fmts = (msg, total_its, best_dG, new_dG, n, P, current_temp)
            print("{:<18} {:>2} {:>5.2f} {:>5.2f} {:>2} {:4.2f} {:3.2f}".format(*fmts))

    #Number of rejects
    rejects = 0

    #Main loop; n = number of mutations
    for n in range(start_mutations, maxmuts + 1):
        #Temperature control
        start_temp = calc_start_temp(n)
        temp_step = (start_temp - end_temp) / mut_period
        current_temp = start_temp

        #TODO: do ALL mutations
        for i in xrange(mut_period):
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
                        status("Cond_accepted")
                    else:
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

            #Decrease temp.
            current_temp -= temp_step

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
    dG_spacing, bps_mRNA_rRNA, spacing = RBS_predict_calc_dG_spacing(*calc_dG_spacing)
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
        (0.0, [(25, 5), (26, 4), (27, 3), (28, 2)], 4)

    """
    dG_mRNA_rRNA_spc = 1e12
    dG_mRNA_rRNA     = 1e12
    bps_mRNA_rRNA    = []
    al_spacing       = 0
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
            al_spacing = spacing

    return dG_mRNA_rRNA_spc, bps_mRNA_rRNA, al_spacing


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


def rbs_method(s):
    """RBS library method."""
    s = s.lower()
    if s in ["exp", "fuzzy"]:
        return s
    else:
        raise ValueError("Unknown RBS_library method: {}".format(s))


class RBSLibrary(BaseOperation):

    """
    ``RBS_library``: Create a library of different RBS expression levels.

    Options and default values:

    - ``target=5000000`` Target expression level to reach for in AU. If target
      is reached, computation is stopped, and library will be created.  if
      target is not reached within the specified number of mutations, a library
      of expression levels closest to target as possible will be created.
    - ``n=10`` Number of library sequences to create.
    - ``max_mutations=10`` Maximum number of mutations to attempt.
    - ``method=exp`` How to create the library. Two methods are available:
      ``exp`` and ``fuzzy``. ``exp`` creates a library where each new sequence
      is an ``m``-fold improvement over the last. ``m`` can either be supplied
      via the ``m``-parameter, or calculated automatically to create an evenly
      spaced library from wt level to target. The ``exp`` method runs multiple
      Monte Carlo simulations to reach each target, however, it uses
      information from previous runs to more quickly reach subsequent targets.
      ``fuzzy`` tries to replicate the ``exp`` library, only the Monte Carlo
      simulation is only run once, and inbetween are collected along the way.
      This method yields a less precise library, but is quicker. Additionally,
      ``fuzzy`` enables picking out the best sequences below a certain mutation
      count by using the ``m`` parameter.  Fx. using ``m=6``, ``fuzzy`` will
      collect the best possible sequences with a maximum of 1, 2, .. 6
      mutations. It will then try to fill out the rest of the library with
      evenly spaced sequences.
    - ``m=0`` see ``method`` for explanation on ``m``.

    This operation will run an Monte-Carlo simulation in an attempt to reach
    the specified target value within a number of mutations. Lower numbers of
    mutations are tried first and are always prioritised over similar
    expression levels with a higher mutation count.

    """

    default_options = {"target": (float, 5000000.),
                       "n": (int, 10),
                       "max_mutations": (int, 10),
                       "method": (rbs_method, "exp"),
                       "m": (float, 0)}
    required = ()
    genome_allowed = False
    op_str = "RBS_library"
    heavy = True

    def post_init(self):
        #Correct very low target values
        if self.options["target"] < 0.1:
            val = self.options["target"]
            self.options["target"] = 0.1
            self.create_opt_str()
            log.debug("{} adjusting very low target value '{}' to 0.1".format(self, val))

        if str(self.gene.cds[0:3]) not in ('ATG', 'GTG', 'TTG', 'CTG'):
            self.error('Invalid start codon: ' + str(self.gene.cds[0:3]))

    def run(self):
        method = self.options["method"]
        target = self.options["target"]
        n = self.options["n"]
        m = self.options["m"]
        max_mutations = self.options["max_mutations"]
        if method == "exp":
            muts = RBS_library(self.gene, target, n, max_mutations, m)
        elif method == "fuzzy":
            muts = RBS_library_fuzzy(self.gene, target, n, max_mutations, m)

        if not muts:
            return None

        muts_out = list()
        for i, m in enumerate(muts):
            code = "RBSlib{}_{:.1f}/{:.1f}({})".format(i, m._AU, m._orgAU, m._n)
            l_op = str(self) + " {:.3f} (wt: {:.2f})".format(m._AU, m._orgAU)
            muts_out.append((m, code, l_op, [m._orgAU, m._AU]))

        return muts_out


OPERATIONS = dict()

def register_operation(op):
    OPERATIONS[op.op_str] = op

register_operation(StartCodonOptimal)
register_operation(TranslationalKnockout)
register_operation(RBSLibrary)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
