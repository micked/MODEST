#!/usr/bin/env python

"""
Mutation generation routines related to translation initiation
"""

from __future__ import print_function

import random
import math

from mutation_tools import find_mutation_box
from mutation_tools import compare_seqs
from oligo_design import Mutation
from helpers import degenerate_nucleotides
from helpers import valid_rna
from ViennaRNA import ViennaRNA
from ViennaRNA import brackets_to_basepairing
from ViennaRNA import basepairing_to_brackets


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
                    #print(start_offset, i, parent, child, mutation, gene.cds[0:50])
                    mutation.pos += i + start_offset
                    #Apply mutation
                    new_mut = gene.do_mutation(mutation)
                    new_mut._codon_offset = (i+start_offset)/3
                    return new_mut


def RBS_single_mutation(gene, maximise=True, insert=False, delete=False, top=3):
    """Maxi/minimise RBS binding with single mutations (Bruteforce approach)"""
    mutation_list = list()
    
    old_value = RBS_predict(gene.leader, gene.cds)
    mutation_list.append((old_value, 0, ""))

    mutations = {"A": ["T", "G", "C"],
                 "T": ["A", "G", "C"],
                 "G": ["A", "T", "C"],
                 "C": ["A", "T", "G"]}

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
            # new_seq = new_leader + str(gene.cds)
            # RBS_calc = RBS_Calculator(new_seq, [34+len(m), 34+len(m)], "")
            # RBS_calc.calc_dG()

            tmp_mut = "{}={}".format(n, m.lower())
            # if not len(RBS_calc.dG_total_list):
            #     continue
            new_value = RBS_predict(new_leader, gene.cds)
            mutation_list.append((new_value, i, tmp_mut))

    muts = list()
    for m in sorted(mutation_list, key=lambda x: x[0], reverse=maximise)[0:top]:
        if not m[2]:
            break
        mut = Mutation("eq", m[2], -len(leader)+m[1])
        mut._adjustment = m[0] / old_value
        muts.append(gene.do_mutation(mut))

    return muts


def generate_RBS_library(gene, target, n, max_mutations, passes):
    """Generate an RBS library

    target          -- Target translation initiation rate (AU).
    n               -- Number of mutations (oligos) to output.
    max_mutations   -- Maximum number of mutations to perform.
    passes          -- Number of passes to perform to increase resolution.

    """
    #Convert target to dG
    target = AU_to_dG(target)

    MC = RBSMonteCarlo(gene, target)
    #MC.verbose = True
    lib = MC.create_library(n, max_mutations, passes)

    if not lib:
        return False

    org_leader = gene.leader
    org_expr = dG_to_AU(MC.original_dG)

    muts_out = list()
    for expr, muts, seq, prox in lib:
        muts = sorted(muts)
        before = org_leader[muts[0]:muts[-1]+1]
        after = list()
        for i in range(muts[0], muts[-1]+1):
            if i in muts:
                after.append(seq[i].lower())
            else:
                after.append(org_leader[i].upper())
        after = "".join(after)

        mut = "{}={}".format(before, after)
        mut = Mutation("eq", mut, -len(org_leader)+muts[0])
        mut._expr = expr
        mut._org_expr = org_expr
        mut._n_muts = len(muts)
        muts_out.append(gene.do_mutation(mut))

    return muts_out

"""
RBS Monte Carlo functions
"""

#Constants
RBS_RT_EFF = 2.222
RBS_LOG_K =  7.824
RBS_K = 2500.0

class RBSMonteCarlo:
    """Monte Carlo simulations

    This class works by ..

    Two periods ..

    """
    #start_temp gets colder when more mutations are introduced
    start_temp = lambda self, n_muts: (0.6 / n_muts) * 3
    end_temp = 0.01
    min_mutation_period = 1000
    max_mutation_period = 4000
    #Rejects before bumping to an extra mutation
    stationary = 500
    #Tolerance to target (dG)
    tol = 0.25

    #Maximum attempted moves before restarting simulation
    max_attempted_moves = 100

    #Keeping track of mutations
    mutations = set()

    # start_codons = ["ATG", "GTG", "TTG"]
    start_codons = ["ATG"]
    prevent_start_codons = False

    #Maximise/minimise translation rates
    high = -17. # ~5M AU
    low = 32. # ~0.001 AU

    mutation_lib = list()

    #Verbose mode
    verbose = False

    total_its = 0

    def __init__(self, gene, target):
        """Set up Monte Carlo simulations

        Supply target in dG, or "high" to maximise, "low" to minimise

        """
        self.RNAfold = ViennaRNA()

        self.original_leader = str(gene.leader).upper()
        self.leader_len = len(gene.leader)
        self.cds = str(gene.cds).upper()
        self.gene = gene

        if target == "high":
            self.target = self.high
        elif target == "low":
            self.target = self.low
        else:
            try:
                self.target = float(target)
            except ValueError:
                raise ValueError("Could not parse: {} as target dG (use only lower-case for high/low)".format(target))

        self.original_dG = RBS_predict(self.original_leader, self.cds, RNAfold=self.RNAfold)

        self.wobble = list()
        #Wobble sequences
        for N in gene.leader_wobble.upper():
            if N in ["A", "T", "G", "C"]:
                self.wobble.append(False)
            else:
                self.wobble.append(degenerate_nucleotides[N])

        if not any(self.wobble):
            raise ValueError("Wobble sequnce [{}] has no possibility for mutation.".format(gene.leader_wobble))

    def create_library(self, n=7, max_mutations=10, passes=1, prioritise_low_count=True):
        """Create an RBS library from start to target with n candidates"""
        self.mutation_lib = list()

        target_expr_lvl = dG_to_AU(self.target)
        #fitness/proximity to target
        f = lambda x: abs(x - target_expr_lvl)
        
        for pas in range(passes):
            self.run(max_mutations=max_mutations, make_library=True)

        mut_lib = [(dG_to_AU(dG), sorted(muts), seq, f(dG_to_AU(dG))) for dG, muts, seq in self.mutation_lib]

        if not mut_lib:
            return False

        selected_library = [min(mut_lib, key=lambda x: x[3])]

        #Loop through each leader with m mutations
        #And pickup the one closest to the target
        if prioritise_low_count:
            for m in range(1, max_mutations+1):
                best_n = [l for l in mut_lib if len(l[1]) == m]
                if best_n:
                    best_n = min(best_n, key=lambda l: l[3])
                    if best_n[2] not in [l[2] for l in selected_library]:
                        selected_library.append(best_n)

        if len(selected_library) >= n:
            return selected_library[0:n]

        n -= len(selected_library) - 1

        best  = min(mut_lib, key=lambda x: x[3])
        worst = max(mut_lib, key=lambda x: x[3])

        stepsize = (best[3] - worst[3]) / n

        #Go through each leader in library and select one closest to target
        for i in range(1, n):
            step = best[3] - (stepsize * i)
            closest = sorted(mut_lib, key=lambda l: abs(l[3] - step))
            for close in closest:
                if close[2] not in [l[2] for l in selected_library]:
                    selected_library.append(close)
                    break

        selected_library.sort(key=lambda l: l[3])

        return selected_library

    def run(self, max_mutations=10, make_library=False):
        """Run simulations"""
        #Initalisation vars
        self.dG = self.original_dG
        self.energy = abs(self.original_dG - self.target)
        self.optimal_energy = self.energy
        self.mutations = set()
        self.leader = list(self.original_leader)

        #Status
        self.target_reached = False
        self.last_P = 0.0

        #Start by mutating away start codons
        if self.prevent_start_codons:
            self.eliminate_start_codons()

        start_mutations = len(self.mutations) + 1
        if start_mutations > max_mutations:
            raise Exception("start_mutations ({}) is larger than max_mutations ({}) due to start codons.".format(start_mutations, max_mutations))

        #Number of rejects
        rejects = 0

        #Main loop
        for n_muts in range(start_mutations, max_mutations+1):
            #Temperature control
            start_temp = self.start_temp(n_muts)
            temp_step = (start_temp - self.end_temp) / self.min_mutation_period
            current_temp = start_temp

            for i in range(self.max_mutation_period):
                #do_move = True
                #attempted_moves = 0
                # while do_move:
                for attempted_move in xrange(self.max_attempted_moves):
                    mutations = self.mutations.copy()
                    candidate = self.leader[:]

                    #Pick a random position to move
                    mp = random.randint(0, self.leader_len-1)
                    while not self.wobble[mp]:
                        mp = random.randint(0, self.leader_len-1)

                    #Pick a new nucleotide to mutate to
                    new_NT = random.choice(self.wobble[mp])
                    while new_NT == self.leader[mp]:
                        new_NT = random.choice(self.wobble[mp])

                    candidate[mp] = new_NT

                    # print("   ", "".join(self.leader))
                    # print(mp, new_NT, "".join(candidate))

                    #Reset mutation list if mutations mutate to original sequence
                    if mp in mutations:
                        if new_NT == self.original_leader[mp]:
                            mutations.remove(mp)
                    else:
                        mutations.add(mp)

                    #Revert a position
                    if len(mutations) > n_muts:
                        rev = random.choice(list(mutations))
                        while rev != mp:
                            rev = random.choice(list(mutations))
                        candidate[rev] = self.original_leader[rev]
                        mutations.remove(rev)

                    #Check for illegal moves
                    if not self.prevent_start_codons or self.check_leader(candidate):
                        break

                else:
                    #Redo this one after a while
                    #Let's see if it gets triggered
                    raise Exception("Too many attempted moves.")

                #Calculate new total free energy
                new_dG = RBS_predict("".join(candidate), self.cds, RNAfold=self.RNAfold)

                #Calculate simulation energy
                new_energy = abs(new_dG - self.target)

                accept_move = False
                if new_energy < self.energy:
                    accept_move = True
                    rejects = 0
                    if self.verbose:
                        self.status("Accepted", new_dG, n_muts)
                elif self.energy == new_energy:
                    rejects += 1
                    if random.random() > 0.8:
                        #Randomly accept
                        accept_move = True
                        if self.verbose:
                            self.status("Randomly_accepted:", new_dG, n_muts)
                    elif self.verbose:
                        self.status("Randomly_rejected:", new_dG, n_muts)
                else:
                    try:
                        P = math.exp((self.energy - new_energy)/current_temp)
                    except OverflowError:
                        if self.verbose:
                            self.status("OverflowError: ", new_dG, n_muts)
                    else:
                        if P > random.random():
                            #Conditionally accept
                            accept_move = True
                            rejects = 0
                            if self.verbose:
                                self.status("Cond_accepted:", new_dG, n_muts, P)
                        else:
                            rejects += 1
                            if self.verbose:
                                self.status("Cond_rejected:", new_dG, n_muts, P)

                #Move is accepted
                if accept_move:
                    self.leader = candidate
                    self.dG = new_dG
                    self.energy = new_energy
                    self.mutations = mutations

                #Library generation
                if make_library:
                    if self.energy < self.optimal_energy:
                        new_addition = (self.dG, list(self.mutations), "".join(self.leader))
                        self.mutation_lib.append(new_addition)
                        self.optimal_energy = self.energy

                #Target reached, STOP
                if new_energy < self.tol:
                    self.target_reached = True
                    return "".join(self.leader)

                #Decrease temp while in "hot period"
                if i < self.min_mutation_period:
                    current_temp -= temp_step
                #Look for stationary phases in "cold period"
                elif rejects >= self.stationary:
                    #Bump number of mutations if we are stationary
                    break

                self.total_its += 1

        return "".join(self.leader)

    def status(self, msg, dG, n, P=None):
        """Print status"""
        if P is None: P = self.last_P
        self.last_P = P
        print("{:<18} {:>2} {:>5.2f} {:>5.2f} {:>2} {:4.2f}".format(msg, self.total_its, self.dG, dG, n, P))

    def check_leader(self, new_seq):
        """Checks the leader sequence for illegal moves.

        Returns False on invalid sequence

        """
        leader_str = "".join(new_seq)
        for stc in self.start_codons:
            if leader_str.find(stc) != -1:
                return False

        return True

    def eliminate_start_codons(self):
        """Start by mutating away other start codons"""
        leader_str = "".join(self.leader)
        for stc in self.start_codons:
            fnd = leader_str.find(stc)
            if fnd != -1:
                self.leader[fnd+2] = random.choice(["A", "T", "C"])
                self.mutations.add(fnd+2)


def AU_to_dG(AU):
    """Convert a dG value to expression level (AU)"""
    if AU == 0:
        return 1e12
    return RBS_RT_EFF * (RBS_LOG_K - math.log(AU))


def dG_to_AU(dG):
    """Convert an expression level to dG"""
    return RBS_K * math.exp(-dG / RBS_RT_EFF)

"""
RBS Calculator
"""

def RBS_predict(leader, cds, nt_cutoff=35, RNAfold=ViennaRNA(), verbose=False):
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
    leader = str(leader).upper().replace("T", "U")
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