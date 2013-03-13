#!/usr/bin/env python

"""
Mutation generation routines related to translation initiation
"""

from __future__ import print_function

import random
import math

from mutation_tools import find_mutation_box
from mutation_tools import compare_seqs
from RBS.RBS_Calculator import RBS_Calculator
from RBS.RBS_Calculator import NoRBSError
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

"""
RBS Monte Carlo functions
"""

#Constants
RT_EFF = 2.222
LOG_K =  7.824
K = 2500.0

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

    #Default possible moves
    moves = {"A": ["T", "G", "C"],
             "T": ["A", "G", "C"],
             "G": ["A", "T", "C"],
             "C": ["A", "T", "G"]}

    start_codons = ["ATG", "GTG", "TTG"]

    #Maximise/minimise translation rates
    high = -17. # ~5M
    low = 17. # ~1

    #Verbose mode
    verbose = False
    live_plot = False

    #Live plot sillyness
    plot_dG = list()
    plot_n_mut = list()
    plot_temp = list()

    def __init__(self, gene, target):
        """Set up Monte Carlo simulations

        Supply target in dG, or "high" to maximise, "low" to minimise

        """
        self.original_leader = str(gene.leader).upper()
        self.leader_len = len(gene.leader)
        self.leader = list(self.original_leader)
        self.cds = str(gene.cds).upper()
        self.gene = gene

        #To keep track of mutations
        self.mutations = set()

        if target == "high":
            self.target = self.high
        elif target == "low":
            self.target = self.low
        else:
            try:
                self.target = float(target)
            except ValueError:
                raise ValueError("Could not parse: {} as target dG (use only lower-case for high/low)".format(target))

        self.original_dG = RBSPredict(self.original_leader, self.cds)["dG"]
        self.dG = self.original_dG
        self.energy = abs(self.original_dG - self.target)

        #Status
        self.target_reached = False
        self.last_P = 0.0

    def run(self, max_mutations=7, iterations=1, max_rounds=10000):
        """Run simulations"""
        self.max_mutations = max_mutations

        #Start by mutating away start codons
        self.eliminate_start_codons()

        start_mutations = len(self.mutations) + 1
        if start_mutations >= max_mutations:
            raise Exception("start_mutations ({}) is larger than max_mutations ({}) due to start codons.".format(start_mutations, max_mutations))

        #Total iterations
        self.total_its = 0
        rejects = 0

        #Main loop
        for n_muts in range(start_mutations, max_mutations+1):
            #Temperature control
            start_temp = self.start_temp(n_muts)
            temp_step = (start_temp - self.end_temp) / self.min_mutation_period
            current_temp = start_temp

            for i in range(self.max_mutation_period):
                do_move = True
                while do_move:
                    mutations = self.mutations.copy()
                    candidate = self.leader[:]

                    #Pick a random position to move
                    mp = random.randint(0, self.leader_len-1)
                    candidate[mp] = random.choice(self.moves[self.leader[mp]])

                    #Reset mutation list if mutations mutate to original sequence
                    if mp in mutations:
                        if candidate[mp] == self.original_leader[mp]:
                            mutations.remove(mp)
                    else:
                        mutations.add(mp)

                    #Revert a position
                    if len(mutations) > n_muts:
                        rev = random.choice(list(mutations))
                        candidate[rev] = self.original_leader[rev]
                        mutations.remove(rev)

                    #Check for illegal moves
                    do_move = not self.check_leader(candidate)

                #Calculate new total free energy
                new_dG = RBSPredict("".join(candidate), self.cds)["dG"]

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
                    P = math.exp((self.dG - new_dG)/current_temp)
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

                if self.live_plot: 
                    self.plot_dG.append(self.dG)
                    self.plot_n_mut.append(n_muts)
                    self.plot_temp.append(current_temp)
                    if self.total_its % 10 == 0:
                        self.update_live_plot()

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

    def update_live_plot(self):
        """Silly function"""
        import matplotlib.pyplot as plt
        plt.ion()
        plt.clf()
        plt.plot(self.plot_dG, "-")
        plt.plot(self.plot_temp, "-g")
        plt.ylabel("dG / temp")
        plt.axis([0, len(self.plot_dG), self.target, 4])
        plt.twinx()
        plt.ylabel("Mutations")
        plt.plot(self.plot_n_mut, "-r")
        plt.axis([0, len(self.plot_dG), 0, self.max_mutations+1])
        plt.draw()
        plt.ioff()


def AU_to_dG(AU):
    """Convert a dG value to expression level (AU)"""
    return RT_EFF * (LOG_K - math.log(AU))


def dG_to_AU(dG):
    """Convert an expression level to dG"""
    return K * math.exp(-dG / RT_EFF)

"""
RBS Calculator
"""

def RBSPredict(leader, cds):
    """TODO: Create a more solid RBS calculator"""
    start_range = [len(leader), len(leader)]
    mRNA = leader.upper() + cds.upper()
    tmp_calc = RBS_Calculator(mRNA, start_range, "")
    try:
        tmp_calc.calc_dG()
    except NoRBSError:
        tmp_calc.dG_total_list.append(1e20)

    if not tmp_calc.dG_total_list:
        tmp_calc.dG_total_list.append(1e20)

    expr_lvl = tmp_calc.calc_expression_level(tmp_calc.dG_total_list[0])
    return {"expr_lvl": expr_lvl, "dG": tmp_calc.dG_total_list[0]}