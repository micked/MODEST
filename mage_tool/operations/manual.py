#!/usr/bin/env python

"""
Custom mutations
"""

from __future__ import print_function
from __future__ import absolute_import

import re
import math
import random
import logging
from itertools import product

from mage_tool.oligo_design import Mutation
from mage_tool.helpers import degenerate_nucleotides
from mage_tool.helpers import dgn_to_nts
from mage_tool.operations import BaseOperation


#Define a log
log = logging.getLogger("MODEST.man")
log.addHandler(logging.NullHandler())


"""
Default values/tables
"""

default_codon_table = {'CTT': 'L', 'ATG': 'M', 'ACA': 'T', 'ACG': 'T',
                       'ATC': 'I', 'AAC': 'N', 'ATA': 'I', 'AGG': 'R',
                       'CCT': 'P', 'CTC': 'L', 'AGC': 'S', 'AAG': 'K',
                       'AGA': 'R', 'CAT': 'H', 'AAT': 'N', 'ATT': 'I',
                       'CTG': 'L', 'CTA': 'L', 'ACT': 'T', 'CAC': 'H',
                       'AAA': 'K', 'CCG': 'P', 'AGT': 'S', 'CCA': 'P',
                       'CAA': 'Q', 'CCC': 'P', 'TAT': 'Y', 'GGT': 'G',
                       'TGT': 'C', 'CGA': 'R', 'CAG': 'Q', 'CGC': 'R',
                       'GAT': 'D', 'CGG': 'R', 'TTT': 'F', 'TGC': 'C',
                       'GGG': 'G', 'TAG': '$', 'GGA': 'G', 'TAA': '$',
                       'GGC': 'G', 'TAC': 'Y', 'GAG': 'E', 'TCG': 'S',
                       'TTA': 'L', 'GAC': 'D', 'CGT': 'R', 'GAA': 'E',
                       'TGG': 'W', 'GCA': 'A', 'GTA': 'V', 'GCC': 'A',
                       'GTC': 'V', 'GCG': 'A', 'GTG': 'V', 'TTC': 'F',
                       'GTT': 'V', 'GCT': 'A', 'ACC': 'T', 'TGA': '$',
                       'TTG': 'L', 'TCC': 'S', 'TCA': 'S', 'TCT': 'S'}


default_dgn_table = {'A': 'GCN', 'C': 'TGY', 'E': 'GAR', '$': 'TRR',
                     'G': 'GGN', 'F': 'TTY', 'I': 'ATH', 'H': 'CAY',
                     'K': 'AAR', 'M': 'ATG', 'L': 'YTN', 'N': 'AAY',
                     'Q': 'CAR', 'P': 'CCN', 'S': 'WSN', 'R': 'MGN',
                     'T': 'ACN', 'W': 'TGG', 'V': 'GTN', 'Y': 'TAY',
                     'D': 'GAY'}

#E. coli default codon usage
default_codon_usage = {'CTT': ['L', 0.12], 'ATG': ['M', 1.0],
                       'ACA': ['T', 0.13], 'ACG': ['T', 0.24],
                       'ATC': ['I', 0.35], 'ATA': ['I', 0.07],
                       'AGG': ['R', 0.03], 'CCT': ['P', 0.17],
                       'AGC': ['S', 0.33], 'AGA': ['R', 0.02],
                       'ATT': ['I', 0.58], 'CTG': ['L', 0.46],
                       'CTA': ['L', 0.05], 'ACT': ['T', 0.16],
                       'CCG': ['P', 0.55], 'AGT': ['S', 0.14],
                       'CCA': ['P', 0.14], 'CCC': ['P', 0.13],
                       'TAT': ['Y', 0.53], 'GGT': ['G', 0.29],
                       'CGA': ['R', 0.07], 'CGC': ['R', 0.44],
                       'CGG': ['R', 0.07], 'GGG': ['G', 0.12],
                       'TAG': ['$', 0.0], 'GGA': ['G', 0.13],
                       'TAA': ['$', 0.64], 'GGC': ['G', 0.46],
                       'TAC': ['Y', 0.47], 'CGT': ['R', 0.36],
                       'GTA': ['V', 0.17], 'GTC': ['V', 0.18],
                       'GTG': ['V', 0.4], 'GAG': ['E', 0.3],
                       'GTT': ['V', 0.25], 'GAC': ['D', 0.35],
                       'GAA': ['E', 0.7], 'AAG': ['K', 0.27],
                       'AAA': ['K', 0.73], 'AAC': ['N', 0.53],
                       'CTC': ['L', 0.1], 'CAT': ['H', 0.55],
                       'AAT': ['N', 0.47], 'CAC': ['H', 0.45],
                       'CAA': ['Q', 0.3], 'CAG': ['Q', 0.7],
                       'TGT': ['C', 0.42], 'TCT': ['S', 0.11],
                       'GAT': ['D', 0.65], 'TTT': ['F', 0.57],
                       'TGC': ['C', 0.58], 'TGA': ['$', 0.36],
                       'TGG': ['W', 1.0], 'TTC': ['F', 0.43],
                       'TCG': ['S', 0.16], 'TTA': ['L', 0.15],
                       'TTG': ['L', 0.12], 'TCC': ['S', 0.11],
                       'ACC': ['T', 0.47], 'TCA': ['S', 0.15],
                       'GCA': ['A', 0.21], 'GCC': ['A', 0.31],
                       'GCG': ['A', 0.38], 'GCT': ['A', 0.11]}


def gene_mutation(gene, mutation):
    """Do specified custom mutations.

    mutation can be either be specified as a specific mutation at a specific
    location using the ``[before=after].pos`` notation:

        >>> from oligo_design import Gene
        >>> gene = Gene("ficX", 100, 1, "ATGATTATAGCACACTA", "AGCAGCATAC")
        >>> gene_mutation(gene, "[ATT=GCC].4")
        Mutation: [ATT=gcc] at pos 103

    Position is 1-indexed relative to gene, which means negative positions are
    allowed:

        >>> gene_mutation(gene, "[TAC=CCA].-3")
        Mutation: [TAC=cca] at pos 97
        >>> gene_mutation(gene, "[ATACATG=GGG].-4")
        Mutation: [ATACATG=ggg] at pos 96

    In case ``before`` does not match sequence, function returns False:

        >>> gene_mutation(gene, "[GGG=].1")
        False

    Another alternative is the lazy method, which automatically finds the
    mutation position based on the specified sequence. ``mutation`` is given as
    a sequence with an embedded mutation:

        pre_seq[before=after]post_seq

    For example, the mutation ``AAAA[GGG=CC]TTTT`` will mutate the sequence
    ``AAAAGGGTTTT`` to ``AAAACCTTTT``.

        >>> gene_mutation(gene, "ATT[AT=GG]AGC")
        Mutation: [AT=gg] at pos 106
        >>> gene_mutation(gene, "ATA[CATG=AGG]ATT")
        Mutation: [CATG=agg] at pos 99

    pre_seq and post_seq are chosen to make the lookup sequence unique. If more
    than one match is found (or none at all), function returns False.

        >>> gene_mutation(gene, "A[GC=T]A")
        False

    For detailed information on why a mutation fail, enable logging (See
    official documentation on Python logs).

    """
    #Type '1', mutation at specific location
    t1 = re.match(r"^\[(\w*)=(\w*)\]\.(-?\d*)$", mutation)

    #Type '2', lazy mutation
    t2 = re.match(r"^(\w*)\[(\w*)=(\w*)\](\w*)$", mutation)

    if t1:
        m = t1.groups()
        before = m[0]
        after = m[1]
        offset = int(m[2])
        if offset > 1: offset -= 1
        if str(gene[offset, offset+len(m[0])]) != str(before):
            log.debug("gene_mutation: {} not found in {}, found {} instead."
                      "".format(before, gene, gene[offset, offset+len(m[0])]))
            return False
    elif t2:
        m = t2.groups()
        before = m[1]
        after = m[2]
        original_string = m[0]+m[1]+m[3]
        #Find mutation offset on gene.
        offset = len(m[0])
        gene_str = str(gene.leader) + str(gene.cds)
        seq_count = gene_str.count(original_string)
        if not seq_count:
            log.debug("gene_mutation: {} not found in {}."
                      "".format(original_string, gene))
            return False
        elif seq_count > 1:
            log.debug("gene_mutation: {} found more than once in {}."
                      "".format(original_string, gene))
            return False
        else:
            offset += gene_str.find(original_string) - len(gene.leader)
    else:
        return False

    #Do mutation
    #mut = "{}={}".format(before, after)
    mutation = Mutation(before, after, offset)
    #Return gene mutation.
    return gene.do_mutation(mutation)


def dna_mut(s):
    """Custom type to verify a DNA mutation."""
    t1 = re.match(r"^\[(\w*)=(\w*)\]\.(-?\d*)$", s)
    t2 = re.match(r"^(\w*)\[(\w*)=(\w*)\](\w*)$", s)
    if not t1 and t2:
        raise ValueError("Invalid DNA mutation: {}".format(s))
    return s


class DNAMutation(BaseOperation):

    """
    ``dna_mutation``: Allows for a desired mutation.

    Options:

    - ``mut=upstream[before=after]downstream``

    I.e. ``mut=TATCAACGCC[GCTC=A]GCTTTCATGACT`` changes
    TATCAACGCC\ **GCTCG**\ CTTTCATGACT to TATCAACGCC\ **A**\ GCTTTCATGACT.

    """

    default_options = {"mut": (dna_mut, None)}
    required = ("mut",)
    genome_allowed = True
    op_str = "dna_mutation"

    def run(self):
        mut = self.options["mut"]
        mut = gene_mutation(self.gene, mut)
        if not mut:
            log.error(str(self) + " Not mutating, did not find mutation box")
            return None

        code = "DNAMut"
        return [(mut, code, str(self), [])]


def residue_mutation(gene, mutations, codon_table=default_codon_table,
                     dgn_table=default_dgn_table, usage_table=default_codon_usage):
    """Substitue residue.

    Substitutions are given as a list of amino acid mutations:

        >>> from oligo_design import Gene
        >>> gene = Gene("ficX", 100, 1, "ATGGCAACAATAAAAGATGTAGCGAAACGA", "")

        #Mutate K5 to A and V7 to A
        >>> residue_mutation(gene, ["K5A", "V7A"])
        Mutation: [AAAGATGT=gcAGATGc] at pos 112

    Deletions are denoted with a star:

        >>> residue_mutation(gene, ["K5*"])
        Mutation: [AAA=] at pos 112

    Insertions are denoted with a star and a lower letter index:

        >>> residue_mutation(gene, ["*1aA", "*2aK"])
        Mutation: [GCA=gcgGCAaaa] at pos 103

    Stop-codons are denoted by $:

        >>> residue_mutation(gene, ["K10$"])
        Mutation: [C=t] at pos 127

    """
    new_dna_string = str(gene.cds)
    mutations = sorted(mutations, key=lambda x: x[1:-1])
    seq = gene.cds.copy()
    offset = 0

    for mut in mutations:
        #Find desired residue substitution.
        m = re.match("^([A-Z*$])(\d+)([a-z]?)([A-Z*$])$", mut)

        if not m:
            log.debug("Invalid residue mutation: {}.".format(mut))
            return None

        old_AA, pos, pos_letter, new_AA = m.groups()
        pos = int(pos)

        dna_pos = (pos - 1) * 3 + offset
        #Deletion
        if new_AA == "*":
            for i in range(3):
                seq.delete(dna_pos, in_place=True)
            offset -= 3
        else:
            #New degenerate codon
            new_dgn = dgn_table[new_AA]
            #All possible new codons
            new_codons = product(*[list(dgn_to_nts[nt]) for nt in new_dgn])
            new_codons = ["".join(cdn) for cdn in new_codons]
            #Insertion, choose codon with max usage
            if old_AA == "*":
                usage, new_codon = max([(usage_table[cdn][1], cdn) for cdn in new_codons])
            #Prioritise number of muts, then usage
            else:
                min_muts = 3
                new_codon = new_codons[0]
                for cdn in new_codons:
                    muts = len([1 for i in range(3) if cdn[i] != seq[dna_pos+offset+i]])
                    #True if we have fever muts
                    if muts < min_muts:
                        min_muts = muts
                        new_codon = cdn
                    #True if number muts is equal, but usage is better
                    elif muts == min_muts and usage_table[cdn] > usage_table[new_codon]:
                        new_codon = cdn

            #Do mutation
            for i, nt in enumerate(new_codon):
                if old_AA == "*":
                    seq.insert(nt, dna_pos + 3 + i, in_place=True)
                else:
                    seq.mutate(nt, dna_pos + i, in_place=True)

            if old_AA == "*": offset += 3

    mutation = seq.get_mutation()
    if mutation:
        return gene.do_mutation(mutation)
    return None


def residue_mutlist(s):
    """Several residue mutations."""
    muts = s.split(";")
    for mut in muts:
        m = re.match("^([A-Z*$])(\d+)([a-z]?)([A-Z*$])$", mut)
        if not m:
            raise ValueError("Invalid mutation: {}".format(mut))
    return muts


class ResidueMutation(BaseOperation):

    """
    ``residue_mutation``: Mutating a residue.

    Options:

    - ``mut=``

    I.e. ``mut=`` changes

    """

    default_options = {"mut": (residue_mutlist, None)}
    required = ("mut",)
    genome_allowed = False
    op_str = "residue_mutation"

    def run(self):
        cdn_tbl = self.config["codon_table"]
        dgn_tbl = self.config["dgn_table"]
        cdn_usage = self.config["codon_usage"]
        mut = self.options["mut"]
        mut = residue_mutation(self.gene, mut, cdn_tbl, dgn_tbl, cdn_usage)
        if not mut:
            log.error(str(self) + " Not mutating.")
            return None

        code = "ResMut"
        return [(mut, code, str(self), [])]


OPERATIONS = {}
def register_operation(op):
    """Register an operation with the global OPERATIONS dict."""
    OPERATIONS[op.op_str] = op


register_operation(DNAMutation)
register_operation(ResidueMutation)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
