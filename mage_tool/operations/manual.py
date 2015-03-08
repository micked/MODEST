#!/usr/bin/env python

"""
Custom mutations
"""

from __future__ import print_function
from __future__ import absolute_import

import re
import logging
from itertools import product

from mage_tool.helpers import dgn_to_nts
from mage_tool.oligo_design import Mutation
from mage_tool.operations import BaseOperation
from mage_tool.helpers import default_dgn_table
from mage_tool.helpers import default_codon_usage
from mage_tool.helpers import default_codon_table
from mage_tool.helpers import reverse_complement


#Define a log
log = logging.getLogger("manual")
log.addHandler(logging.NullHandler())


class FindMutation(BaseOperation):

    """
    ``find_mutation``: Automatically find a mutation in a sequence.

    Options:

    - ``mutation=UPSTREAM[BEFORE->AFTER]DOWNSTREAM``

    Example::

        rpoD find_mutation mutation=GAGCAA[AAC->TAG]CCG

    Finds  GAGCAA\ *AAC*\ CCG in the sequence an changes it to
    GAGCAA\ *TAG*\ CCG by creating the mutation AAC->TAG at position 10.

    """

    default_options = {'mutation': (str, '')}
    required = ('mutation')
    genome_allowed = True
    op_str = 'find_mutation'

    def post_init(self):
        check = re.match(r"^(\w*)\[(\w*)->(\w*)\](\w*)$", self.options['mutation'])
        if not check:
            self.error('Invalid find_mutation format. Expected UPSTREAM[BEFORE->AFTER]DOWNSTREAM, got \'{}\'.'.format(self.options['mutations']))
            return

        self.mut = check.groups()

    def run(self):
        rev_com = False
        upstream, before, after, downstream = self.mut
        seq = str(self.gene.leader) if self.gene.leader else ''
        seq = (seq + str(self.gene.cds)).upper()

        leader_len = len(self.gene.leader) if self.gene.leader else 0
        offset = len(upstream) - leader_len

        #Original sequence
        orig = (upstream + before + downstream).upper()
        #Find occurences
        count = seq.count(orig)
        if count == 0:
            #Try reverse complement.
            upstream, after, before, downstream = [reverse_complement(e) for e in self.mut][::-1]
            orig = (upstream + before + downstream).upper()
            count = seq.count(orig)
            offset = len(upstream)
            rev_com = True
            log.info(str(self) + ' \'{}\' not found in sequence. Trying reverse complement'.format(self.options['mutation']))
            if count == 0:
                log.error(str(self) + ' Mutation not found in sequence.'.format(upstream, before, after, downstream))
                return None
        if count > 1:
            if rev_com:
                log.error(str(self) + ' Ambiguous mutation. \'{}[{}->{}]{}\' found more than once.'.format(upstream, before, after, downstream))
            else:
                log.error(str(self) + ' Ambiguous mutation. \'{}\' found more than once.'.format(self.options['mutation']))
            return None

        pos = seq.find(orig) + offset
        mut = self.gene.do_mutation(Mutation(before, after, pos))
        return [(mut, 'fmut', str(self), [])]


class DNAMutation(BaseOperation):

    """
    ``mutation``: Any kind of nucleotide substitution.

    Options:

    - ``position``: Position of first nucleotide.
    - ``mutation=BEFORE->AFTER``: mutate BEFORE to AFTER in sequence.

    Substitute BEFORE nucleotides in sequence to AFTER at position.

    Example sequence: ``ATGGCTGAA``.

    Mutate ``GCT`` to ``AAA`` at position 4::

        araC mutation mutation=GCT->AAA,position=4
        #Sequence after mutation: ATGaaaGAA.

    Insert ``AAA`` as position 4::

        araC mutation mutation=->AAA,position=4
        #Sequence after mutation: ATGaaaGCTGAA.

    Delete ``GCT`` at position 4::

        araC mutation mutation=GCT->,position=4
        #Sequence after mutation: ATGGAA.

    Substitute ``GCT`` with ``AAATTT``::

        araC mutation mutation=GCT->AAATTT,position=4
        #Sequence after mutation: ATGaaatttGAA.

    """

    default_options = {'mutation': (str, ''),
                       'position': (int, None)}
    required = ('mutation', 'position')
    genome_allowed = True
    op_str  = 'mutation'

    def post_init(self):
        self.mut = None

        mut_str = self.options['mutation']
        if mut_str.count('->') != 1:
            self.error('Invalid mutation format. Use \'BEFORE->AFTER\'')
            return

        before, after = mut_str.split('->')
        pos = self.options['position'] - 1
        found = str(self.gene[pos:pos+len(before)])

        if found != before:
            self.error('Trying to mutate \'{}\', but found \'{}\' in sequence at position {}.'.format(before, found, pos+1))
            return

        self.mut = (before, after, pos)

    def run(self):
        mut = self.gene.do_mutation(Mutation(*self.mut))
        return [(mut, 'mutation', str(self), [])]


class Deletion(BaseOperation):

    """
    ``deletion``: Delete nucleotides.

    Options:

    - ``delete``: What to delete. An integer (number of nucleotides) or a
      string of nucleotides.
    - ``position``: Position of deletion.

    Delete the nucleotides specified in ``delete`` from a DNA sequence.

    Delete AAC from position 13 in thiD::

        thiD deletion delete=AAC,position=13

    Delete 3 nucleotides from position 13 in thiD::

        thiD deletion delete=3,position=13

    Specifying nucleotides to delete enables error-checking, but when deleting
    large stretches, an integer is more handy.

    """

    default_options = {'delete': (str, ''),
                       'position': (int, None)}
    required = ('delete', 'position')
    genome_allowed = True
    op_str = 'deletion'

    def post_init(self):
        pos = self.options['position'] - 1
        self.mut = None

        try:
            #Number of nucleotides to delete
            n_nts = int(self.options['delete'])
            #Nucleotides to delete
            nts = str(self.gene[pos:pos+n_nts])
        except ValueError:
            nts = self.options['delete']
            n_nts = len(nts)

        #Error Checking
        if  pos < 0 or pos + n_nts > len(self.gene.cds):
            self.error('It is not (yet) possible to delete nucleotides beyond gene CDS, use genome deletion instead.')
            return
        found = self.gene[pos:pos+n_nts]
        if nts != str(found):
            self.error("Tried to delete '{}', but found '{}' in sequence.".format(nts, found))
            return

        self.mut = (nts, '', pos)

    def run(self):
        if self.mut is None:
            return None

        #Perform mutation
        mut = self.gene.do_mutation(Mutation(*self.mut))
        return [(mut, 'deletion', str(self), [])]


class Insertion(BaseOperation):

    """
    ``insertion``: Insert nucleotides.

    Options:

    - ``insert``: String of nucleotides to insert.
    - ``position``: Insertion position.

    Examples::

        thiD insertion insert=TTT,position=7
        #Insert TTT as position 7 in thiD.

    Sequence before is ATGAAACGA, and after insertion it is ATGAAAtttCGA.

    """

    default_options = {'insert': (str, ''),
                       'position': (int, None)}
    required = ('insert', 'position')
    genome_allowed = True
    op_str = 'insertion'

    def run(self):
        ins = self.options['insert'].upper()
        pos = self.options['position'] - 1
        mut = self.gene.do_mutation(Mutation('', ins, pos))
        return [(mut, 'insertion', str(self), [])]


def residue_mutation(gene, mutations, codon_table=default_codon_table,
                     dgn_table=default_dgn_table, usage_table=default_codon_usage):
    """Substitue residue.

    Substitutions are given as a list of amino acid mutations:

        >>> from mage_tool.oligo_design import Gene
        >>> gene = Gene("ficX", 100, 1, "ATGGCAACAATAAAAGATGTAGCGAAACGA", "")

        #Mutate K5 to A and V7 to A
        >>> residue_mutation(gene, ["K5A", "V7A"])
        Mutation: [AAAGATGT->gcAGATGc] at pos 112

    Deletions are denoted with a star:

        >>> residue_mutation(gene, ["K5*"])
        Mutation: [AAA->] at pos 112

    Insertions are denoted with a star or `@` and a lower letter index:

        >>> residue_mutation(gene, ["*1aA", "*2aK"])
        Mutation: [GCA->gcgGCAaaa] at pos 103

        >>> residue_mutation(gene, ["@1aA", "@2aK"])
        Mutation: [GCA->gcgGCAaaa] at pos 103

    Stop-codons are denoted by $:

        >>> residue_mutation(gene, ["R10$"])
        Mutation: [C->t] at pos 127

    """
    mutations = sorted(mutations, key=lambda x: x[1:-1])
    seq = gene.cds.copy()
    offset = 0

    for mut in mutations:
        #Find desired residue substitution.
        m = re.match("^([A-Z*@$])(\d+)([a-z]?)([A-Z*$])$", mut)

        if not m:
            #log.debug("Invalid residue mutation: {}.".format(mut))
            raise Exception("Invalid residue mutation: {}.".format(mut))

        old_AA, pos, pos_letter, new_AA = m.groups()
        pos = int(pos)

        dna_pos = (pos - 1) * 3 + offset

        if not pos_letter and codon_table[seq[dna_pos:dna_pos+3]] != old_AA :
            #log.debug("Invalid residue mutation: {}. Old Residue: {} does not match {} found in cds".format(mut, old_AA, codon_table[seq[dna_pos:dna_pos+3]]))
            raise Exception("Invalid residue mutation: {}. Old Residue: {} does not match {} found in cds".format(mut, old_AA, codon_table[seq[dna_pos:dna_pos+3]]))

        #Deletion
        if new_AA == "*":
            for i in range(3):
                seq.delete(dna_pos, in_place=True)
            offset -= 3
        else:
            #New degenerate codon
            new_dgn = dgn_table[new_AA]
            #All possible new codons
            #new_codons = product(*[list(dgn_to_nts[nt]) for nt in new_dgn])
            #new_codons = ["".join(cdn) for cdn in new_codons]
            new_codons = [c for c, r in codon_table.items() if r==new_AA]

            #Insertion, choose codon with max usage
            if old_AA in {'*', '@'}:
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
                if old_AA in {'*', '@'}:
                    seq.insert(nt, dna_pos + 3 + i, in_place=True)
                else:
                    seq.mutate(nt, dna_pos + i, in_place=True)

            if old_AA in {'*', '@'}: offset += 3

    mutation = seq.get_mutation()
    if mutation:
        return gene.do_mutation(mutation)

    raise Exception('Unknown error occured.')
    return None


class ResidueMutation(BaseOperation):

    """
    ``residue_mutation``: Mutating a residue.

    Options:

    - ``mut`` Amino acid substitions. Multiple substitions can be separated by a semicolon (;)

    Examples::

        thiD residue_mutation mut=N5Q
        #substitue N with Q at residue 5 in thiD.

    Amino acid sequence before is MKRINALTIA, and after substitution it is MKRIQALTIA.

    Deletions are denoted with a ``*``::

        thiD residue_mutation mut=N5*

    Insertions are denoted as a an ``@`` followed by an insertion number,
    lower-case suffix letter and then the inserted amino acid::

        thiD residue_mutation mut=@5aA
        #insert alanine in thiD after position 5.

        thiD residue_mutation mut=@5aA;@5bA
        #insert two alanines in thiD after position 5.

    The symbol for stop codons is $::

        thiD residue_mutation mut=N5$
        #substitue N with a stop codon at
        #residue 5 in thiD.

    """

    default_options = {"mutation": (str, None)}
    required = ("mutation",)
    genome_allowed = False
    op_str = "residue_mutation"

    def post_init(self):
        self.muts = self.options['mutation'].split(';')
        for mut in self.muts:
            m = re.match(r'^([ACDEFGHIKLMNPQRSTVWY*$])(\d+)([a-z]?)([ACDEFGHIKLMNPQRSTVWY*$])$', mut)
            if not m:
                self.error('Invalid mutation \'{}\' of format \'A10B\' in mut={}'.format(mut, self.options['mutation']))
            else:
                m = m.groups()
                if not m[2]:
                    gene_pos = (int(m[1])-1)*3
                    try:
                        codon_table = self.config['codon_table']
                    except KeyError:
                        codon_table = default_codon_table
                    expected_codon = codon_table[self.gene.cds[gene_pos:gene_pos+3]]
                    if str(m[0]) != str(expected_codon):
                        self.error('Invalid residue mutation: {}. Old Residue: {} does not match {} found in cds'.format(mut, m[0], expected_codon))

    def run(self):
        cdn_tbl = self.config["codon_table"]
        dgn_tbl = self.config["dgn_table"]
        cdn_usage = self.config["codon_usage"]
        mut = self.muts
        try:
            mut = residue_mutation(self.gene, mut, cdn_tbl, dgn_tbl, cdn_usage)
        except Exception as ex:
            self.error(ex)
            return

        code = "ResMut"
        return [(mut, code, str(self), [])]


OPERATIONS = {}
def register_operation(op):
    """Register an operation with the global OPERATIONS dict."""
    OPERATIONS[op.op_str] = op


register_operation(DNAMutation)
register_operation(ResidueMutation)
register_operation(Deletion)
register_operation(Insertion)
register_operation(FindMutation)

if __name__ == "__main__":
    import doctest
    doctest.testmod()
