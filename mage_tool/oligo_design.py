
"""
Module for designing DNA oligos from mutation objects
"""

from __future__ import print_function

import random
import logging
from copy import deepcopy

from mage_tool import ViennaRNA
import mage_tool.run_control as rc
from mage_tool.helpers import is_inside
from mage_tool.helpers import valid_dna
from mage_tool.helpers import valid_rna
from mage_tool.helpers import valid_dgn
from mage_tool.helpers import dgn_to_nts
from mage_tool.helpers import nts_to_dgn
from mage_tool.helpers import reverse_complement
from mage_tool.helpers import extract_circular
from mage_tool.helpers import make_primer
from mage_tool.helpers import make_rev_primer

#Define a log
log = logging.getLogger("MODEST.oligo")
log.addHandler(logging.NullHandler())

class Oligo:
    """Oligo Object"""

    def __init__(self, mut, gene, project=None, number=0, oligo_len=rc.CONF["oligo_length"]):
        self.mut = mut
        self.gene = gene
        self.oligo_len = oligo_len
        self.pos = mut.pos
        self.oligo = None
        self.optimised = 0
        self.dG_fold = 0
        self.replichore = None
        self.strand = 1
        self.project = project.replace(" ", "_") if project else project
        self.number = number
        self.custom_id = None

        self.barcodes_forward = list()
        self.barcodes_reverse = list()
        self.barcode_ids = list()
        self.code = ""
        self.operation = ""
        self.operation_values = list()

    def set_oligo(self, genome, optimise=True, threshold=-20):
        if optimise and len(self.mut.after) <= self.oligo_len - 2*rc.CONF['min_homology']:
            self.oligo, self.dG_fold, self.optimised = self.optimise_folding(genome, threshold)
        else:
            self.oligo = self.make_oligo(genome, offset=0)
            self.dG_fold = ViennaRNA.mfe(str(self.oligo))

    def make_oligo(self, genome, offset=0):
        """Make oligo from mutation"""
        #Make sure what is being mutated is actually being mutated
        if str(genome[self.mut.pos:self.mut.pos+len(self.mut.before)]) != str(self.mut.before):
            found = genome[self.mut.pos:self.mut.pos+len(self.mut.before)]
            mutpos = self.mut.pos
            len_bef = len(self.mut.before)
            extended = "{}({}){}".format(genome[mutpos-5:mutpos], found,
                                         genome[mutpos+len_bef:mutpos+len_bef+5])
            raise Exception("Trying to mutate {}, but found {} in genome. "
                            "{}".format(self.mut.before, found, extended))

        #In case an insertion is longer than the oligo size
        if len(self.mut.after) > self.oligo_len - 2*rc.CONF['min_homology']:
            post_seq_len = rc.CONF['min_homology']
            pre_seq_len = rc.CONF['min_homology']
        else:
            #Calculate flanking sequence lengths
            post_seq_len = (self.oligo_len - len(self.mut.after))/2
            pre_seq_len = self.oligo_len - post_seq_len - len(self.mut.after) + offset
            post_seq_len -= offset

        #Fetch pre sequence
        pre_start = self.mut.pos - pre_seq_len
        pre_end = self.mut.pos
        pre_seq = extract_circular(genome, pre_start, pre_end)

        #Fetch post sequence
        post_start = self.mut.pos + len(self.mut.before)
        post_end = post_start + post_seq_len
        post_seq = extract_circular(genome, post_start, post_end)

        return pre_seq + self.mut.after + post_seq

    def optimise_folding(self, genome, threshold=-20.0):
        """Optimise oligo folding if dG below threshold"""
        if not self.oligo:
            oligo = self.make_oligo(genome, offset=0)
        else:
            oligo = self.oligo

        #Calculate default dG
        optimised_oligo = (oligo, ViennaRNA.mfe(str(oligo)), 0)
        candidates = [optimised_oligo]

        if ViennaRNA.mfe(str(oligo)) < threshold:
            #Offset towards 3'-end
            three_end_cap = rc.CONF['min_homology']
            end = self.oligo_len/2 - three_end_cap - len(self.mut.after)/2
            for i in range(1, end):
                candidate = self.make_oligo(genome, offset=i)
                candidates.append((candidate, ViennaRNA.mfe(str(candidate)), i))

            #Pick the best
            optimised_oligo = max(candidates, key=lambda x: x[1])

            #We are still below threshold. Panic, then offset towards 5'-end
            if optimised_oligo[1] < threshold:
                five_end_cap = rc.CONF['min_homology']
                start = -self.oligo_len/2 + five_end_cap + len(self.mut.after)/2
                for i in range(start, 0):
                    candidate = self.make_oligo(genome, offset=i)
                    candidates.append((candidate, ViennaRNA.mfe(str(candidate)), i))

                #Pick the new best (old candidates can still apply)
                optimised_oligo = max(candidates, key=lambda x: x[1])

        return optimised_oligo

    def target_lagging_strand(self, ori, ter):
        """Calculate replichore and target the lagging strand.

        ori and ter are ranges: [start, end]
        """
        #Check whether pos is not inside oriC or Ter
        if ori[0] <= self.pos <= ori[1]:
            self.replichore = 0
            log.warning("Oligo: {} inside oriC".format(self.short_id()))
        elif ter[0] <= self.pos <= ter[1]:
            self.replichore = 0
            log.warning("Oligo: {} inside termination region".format(self.short_id()))
        #OriC is on the second half of the genome
        elif ori[1] > ter[0]:
            if not ter[1] < self.pos < ori[0]:
                self.replichore = 1
            else:
                self.replichore = 2
        #OriC is on the first half of the genome
        else:
            if ori[1] < self.pos < ter[0]:
                self.replichore = 1
            else:
                self.replichore = 2

        #Scale replichore to target (strand 2 -> 1 and 1 -> -1)
        rp = int((self.replichore - 1.5) * 2)

        #Mismatch
        if self.replichore and rp != self.strand:
            self.oligo = reverse_complement(self.oligo)
            self.strand = self.strand * -1

    def add_barcodes(self, barcode_ids, barcoding_lib):
        """Add barcodes to an oligo, from barcoding_lib."""
        barcode_ids = barcode_ids.split('+')
        barcode_ids.reverse()
        for barcode in barcode_ids:
            self.barcodes_forward.insert(0, barcoding_lib[barcode]["forward"].lower())
            self.barcodes_reverse.append(barcoding_lib[barcode]["reverse"].lower())
            self.barcode_ids.insert(0, barcode)

    def output(self):
        """Return the oligo with barcodes"""
        out = self.barcodes_forward + [self.oligo] + self.barcodes_reverse
        out = map(str, out)
        return "".join(out)

    def id(self, id_type="auto"):
        """Calculate an ID based the various self-contained variables.

        id_type can be either full or auto.

        """
        id_type = id_type.lower()

        #Test id_type
        if id_type not in ["full", "auto"]:
            raise Exception("Unknown id_type: {}".format(id_type))

        #Auto
        if id_type == "auto" and self.custom_id:
            return self.custom_id

        idlst = [self.short_id()]

        if self.barcode_ids:
            idlst.append("BC:" + "+".join(self.barcode_ids))
        #Other stuff
        idlst.append(self.mut.__str__(idx=1))
        idlst.append(str(self.gene))
        idlst.append(self.code)
        idlst.append("RP" + str(self.replichore))
        #Custom id
        if self.custom_id:
            idlst.append("["+self.custom_id+"]")
        #Optimised folding
        if self.optimised:
            idlst.append("OPT:{}({:.1f})".format(self.optimised, self.dG_fold))

        return "_".join(idlst)

    def short_id(self):
        """Return id in the form of projectXXX"""
        return "{prj}{i:0>4}".format(
            prj = self.project if self.project else "",
            i = self.number)

    def set_custom_id(self, cid):
        """Set a custom ID."""
        self.custom_id = cid

    def copy(self):
        """Return a copy"""
        return deepcopy(self)


class Mutation:
    def __init__(self, before, after, pos, ref_seq=None):
        """A Mutation.

        Mutation objects are internally 0-indexed.

        """
        #Parse pos
        try:
            self.pos = int(pos)
        except ValueError:
            raise ValueError("Invalid position: '{}'.".format(pos))

        before = str(before)
        after  = str(after)
        #Parse and set before sequence
        if len(before) and not valid_dna(before):
            raise ValueError("Invalid 'before' sequence: '{}'.".format(before))
        self.before = before.upper()
        if ref_seq:
            if self.before != str(ref_seq[self.pos:self.pos+len(self.before)]):
                raise ValueError("Trying to mutate '{}', but found '{}' in reference sequence."
                                 "".format(self.before, ref_seq[self.pos:self.pos+len(before)]))

        #Parse and set after sequence
        if len(after) and not valid_dgn(after):
            raise ValueError("Invalid 'after' sequence: '{}'.".format(after))
        self.after = after

        #Automatically adjust after to all-lower if no lowercase chars are found.
        for n in "atgcrymkswbdhvn":
            if n in self.after:
                return
        self.after = self.after.lower()

    def MASC_primers(self, ref_genome, lengths=[200], temp=62.0, primer_c=2e-7, salt_c=0.05):
        """Design Multiplex Allele Specific Colony PCR primers.

        lengths is a list of PCR fragment lengths, and temp is the target tm.

        primer_c is primer concentration in mol/L, salt_c is monovalent cation
        concentration in mol/L.

        """
        if self.after and not valid_dna(self.after):
            raise Exception("Only strictly valid DNA (ATGC) can make MASC primers.")
        fpwt = ""
        fpmut = ""
        fp_offset = 0

        wtpos = self.pos+len(self.before)

        #In loop to shift frame if forward primers are identical.
        while str(fpwt) == str(fpmut):
            #fpwt: forward primer(wt)
            fpwt = make_primer(ref_genome, temp, wtpos-3+fp_offset, wtpos+fp_offset, salt_c, primer_c, 200)
            #fpmut: forward primer(mut)
            if self.after:
                if fp_offset > 0:
                    fpmut_ref = extract_circular(ref_genome, self.pos-200, self.pos)+self.after.upper()+extract_circular(ref_genome, self.pos, self.pos+fp_offset)
                else:
                    fpmut_ref = extract_circular(ref_genome, self.pos-200, self.pos)+self.after.upper()
                    #if len(self.after) > len(self.before):
                    #    fpmut_ref = fpmut_ref + ref_genome[self.pos+len(self.before)]
            else:
                fpmut_ref = extract_circular(ref_genome, self.pos-200+fp_offset, self.pos+fp_offset)+ref_genome[self.pos+fp_offset+len(self.before)]

            fpmut = make_primer(fpmut_ref, temp, len(fpmut_ref)-3, len(fpmut_ref), salt_c, primer_c, 200)
            #If primers are identical, shift one nt.
            fp_offset += 1

        avg_fp_len = int(float(len(fpwt[0])+len(fpmut[0]))/2)
        primers = {"fpwt": fpwt, "fpmut": fpmut}
        lengths = [lengths] if type(lengths) is int else lengths
        for l in lengths:
            ref = extract_circular(ref_genome, wtpos-1-avg_fp_len, wtpos-1-avg_fp_len+l)
            ref = reverse_complement(ref)
            pr = make_rev_primer(ref, temp, 0, 3, salt_c, primer_c, 200)

            primers[l] = pr

        return primers

    def __repr__(self):
        return "Mutation: [{}->{}] at pos {}".format(
            self.before,
            self.after,
            self.pos)

    def __str__(self, idx=0):
        return "[{}->{}].{}".format(self.before, self.after, self.pos+idx)

    def small_str(self, idx=0):
        before = self.before
        after = self.after
        if len(before) > 6:
            before = before[0:2] + ".." + before[-2:]
        if len(after) > 6:
            after = after[0:2] + ".." + after[-2:]

        return "[{}->{}].{}".format(before, after, self.pos+idx)

    def copy(self):
        """Return a copy"""
        return deepcopy(self)


class Gene:
    """Defines a single gene and all relevant information

    Instantiation is straightforward:

        >>> gene = Gene("ficX", 20, 1, "ATGTCTGCAACAAAACTG", "TTTTTGGAATGAGCT")
        >>> gene
        ficX

    Getting a part of a gene sequence can be done in a couple of ways. Using the
    normal int or slice to get a genes cds:

        >>> gene[0:3]
        'ATG'
        >>> gene[0]
        'A'

    However, a part of a gene can also be obtained by using a tuple instead of a
    slice, which can include the leader:

        >>> gene[0,3]
        'ATG'
        >>> gene[-3,0]
        'GCT'
        >>> gene[-6,-2]
        'TGAG'
        >>> gene[-3,3]
        'GCTATG'

    None can be used to get the 'rest':

        >>> gene[0,]
        'ATGTCTGCAACAAAACTG'
        >>> gene[None, 3]
        'TTTTTGGAATGAGCTATG'

    """

    def __init__(self, name, pos, strand, cds, leader=None, promoter=None,
                 promoter_pos=None):
        self.name = name
        self.pos = pos
        self.strand = strand

        #Coding sequence
        self.cds = cds
        if not isinstance(cds, Sequence):
            self.cds = Sequence(cds)

        #Leader sequence
        self.leader = leader
        if leader and not isinstance(leader, Sequence):
            self.leader = Sequence(leader)

        self.leader_pos = pos-len(leader) if leader else None

        #Promoter sequence
        self.promoter = promoter
        if promoter and not isinstance(promoter, Sequence):
            self.promoter = Sequence(promoter)

        self.promoter_pos = promoter_pos
        self.in_operon = False

    def __str__(self):
        return self.name

    def __repr__(self):
        return self.__str__()

    def do_mutation(self, mut):
        """Return Mutation in genome context (positive strand)"""
        mut = mut.copy()

        if str(mut.before) != str(self[mut.pos,mut.pos+len(mut.before)]):
            extended = "{}({}){}".format(self[mut.pos-3,mut.pos],
                                         self[mut.pos,mut.pos+len(mut.before)],
                                         self[mut.pos+len(mut.before), mut.pos+len(mut.before)+3])
            raise MutationError("Trying to mutate {} but found in gene: {}".format(str(mut.before), extended))

        if self.strand == 1:
            mut.pos = mut.pos+self.pos
            return mut
        else:
            mut.pos = self.pos - mut.pos - len(mut.before)
            mut.before = reverse_complement(mut.before)
            mut.after = reverse_complement(mut.after)

        return mut

    def __getitem__(self, slc):
        """Get a part of the gene sequence"""
        if type(slc) is slice or type(slc) is int:
            return self.cds[slc]
        elif type(slc) is tuple:
            #Parse tuple
            if len(slc) == 1:
                st, end = slc[0], len(self.cds)
            elif len(slc) == 2:
                st, end = slc
                if st is None: st = -len(self.leader)
            else:
                raise ValueError("Invalid tuple for Gene: {}".format(slc))

            #Aquire leader/cds
            if st < 0 and end <= 0:
                if end == 0:
                    return self.leader[st:]
                return self.leader[st:end]
            elif st < 0 < end:
                return self.leader[st:] + self.cds[:end]
            elif st >= 0 and end > 0:
                return self.cds[st:end]
            elif st == end:
                return ''
            else:
                raise ValueError("Invalid tuple for Gene: {}".format(slc))
        else:
            raise TypeError("Getting a sequence from a gene can be a slice, int or tuple.")

    def __len__(self):
        return len(self.leader) + len(self.cds)

    def copy(self):
        """Return a copy"""
        return deepcopy(self)


class Sequence:
    """A DNA sequence object.

    A Sequence is a simple DNA sequence object with various mutation methods.

        >>> s = Sequence("ATAGCACATA")

    The mutation methods keep track of what positions are mutated, and attempts
    to redo and create the best possible 'alignment' while the sequence is
    being mutated.

    TODO: mut example
    TODO: delete example

    Insertions:

        >>> seq = s.copy()
        >>> seq.insert("G", 1)
        AGTAGCACATA
        >>> seq.insert("T", 1)
        ATGTAGCACATA

    Insertions at either end is currently not supported:

        >>> seq.insert("G", 0)
        False
        >>> seq.insert("G", 12)
        False

    """
    def __init__(self, s):
        seq = str(s).upper()
        self.org_seq = seq
        self.org_len = len(seq)

        self.type = "DNA"
        if "U" in seq:
            self.type = "RNA"
            seq = seq.replace("U", "T")

        if not valid_dna(seq):
            raise ValueError("Invalid {}: {}.".format(self.type, s))

        self.seq = {i: [[nt, nt]] for i, nt in enumerate(seq)}
        self.seql = [(i, 0) for i in range(len(seq))]
        self.wobbles = list()
        self.mutations = {i: ["n"] for i in range(len(seq))}

    """
    Wobble methods
    ~~~~~~~~~~~~~~
    """

    def add_wobble(self, wbseq, pos):
        """Add a wobble sequence to the sequence."""
        wbl = SequenceWobble(wbseq, pos, len(self))
        if not self.test_wobble(wbl):
            msg="Wobble error output: \n{:<4d} ".format(wbl.start)
            p_msg, w_msg, e_msg = "", "", ""
            for i in wbl.range():
                p_msg += self[i]
                w_msg += nts_to_dgn[wbl[i]]
                e_msg += " " if self[i] in wbl[i] else "*"
                if (i+1-wbl.start) % 60 == 0:
                    msg += "{}\n     {}\n     {}\n\n{:<4d} ".format(p_msg, w_msg, e_msg, i+1)
                    p_msg, w_msg, e_msg = "", "", ""
            msg += "{}\n     {}\n     {}".format(p_msg, w_msg, e_msg, i+1)
            log.debug(msg)
            raise WobbleError("Wobble does not match seq. Check log for more output.")

        for i in range(len(self.wobbles)):
            o_wbl = self.wobbles[i]
            if is_inside(o_wbl.start, o_wbl.end, wbl.start, wbl.end):
                self.wobbles[i] += wbl
                return

        self.wobbles.append(wbl)
        self.wobbles.sort()

    def test_wobble(self, wbl):
        """Test a wobble sequence against the sequence."""
        for i in wbl.range():
            if self[i] not in wbl[i]:
                return False
        return True

    def get_wobble(self, i):
        """Returns a list of possible codons at pos i."""
        for wbl in self.wobbles:
            if i in wbl:
                return wbl[i]

        #Run-out, return False to distinguish between N and no wbl.
        return False

    def get_wobble_str(self, not_wbl="N", dels=""):
        """Get wobbles represented as a string.

        Free positions marked as not_wbl, 'N' by default.

            >>> s = Sequence("AACGGCACA")
            >>> s.add_wobble( "ACNGC", 1)
            >>> s.get_wobble_str()
            'NACNGCNNN'
            >>> s.get_wobble_str("*")
            '*ACNGC***'

        Insertions can be marked as well, with dels:

            >>> s.delete(7)
            AACGGCAA
            >>> s.get_wobble_str(" ", "-")
            ' ACNGC - '

        """
        wbl_str = list()

        for a in range(self.org_len):
            for b in range(len(self.seq[a])):
                if self.mutations[a][b] == "d":
                    wbl_str.append(dels)
                else:
                    c = self.seql.index((a, b))
                    for wbl in self.wobbles:
                        if c in wbl:
                            wbl_str.append(str(wbl)[c])
                            break
                    else:
                        wbl_str.append(not_wbl)

        return "".join(wbl_str)

    def reverse_complement(self):
        """Reverse complement a Sequence including wobbles."""
        if self.get_mutated_positions():
            raise NotImplementedError("A mutated sequence can not (yet) be reverse complemented.")

        new_seq = reverse_complement(str(self.org_seq))
        new = Sequence(new_seq)

        for wbl in self.wobbles:
            new_wbl_seq = reverse_complement(wbl.seq)
            new_wbl_start = self.org_len - len(wbl.seq) - wbl.start
            new.add_wobble(new_wbl_seq, new_wbl_start)

        return new

    """
    Mutation methods
    ~~~~~~~~~~~~~~~~
    """

    def mutate(self, nt, a, b=None, in_place=True):
        """Do a point mutation.

        Supply only ``a`` to mutate an absolute position:

            >>> s = Sequence("ATGC")
            >>> s.insert("A", 3, in_place=True)
            ATGAC
            >>> s.mutate("T", 3, in_place=False)
            ATGTC

        Supply both ``a`` and ``b`` to mutate in a specific position:

            >>> s.mutate("T", 3, 0, in_place=False)
            ATGAT

        Mutation methods return False if mutation is not possible, returns the
        new sequence otherwise:

            >>> s.mutate("A", 0, in_place=False)
            False

        """
        if nt not in ["A", "T", "G", "C"] or len(nt) != 1:
            raise MutationError("{} is not a single DNA nucleotide!".format(nt))

        if b is None:
            c = a
            a, b = self.seql[a]
        else:
            #Nonexistant
            if (a, b) not in self.seql:
                return False
            c = self.seql.index((a, b))

        if self.seq[a][b][0] == nt:
            return False

        new = self
        if not in_place:
            new = self.copy()

        for wbl in new.wobbles:
            if c in wbl and nt not in wbl[c]:
                return False

        new.seq[a][b][0] = nt
        if b == 0:
            if new.seq[a][b][1] == nt:
                new.mutations[a][b] = "n"
            else:
                new.mutations[a][b] = "m"

        return new

    def delete(self, a, b=None, in_place=True):
        """Delete a position.

        Behaves similarly to mutate()

        """
        if b is None:
            c = a
            a, b = self.seql[a]
        else:
            #Nonexistant
            if (a, b) not in self.seql:
                return False
            c = self.seql.index((a, b))

        #Inside wobble, abort
        for wbl in self.wobbles:
            if c in wbl:
                return False

        #Already deleted
        if self.seq[a][b][0] == "-":
            return False

        #Maybe take a copy
        new = self.copy() if not in_place else self

        for wbl in new.wobbles:
            if wbl.lt(c): wbl.offset -= 1

        #Deleting insertion
        if len(new.seq[a]) > 1:
            org_nt = new.seq[a][b][1]
            #delete sequence and mutation position
            del(new.seq[a][b])
            del(new.mutations[a][b])
            if b == 0:
                #Convert old insertion to mutation
                if new.seq[a][b][0] != org_nt:
                    new.mutations[a][b] = "m"
                #Insertion is equal to origital nt
                else:
                    new.mutations[a][b] = "n"
                new.seq[a][b][1] = org_nt

            #Delete sequence list position
            #while (a, b) in new.seql:
            #    b += 1
            #    c += 1
            #del(new.seql[c-1])
            new.recreate_seql()
            return new

        #Prev mut is an insertion
        #Delete insertion and convert current position to a mutation
        if a > 1 and len(new.seq[a-1]) > 1:
            e = len(new.seq[a-1]) - 1
            nt = new.seq[a-1][e][0]
            del(new.seq[a-1][e])
            #del(new.seql[c-1])
            del(new.mutations[a-1][e])
            m = "m"
            #Reverted to original
            if nt == new.seq[a][b][1]:
                m = "n"
            new.seq[a][b][0] = nt
            new.mutations[a][b] = m
            new.recreate_seql()
            return new

        #Normal deletion
        new.seq[a][b][0] = "-"
        #del(new.seql[c])
        new.mutations[a] = ["d"]
        new.recreate_seql()
        return new

    def insert(self, nt, a, b=None, in_place=True):
        """Insert a nt as pos.

        Behaves similarly to mutate()

        """
        if nt not in ["A", "T", "G", "C"] or len(nt) != 1:
            raise MutationError("{} is not a single DNA nucleotide!".format(nt))

        if b is None:
            c = a
            try:
                a, b = self.seql[a]
            except IndexError:
                if c in range(self.org_len):
                    a += 1
                    b = 0
                # elif c == len(self):
                #     a, b = self.seql[-1]
                #     b += 1
                else:
                    return False
        else:
            raise NotImplementedError("Coordinates not possible for insertions.")

        #Inserting as very first codon
        if c == 0 and self.seql[0] == (0, 0):
            return False

        #Inside wobble, abort
        for wbl in self.wobbles:
            if c in wbl and c-1 in wbl:
                return False

        new = self
        if not in_place:
            new = self.copy()

        for wbl in new.wobbles:
            if wbl.lt(c) or wbl.lt(c - 1): wbl.offset += 1

        #Inserting into deletion
        if a and new.seq[a-1][0][0] == "-":
            e = a - 1
            while e != 0 and new.seq[e-1][0][0] == "-":
                e -= 1
            #Convert deletion to insertion
            new.seq[e][0][0] = nt
            if new.seq[e][0][1] != nt:
                new.mutations[e] = ["m"]
            else:
                new.mutations[e] = ["n"]
            #new.seql.insert(c, (e, 0))
            new.recreate_seql()
            return new

        #Inserting into insertions
        elif b > 0:
            e = len(new.seq[a])
            new.seq[a].insert(b, [nt, "-"])
            new.mutations[a].insert(b, "i")
            #new.seql.insert(c-b+e, (a, e))
            new.recreate_seql()
            return new

        #No other insertions
        elif len(new.seq[a-1]) == 1:
            new.seq[a-1].append([nt, "-"])
            new.mutations[a-1].append("i")
            #new.seql.insert(c, (a-1, 1))
            new.recreate_seql()
            return new

        else:
            #Inserting into end of insertions
            e = len(new.seq[a-1])
            new.seq[a-1].append([nt, "-"])
            new.mutations[a-1].append("i")
            #new.seql.insert(c, (a-1, b+e))
            new.recreate_seql()
            return new

    def do_mutations(self, *muts):
        """Perform a list of mutations in the form of ("[m|i|d]", p, nt).

        Any number of mutations can be supplied:

            >>> s = Sequence("ATGC")
            >>> f = s.do_mutations(("m", 0, "T"), ("m", 1, "T"), ("m", 2, "T"))
            >>> str(s)
            'TTTC'

        Returns number of failed mutations:

            >>> f
            1

        """
        f = 0
        for mut in muts:
            if mut[0] == "m":
                if not self.mutate(mut[2], mut[1], in_place=True): f += 1
            elif mut[0] == "i":
                if not self.insert(mut[2], mut[1], in_place=True): f += 1
            elif mut[0] == "d":
                if not self.delete(mut[1], in_place=True): f += 1
            else:
                raise MutationError("Invalid mutation code: {}".format(mut[0]))

        return f

    def optimise_mutations(self, penalty_mid=(2,3,3), in_place=True):
        """Try to remove insertions and deletions for a better alignment.

        When doing lots of deletions and insertions, one may end up with a
        sequence that has a lot of mutations in it, rather than insertions/
        deletions:

            >>> s = Sequence("AACGGCACAGTTGCTGGAAATTGCAGGAGTTGGCG")
            >>> s.do_mutations(("i", 18, "A"), ("i", 18, "T"), ("i", 18, "C"),
            ...                ("m", 21, "T"), ("m", 22, "T"), ("d", 23), ("d", 23))
            0
            >>> str(s)
            'AACGGCACAGTTGCTGGACTATTGCAGGAGTTGGCG'
            >>> s.get_mutation()
            Mutation: [AATT->ctatt] at pos 18
            >>> s = s.optimise_mutations()
            >>> str(s)
            'AACGGCACAGTTGCTGGACTATTGCAGGAGTTGGCG'
            >>> s.get_mutation()
            Mutation: [A->ct] at pos 18

        """
        insr = self.get_insertions()
        dels = self.get_deletions()
        if not insr or not dels:
            return self

        st = min(insr[0][0], dels[0][0])
        #TODO: This is an inelegant and cheapo fix
        if st == 0: st = 1
        nd = max(insr[-1][0], dels[-1][0]) + 1
        mut_seq = ["0"] + [nt[0] for a in range(st, nd) for nt in self.seq[a] if nt[0] != "-"]
        org_seq = ["0"] + [nt[1] for a in range(st, nd) for nt in self.seq[a] if nt[1] != "-"]

        p_mis, p_ins, p_del = penalty_mid

        mat = list()#[[(0, -1)]]
        for i, m in enumerate(mut_seq):
            row = list()
            mat.append(row)
            for j, o in enumerate(org_seq):
                if i > 0 and j > 0:
                    score = 0 if m == o else p_mis
                    c_mis = (mat[i-1][j-1][0] + score, 0)
                    c_ins = (mat[i][j-1][0]   + p_ins, 1)
                    c_del = (mat[i-1][j][0]   + p_del, 2)
                    row.append(min(c_mis, c_ins, c_del))
                elif i == 0 and j == 0:
                    row.append((0, 0))
                elif i == 0:
                    row.append((mat[i][j-1][0]   + p_ins, 1))
                elif j == 0:
                    row.append((mat[i-1][j][0]   + p_del, 2))

        new = self
        if not in_place:
            new = self.copy()

        for a in range(st, nd):
            new.seq[a] = [["-", "-"]]
            new.mutations[a] = ["n"]
        b = 0
        i, j = len(mat)-1, len(mat[0])-1

        while i or j:
            score, v = mat[i][j]
            a = nd-len(org_seq)+j
            if v == 0:
                b = 0
                new.seq[a][b] = [mut_seq[i], org_seq[j]]
                if mut_seq[i] != org_seq[j]:
                    new.mutations[a][b] = "m"
                i -= 1
                j -= 1
            elif v == 1:
                b = 0
                new.seq[a][b] = ["-", org_seq[j]]
                new.mutations[a][b] = "d"
                j -= 1
            elif v == 2:
                b += 1
                new.seq[a].insert(1, [mut_seq[i], "-"])
                new.mutations[a].insert(1, "i")
                i -= 1

        #new.seql = [(a, b) for a in range(new.org_len)
        #                   for b in range(len(new.seq[a]))
        #                   if new.mutations[a][b] != "d"]
        new.recreate_seql()

        return new

    def random_mutation(self, max_mut=None, p_di=(0.2, 0.2)):
        """Perform a random mutation.

        Use p_di to set the probabilities of deletions and insertions. The
        probability of a point mutation is 1 - sum(p_di).

        """
        thr = sum(p_di)
        if thr > 1:
            raise ValueError("Sum of probabilities can not be above 1.")

        #Which pos are we mutating
        pos = random.randint(0, len(self)-1)

        #Determine if we are in wobble
        is_wbl = self.get_wobble(pos)

        #Number to decide deletion, insertion or mutation
        dim = random.random() if not is_wbl else 2

        #Insertion
        if dim < p_di[0]:
            if pos == 0: pos = 1
            nt = random.choice(["A", "T", "G", "C"])
            #print(pos, nt)
            new = self.insert(nt, pos, in_place=False)
            if max_mut:
                muts = self.get_mutated_positions()
                #Remove a random mutation
                while len(muts) > int(max_mut)-1:
                    r = random.choice(range(len(muts)))
                    new.revert_mutation(*muts[r], in_place=True)
                    del(muts[r])
            return new
        #Deletion
        elif dim < thr:
            #print(pos)
            new = self.delete(pos, in_place=False)
            if max_mut:
                muts = new.get_mutated_positions()
                t = self.seql[pos]
                if t in muts: del(muts[muts.index(t)])
                #Remove a random mutation
                while len(muts) > int(max_mut)-1:
                    r = random.choice(range(len(muts)))
                    new.revert_mutation(*muts[r], in_place=True)
                    del(muts[r])
            return new
        #Mutation
        else:
            for i in xrange(len(self)*2):
                candidates = list(is_wbl) if is_wbl else ["A", "T", "G", "C"]
                del(candidates[candidates.index(self[pos])])
                #If there are no possible mutations
                if not candidates:
                    pos = random.randint(0, len(self)-1)
                    is_wbl = self.get_wobble(pos)
                    continue

                nt = random.choice(candidates)
                new = self.mutate(nt, pos, in_place=False)
                if max_mut:
                    muts = new.get_mutated_positions()
                    t = self.seql[pos]
                    if t in muts: del(muts[muts.index(t)])
                    #Remove a random mutation
                    while len(muts) > int(max_mut)-1:
                        r = random.choice(range(len(muts)))
                        new.revert_mutation(*muts[r], in_place=True)
                        del(muts[r])
                return new

            #This is SO going to be optimised.
            return False

    def revert_mutation(self, a, b, in_place=True):
        """Revert a mutation."""
        new = self
        if not in_place:
            new = self.copy()

        if new.mutations[a][b] == "n":
            return new
        elif new.mutations[a][b] == "m":
            new.mutations[a][b] = "n"
            new.seq[a][b][0] = new.seq[a][b][1]
        elif new.mutations[a][b] == "d":
            new.mutations[a][b] = "n"
            new.seq[a][b][0] = new.seq[a][b][1]
            #c = 0
            #while c < len(new.seql) and (a, b) > new.seql[c]:
            #    c += 1
            #new.seql.insert(c, (a, b))
            new.recreate_seql()
        elif new.mutations[a][b] == "i":
            del(new.mutations[a][b])
            del(new.seq[a][b])
            #while (a, b) in new.seql:
            #    b += 1
            #del(new.seql[new.seql.index((a, b-1))])
            new.recreate_seql()

        return new

    def recreate_seql(self):
        """Regenerate the sequence list."""
        self.seql = [(a, b) for a in range(self.org_len)
                            for b in range(len(self.seq[a]))
                            if self.mutations[a][b] != "d"]

    """
    Other
    ~~~~~
    """

    def get_mutation(self):
        """Returns a mutation object relative to the sequence"""
        muts = self.get_mutated_positions()
        if not muts:
            return False

        insertions = len(self.get_insertions())
        deletions  = len(self.get_deletions())

        start = sum(muts[0])
        end   = muts[-1][0]

        #Before mutation
        before = self.org_seq[start:end+1]

        #After mutation. Lowercase all mutations
        after  = ""
        for p in self.seql[start:end+insertions-deletions+1]:
            a, b = p
            if p in muts:
                after += self.seq[a][b][0].lower()
            else:
                after += self.seq[a][b][0].upper()

        #Mutation object
        mut = Mutation(before, after, start)
        return mut

    def get_mutated_positions(self):
        """Return a list coordinates for all mutated positions"""
        return [(a,b) for a in range(self.org_len)
                      for b in range(len(self.seq[a]))
                      if  self.mutations[a][b] != "n"]

    def n_muts(self):
        """Return number of mutations."""
        m = 0
        for a in range(self.org_len):
            for b in range(len(self.seq[a])):
                if self.mutations[a][b] != "n": m += 1
        return m

    def get_deletions(self):
        """Return a list of coordinates for all deletions"""
        return [(a,b) for a in range(self.org_len)
                      for b in range(len(self.seq[a]))
                      if  self.mutations[a][b] == "d"]

    def get_insertions(self):
        """Return a list of coordinates for all insertions"""
        return [(a,b) for a in range(self.org_len)
                      for b in range(len(self.seq[a]))
                      if  self.mutations[a][b] == "i"]

    def get_mutated_str(self):
        """Mutated string with deletions marked with '-'. """
        return "".join([nt[0] for i, nts in sorted(self.seq.items()) for nt in nts])

    def get_mutation_str(self):
        """A string of one-symbol mutations."""
        l = list()
        for a in sorted(self.seq):
            for b in range(len(self.seq[a])):
                if b == 0 and self.mutations[a][b] == "n":
                    l.append(" ")
                else:
                    l.append(self.mutations[a][b])

        return "".join(l)

    def get_original_str(self):
        """Original sequence with insertions marked with '-'. """
        return "".join([nt[1] for i, nts in sorted(self.seq.items()) for nt in nts])

    def __str__(self):
        return "".join([str(nt) for nt in self])

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        for a, b in self.seql:
            yield self.seq[a][b][0]

    def __getitem__(self, a):
        """getitem.

        It is possible to use coordinates:

            >>> s = Sequence("ATCG")
            >>> s[1,0]
            'T'

        Integers:

            >>> s[1]
            'T'

        And slices:

            >>> s[1:3]
            'TC'

        """
        if type(a) is tuple:
            a, b = a
            return self.seq[a][b][0]
        else:
            out = list()
            if type(a) is int:
                a = slice(a, a + 1)
            for a, b in self.seql[a]:
                out.append(self.seq[a][b][0])
            return "".join(out)

    def __len__(self):
        return len(self.seql)

    def copy(self):
        """A function that can be optimised a lot."""
        return deepcopy(self)

    def pprint(self, org=True):
        """Pretty print sequence with wobbles.

        Set org=True, to print original sequence, set to False to print only
        mutated sequence.

        """
        if org:
            mut_str = self.get_mutation_str()
            print("MUT:", mut_str)
            print("ORG:", self.get_original_str())
            print("NEW:", self.get_mutated_str())
            if self.wobbles:
                print("WBL:", self.get_wobble_str(" ", "_"))
        else:
            print(self)
            gp = 0
            for wbl in self.wobbles:
                print(" "*(wbl.start - gp + wbl.offset) + str(wbl), end="")
                gp = wbl.end
            print()


class SequenceWobble:
    def __init__(self, wbseq, pos, seq_len=None):
        self.seq = str(wbseq).upper()
        if self.seq and not valid_dgn(self.seq):
            raise ValueError("Invalid wobble sequence: {}".format(wbseq))

        self.start = pos
        if pos < 0:
            self.seq = self.seq[-pos:]
            self.start = 0

        self.end = pos + len(wbseq)
        if seq_len and self.end > seq_len:
            self.seq = self.seq[:seq_len - pos]
            self.end = seq_len

        self.stop = seq_len if seq_len else self.end
        self.offset = 0
        self.allowed = [dgn_to_nts[nt] for nt in self.seq]

    def slc(self):
        return self.start+self.offset, self.end+self.offset

    def range(self): return range(*self.slc())
    def lt(self, i): return i < self.start + self.offset

    def __lt__(self, other): return self.start <  other.start
    def __gt__(self, other): return self.start >  other.start
    def __le__(self, other): return self.start <= other.start
    def __ge__(self, other): return self.start >= other.start

    def __str__(self):
        """Returns a representation of the wobble in Sequence space.

            >>> w = SequenceWobble("ATC", 3, 9)
            >>> str(w)
            '---ATC---'

        """
        #Span including offset
        st, nd = self.slc()
        return "-"*st + self.seq + "-"*(self.stop - nd)

    def __repr__(self):
        fmtargs = (self.start, self.end, self.offset, self.seq)
        return "wobble({}..{}+{}: {})".format(*fmtargs)

    def __getitem__(self, i):
        return self.allowed[i - self.start - self.offset]

    def __contains__(self, i):
        return self.start + self.offset <= i < self.end + self.offset

    def __add__(self, other):
        """Add two wobbles and find common consensus.

            >>> w1 = SequenceWobble("CDGN", 2)
            >>> w2 = SequenceWobble( "NGCNCA", 3)
            >>> w1 + w2
            wobble(2..9+0: CDGCNCA)

            >>> w3 = SequenceWobble("CACNGCNCA", 0)
            >>> w3
            wobble(0..9+0: CACNGCNCA)
            >>> w1 + w3
            wobble(0..9+0: CACDGCNCA)

        """
        #Wobble start, intersection start
        st, i_st = sorted([self, other], key=lambda w: w.start)
        #Intersection end, wobble end
        i_nd, nd = sorted([self, other], key=lambda w: w.end)

        seq = st.seq[:i_st.start - st.start]

        for i in range(i_st.start, i_nd.end):
            common = self[i] & other[i]
            if not common: raise ValueError("Incompatible wobble sequences: {}, {}".format(self, other))
            seq += nts_to_dgn[common]

        seq += nd.seq[i_nd.end - nd.start:nd.end - nd.start]

        return SequenceWobble(seq, st.start)


class MutationError(ValueError): pass
class WobbleError(ValueError): pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()
