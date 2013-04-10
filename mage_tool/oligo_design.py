#!/usr/bin/env python

"""
Module for designing DNA oligos from mutation objects
"""

from __future__ import print_function

import logging
from copy import deepcopy

import ViennaRNA
from helpers import valid_dna
from helpers import valid_rna
from helpers import valid_dgn
from helpers import reverse_complement
from helpers import dgn_to_nts

#Define a log
log = logging.getLogger("MODEST.oligo")
log.addHandler(logging.NullHandler())


class Oligo:
    """Oligo Object"""

    def __init__(self, mut, gene, project=None, number=0, oligo_len=90):
        self.mut = mut
        self.gene = gene
        self.oligo_len = oligo_len
        self.pos = mut.pos
        self.oligo = None
        self.optimised = 0
        self.dG_fold = 0
        self.replichore = None
        self.project = project.replace(" ", "_")
        self.number = number

        self.barcodes_forward = list()
        self.barcodes_reverse = list()
        self.barcode_ids = list()
        self.code = ""
        self.operation = ""
        self.operation_values = list()

    def set_oligo(self, genome, optimise=True, threshold=-20):
        if optimise:
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
            extended = "{}({}){}".format(genome[mutpos-5:mutpos],
                                         found,
                                         genome[mutpos+len_bef:mutpos+len_bef+5])
            raise Exception("Trying to mutate {}, but found {} in genome. {}".format(self.mut.before, found, extended))

        #Calculate flanking sequence lengths
        post_seq_len = (self.oligo_len - len(self.mut.after))/2
        pre_seq_len = self.oligo_len - post_seq_len - len(self.mut.after) + offset
        post_seq_len -= offset

        #Fetch pre sequence
        pre_seq = genome[self.mut.pos-pre_seq_len:self.mut.pos]

        #Fetch post sequence
        post_seq_start = self.mut.pos+len(self.mut.before)
        post_seq = genome[post_seq_start : post_seq_start+post_seq_len]
        
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
            three_end_cap = 15
            end = self.oligo_len/2 - three_end_cap - len(self.mut.after)
            for i in range(1, end):
                candidate = self.make_oligo(genome, offset=i)
                candidates.append((candidate, ViennaRNA.mfe(str(candidate)), i))

            #Pick the best
            optimised_oligo = max(candidates, key=lambda x: x[1])

            #We are still below threshold. Panic, then offset towards 5'-end
            if optimised_oligo[1] < threshold:
                five_end_cap = 15
                start = -self.oligo_len/2 + five_end_cap
                for i in range(start, 0):
                    candidate = self.make_oligo(genome, offset=i)
                    candidates.append((candidate, ViennaRNA.mfe(str(candidate)), i))

                #Pick the new best (old candidates can still apply)
                optimised_oligo = max(candidates, key=lambda x: x[1])

        return optimised_oligo

    def target_lagging_strand(self, ori, ter):
        """Calculate replichore and target the lagging strand

        Returns False if pos is inside oriC or Ter, True otherwise
        ori and ter are ranges: [start, end]
        """
        #Check whether pos is not inside oriC or Ter
        if ori[0] <= self.pos <= ori[1]:
            self.replichore = 0
            log.warning("Oligo: {} inside oriC".format(self.short_id()))
            return False
        elif ter[0] <= self.pos <= ter[1]:
            self.replichore = 0
            log.warning("Oligo: {} inside termination region".format(self.short_id()))
            return False

        #OriC is on the second half of the genome
        if ori[1] > ter[0]:
            #Replichore 1, reverse complement
            if not ter[1] < self.pos < ori[0]:
                self.replichore = 1
                self.oligo = reverse_complement(self.oligo)
                return True
            else:
                self.replichore = 2
                return True
        #OriC is on the first half of the genome 
        else:
            #Replichore 1, reverse complement
            if ori[1] < self.pos < ter[0]:
                self.replichore = 1
                self.oligo = reverse_complement(self.oligo)
                return True
            else:
                self.replichore = 2
                return True
        
    def add_barcodes(self, barcode_ids, barcoding_lib):
        """TODO"""
        #Add barcode IDs to self.barcode_ids
        #and add barcode sequences to self.barcodes_*
        #forward barcodes should be prepended, and backward barcodes should be appended
        #Maybe add some primer checking too (maybe second function)
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

    def id(self):
        """Calculate an ID based the various self-contained variables"""
        return "{s_id}{barcodes}_{mut}_{gene}_{code}_RP{rp}{opt}".format(
            s_id = self.short_id(),
            barcodes = "_BC:" + "+".join(self.barcode_ids) if self.barcode_ids else "",
            mut = self.mut.__str__(idx=1),
            gene = str(self.gene),
            code = self.code,
            rp = self.replichore,
            opt = "_OPT:{}({:.1f})".format(str(self.optimised), self.dG_fold) if self.optimised else ""
            )

    def short_id(self):
        """Return id in the form of projectXXX"""
        return "{prj}{i:0>4}".format(
            prj = self.project if self.project else "",
            i = self.number
            )

    def copy(self):
        """Return a copy"""
        return deepcopy(self)


class Mutation:
    """Mutation object

    IMPORTANT: Position must be 0-indexed!
    
    mut_format can be:
    arrow: position required
        point mutation: A->T, pos, mut_change="Point_Mutation"
        insertion:      A, pos, mut_change="Insertion"
        deletion        3, pos, mut_change="Deletion"
    eq: position required
        point mutation: A=T, pos
        insertion:      =AT, pos
        deletion:       A=, pos
    genome: search genome for mutation (eq format)
        AATGATA[ATG=GT]ATGATA
        
    """
    def __init__(self, mut_format, mut, pos=0, mut_type="", ref_genome=False):
        if mut_format.lower() == "arrow":
            self._parse_arrow(mut, int(pos), mut_type, ref_genome)
        elif mut_format.lower() == "eq":
            self._parse_eq(mut, int(pos))
        else:
            raise Exception("Format: \"{}\" unknown".format(mut_format))

        for n in ["a", "t", "g", "c"]:
            if n in self.after:
                return
        self.after = self.after.lower()
    
    def __repr__(self):
        return "Mutation: [{}={}] at pos {}".format(
            self.before,
            self.after,
            self.pos)

    def __str__(self, idx=0):
        return "[{}={}].{}".format(self.before, self.after, self.pos+idx)
    
    def small_str(self, idx=0):
        before = self.before
        after = self.after
        if len(before) > 6:
            before = before[0:2] + ".." + before[-2:]
        if len(after) > 6:
            after = after[0:2] + ".." + after[-2:]

        return "[{}={}].{}".format(before, after, self.pos+idx)
    
    def _parse_arrow(self, mut, pos, mut_type, ref_genome=False):
        """Parse arrow format"""
        self.pos = pos
        if mut_type.lower() == "point_mutation":  
            self.before = mut.split("->")[0]
            self.after = mut.split("->")[1]
        elif mut_type.lower() == "insertion":
            self.before = ""
            self.after = mut
        elif mut_type.lower() == "deletion" or mut_type.lower() == "large_deletion":
            if not ref_genome:
                raise Exception("No reference genome supplied")
            self.before = ref_genome[pos-1:pos+int(mut)-1]
            self.after = ""
        else:
            raise Exception("Unknown mut_type: " + mut_type)
            
    def _parse_eq(self, mut, pos):
        self.pos = pos
        mut = mut.strip("[]")
        self.before, self.after = mut.split("=")

    def copy(self):
        """Return a copy"""
        return deepcopy(self)


class Gene:
    """Defines a single gene and all relevant information

    Instantiation is straightforward:

        >>> gene = Gene("ficX", 20, 1, "ATGTCTGCAACAAAACTG", "TTTTTGGAATGAGCT")
        >>> gene
        ficX

    If wobble sequence is not supplied it is generated automatically:

        >>> gene.leader_wobble
        'NNNNNNNNNNNNNNN'

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

    def __init__(self, name, pos, strand, cds, leader, leader_wobble=None,
                 promoter=None, promoter_pos=None):
        self.name = name
        self.pos = pos
        self.strand = strand
        self.cds = cds
        self.leader = leader
        self.leader_pos = pos-len(leader)
        self.promoter = promoter
        self.promoter_pos = promoter_pos

        self.leader_wobble = leader_wobble
        if not leader_wobble:
            self.leader_wobble = "N"*len(self.leader)

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
    """Sequence object"""
    def __init__(self, s):
        seq = str(s).upper()
        self.org_seq = seq

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

    class Wobble:
        def __init__(self, wbseq, pos, seq_len):
            self.seq = str(wbseq).upper()
            if not valid_dgn(self.seq):
                raise ValueError("Invalid wobble sequence: {}".format(wbseq))

            self.start  = pos
            if pos < 0:
                self.seq = self.seq[-pos:]
                self.start = 0

            self.end    = pos + len(wbseq)
            self.offset = 0
            self.stop   = seq_len
            self.allowed = [dgn_to_nts[nt] for nt in self.seq]

        def slc(self):
            return self.start+self.offset, min(self.end+self.offset, self.stop)
        
        def range(self): return range(*self.slc())
        def lt(self, i): return i < self.start + self.offset

        def __lt__(self, other): return self.start <  other.start
        def __gt__(self, other): return self.start >  other.start
        def __le__(self, other): return self.start <= other.start
        def __ge__(self, other): return self.start >= other.start

        def __str__(self): return self.seq
        def __repr__(self):
            fmtargs = (self.start, self.end, self.offset, self.seq)
            return "wobble({}..{}+{}: {})".format(*fmtargs)

        def __getitem__(self, i):
            return self.allowed[i - self.start + self.offset]

        def __contains__(self, i):
            return self.start + self.offset <= i < self.end + self.offset

    def add_wobble(self, wbseq, pos):
        wbl = self.Wobble(wbseq, pos, len(self))
        if not self.test_wobble(wbl):
            raise ValueError("Wobble does not match seq: {}".format(wbseq))

        self.wobbles.append(wbl)
        self.wobbles.sort()

    def test_wobble(self, wbl):
        for i in wbl.range():
            if self[i] not in wbl[i]:
                return False
        return True

    """
    Mutation methods
    ~~~~~~~~~~~~~~~~
    """

    def mutate(self, nt, a, b=None):
        pass
        #TODO: Run random mutations, extract mutations, apply mutations
        #and verify

    def delete(self, a, b=None):
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

        new = self.copy()

        for wbl in new.wobbles:
            if wbl.lt(c): wbl.offset -= 1

        #Deleting insertion
        if len(new.seq[a]) > 1:
            nt = new.seq[a][b][0] 
            del(new.seq[a][b])
            del(new.seql[c])
            #Mutation list change
            del(new.mutations[a][b])
            if b == 0:
                #Convert old insertion to deletion
                if nt != new.seq[a][b][0]:
                    new.mutations[a][b] = "m"
                #Insertion is equal to origital nt
                else:
                    new.mutations[a][b] = "n"
                new.seq[a][b][1] = nt
            return new

        #Prev mut is an insertion
        #Delete insertion and convert current position to a mutation
        if len(new.seq[a-1]) > 0:
            e = len(new.seq[a-1]) - 1
            nt = new.seq[a-1][e][0]
            del(new.seq[a-1][e])
            del(new.seql[c-1])
            del(new.mutations[a-1][e])
            m = "m"
            #Reverted to original
            if nt == new.seq[a][b][1]:
                m = "n"
            new.seq[a][b][0] = nt
            new.mutations[a][b] = m
            return new

        #Normal deletion
        new.seq[a][b][0] = "-"
        del(new.seql[c])
        new.mutations[a] = ["d"]
        return new

    def insert(self, nt, a, b=None):
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

        #Inside wobble, abort
        for wbl in self.wobbles:
            if c in wbl:
                return False

        new = self.copy()

        for wbl in new.wobbles:
            if wbl.lt(c): wbl.offset += 1

        #Inserting into deletion
        if new.seq[a-1][0][0] == "-":
            e = a - 1
            while e != 0 and new.seq[e-1][0][0] == "-":
                e -= 1
            #Convert deletion to insertion
            new.seq[e][0][0] = nt
            if new.seq[e][0][1] != nt:
                new.mutations[e] = ["m"]
            else:
                new.mutations[e] = ["n"]
            new.seql.insert(c, (e, 0))
            return new

        #Inserting into insertions
        if b > 0:
            e = len(new.seq[a])
            new.seq[a].insert(b, [nt, "-"])
            new.mutations[a].insert(b, "i")
            new.seql.insert(c-b+e, (a, e))
            return new

        #No other insertions
        if len(new.seq[a-1]) == 1:
            new.seq[a-1].append([nt, "-"])
            new.mutations[a-1].append("i")
            new.seql.insert(c, (a-1, 1))
            return new

        #Inserting into end of insertions
        e = len(new.seq[a-1])
        new.seq[a-1].append([nt, "-"])
        new.mutations[a-1].append("i")
        new.seql.insert(c, (a-1, b+e))
        return new


    """
    Other
    ~~~~~
    """

    def mutated_str(self):
        return "".join([nt[0] for i, nts in sorted(self.seq.items()) for nt in nts])

    def mutation_str(self):
        l = list()
        for a in sorted(self.seq):
            for b in range(len(self.seq[a])):
                if b == 0 and self.mutations[a][b] == "n":
                    l.append(" ")
                else:
                    l.append(self.mutations[a][b])
        return "".join(l)

    def original_str(self):
        return "".join([nt[1] for i, nts in sorted(self.seq.items()) for nt in nts])

    def __str__(self):
        #return self.org_seq
        return "".join([str(nt) for nt in self])

    def __iter__(self):
        for i, nts in sorted(self.seq.items()):
            for mut, org in nts:
                if mut != "-": yield mut

    def __getitem__(self, a, b=None):
        if b is None:
            a, b = self.seql[a]
        return self.seq[a][b][0]

    def __len__(self):
        return len(self.seql)

    def copy(self):
        return deepcopy(self)

    def pprint(self):
        mut_str = self.mutation_str()
        print(mut_str)
        print(self.original_str())
        print(self.mutated_str())
        gp = 0
        for wbl in self.wobbles:
            inst = sum([n.count("i") for i,n in self.mutations.items() if wbl.lt(i)])
            print(" "*(wbl.start - gp + inst) + str(wbl), end="")
            gp = wbl.end
        print()


class MutationError(ValueError): pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()