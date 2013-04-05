#!/usr/bin/env python

"""
Module for designing DNA oligos from mutation objects
"""

import logging
from copy import deepcopy

from helpers import reverse_complement
import ViennaRNA

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


class MutationError(ValueError): pass


if __name__ == "__main__":
    import doctest
    doctest.testmod()