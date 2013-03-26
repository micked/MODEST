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
        self.project = project
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
            extended = genome[self.mut.pos-5:self.mut.pos+len(self.mut.before)+5]
            raise Exception("Trying to mutate {}, but found {} in genome. [{}]".format(self.mut.before, found, extended))

        #Calculate flanking sequence lengths
        post_seq_len = (self.oligo_len - len(self.mut.after))/2
        pre_seq_len = self.oligo_len - post_seq_len - len(self.mut.after) + offset
        post_seq_len -= offset

        #Fetch pre sequence
        pre_seq = genome[self.mut.pos-pre_seq_len : self.mut.pos]

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
            log.warning("Oligo: {} inside oriC".format(self.id()))
            return False
        elif ter[0] <= self.pos <= ter[1]:
            self.replichore = 0
            log.warning("Oligo: {} inside termination region".format(self.id()))
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
    
    #OLD
    #def add_barcode(self, barcode, barcoding_lib):
    #    """TODO"""
    #    #Add a barcode ID to self.barcode_ids
    #    #and add barcode sequences to self.barcodes_*
    #    #forward should be prepended, and backward should be appended
    #    #Maybe add some primer checking too (maybe second function)
    #    self.barcodes_forward.insert(0, barcoding_lib[barcode]["forward"].lower())
    #    self.barcodes_reverse.append(barcoding_lib[barcode]["reverse"].lower())
    #    self.barcode_ids.insert(0, barcode)
        
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
        return "{prj}{i:0>4}{barcodes}_{mut}_{gene}_{code}_RP{rp}{opt}".format(
            prj = self.project if self.project else "",
            i = self.number,
            barcodes = "_BC:" + "+".join(self.barcode_ids) if self.barcode_ids else "",
            mut = self.mut.small_str(),
            gene = str(self.gene),
            code = self.code,
            rp = self.replichore,
            opt = "_OPT:{}({:.1f})".format(str(self.optimised), self.dG_fold) if self.optimised else ""
            )

    def short_id(self):
        return "{prj}{i:0>4}".format(
            prj = self.project if self.project else "",
            i = self.number
            )

    def copy(self):
        """Return a copy"""
        return deepcopy(self)


class Mutation:
    def __init__(self, mut_format, mut, pos=0, mut_type="", ref_genome=False):
        """Mutation
        
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
    
    #
    # Printing
    #
    
    def __repr__(self):
        return "Mutation: [{}={}] at pos {}".format(
            self.before,
            self.after,
            self.pos+1)

    def __str__(self):
        return "[{}={}].{}".format(self.before, self.after, self.pos+1)
    
    def small_str(self):
        before = self.before
        after = self.after
        if len(before) > 6:
            before = before[0:2] + ".." + before[-2:]
        if len(after) > 6:
            after = after[0:2] + ".." + after[-2:]

        return "[{}={}].{}".format(before, after, self.pos+1)

    #
    # Parsers
    #
    
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
    """Defines a single gene and all relevant information"""

    def __init__(self, name, pos, strand, cds, leader, leader_wobble=None, promoter=None, promoter_pos=None):
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
        if self.strand == 1:
            mut.pos = mut.pos+self.pos
            return mut
        else:
            mut.pos = self.pos - mut.pos - len(mut.before)
            mut.before = reverse_complement(mut.before)
            mut.after = reverse_complement(mut.after)

        return mut

    def copy(self):
        """Return a copy"""
        return deepcopy(self)
