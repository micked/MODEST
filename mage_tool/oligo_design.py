#!/usr/bin/env python

"""
Module for designing DNA oligos from mutation objects
"""

import logging
from copy import deepcopy

from helpers import reverse_complement

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
        self.replichore = None
        self.project = project
        self.number = number
        self.barcodes_forward = []
        self.barcodes_reverse = []
        self.barcode_ids = []

    def make_oligo(self, genome):
        """Make oligo from mutation"""
        #Make sure what is being mutated is actually being mutated
        if str(genome[self.mut.pos:self.mut.pos+len(self.mut.before)]) != str(self.mut.before):
            found = genome[self.mut.pos:self.mut.pos+len(self.mut.before)]
            raise Exception("Trying to mutate {}, but found {} in genome.".format(self.mut.before, found))

        #Calculate flanking sequence lengths
        post_seq_len = (self.oligo_len - len(self.mut.after))/2
        pre_seq_len = self.oligo_len - post_seq_len - len(self.mut.after)
        
        #Fetch pre sequence
        pre_seq = genome[self.mut.pos-pre_seq_len : self.mut.pos]

        #Fetch post sequence
        post_seq_start = self.mut.pos+len(self.mut.before)
        post_seq = genome[post_seq_start:post_seq_start+post_seq_len]
        
        #Set parameters
        self.oligo = pre_seq + self.mut.after.lower() + post_seq

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

    def add_barcode(self, barcode, barcoding_lib):
        """TODO"""
        #Add a barcode ID to self.barcode_ids
        #and add barcode sequences to self.barcodes_*
        #forward should be prepended, and backward should be appended
        #Maybe add some primer checking too (maybe second function)
        pass

    def output(self):
        """Return the oligo with barcodes"""
        out = self.barcodes_forward + [self.oligo] + self.barcodes_reverse
        out = map(str, out)
        return "".join(out)

    def id(self):
        """Calculate an ID based the various self-contained variables"""
        return "{prj}{i:0>4}{barcodes}_{mut}_{gene}_RP{rp}{opt}".format(
            prj = self.project if self.project else "",
            i = self.number,
            barcodes = "_BC:" + ",".join(self.barcode_ids) if self.barcode_ids else "",
            mut = str(self.mut),
            gene = str(self.gene),
            rp = self.replichore,
            opt = "OPT:" + str(self.optimised) if self.optimised else ""
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
    
    #
    # Printing
    #
    
    def __repr__(self):
        return "Mutation: [{}={}] at pos {}".format(
            self.before,
            self.after,
            self.pos)

    def __str__(self):
        return "[{}={}].{}".format(
            self.before,
            self.after,
            self.pos)
    
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
    def __init__(self, name, pos, strand, cds, leader, promoter = None, promoter_pos = None):
        self.name = name
        self.pos = pos
        self.strand = strand
        self.cds = cds
        self.leader = leader
        self.leader_pos = pos-len(leader)
        self.promoter = promoter
        self.promoter_pos = promoter_pos

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