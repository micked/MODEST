import unittest

import random

from .. import oligo_design
from ..IO import Mutation

ref_genome = """GACTAGGGTCCCCGCTCATAAACGATATTAGTTGCAGTCACTAAGGTGCTAATTTTGCTATATT
                GCAAAGGCGAAAACATTTAGTAATCCGCGAGAGGCATTCTGAAAGTTCCTGAATGGTAATCCGC
                CAAACACTGTGTTCACTGTTCACAAAAGTACAGAAAGTCATGTGACGAAAACGTATCGGTCTTG
                GTTAGTAGGTACAACCACGTAAGTAAAAACCTCCGTCGAATTTTTTATAATCTATTAGACA"""
ref_genome = "".join(ref_genome.split())


class TestOligoDesign(unittest.TestCase):
    """Test cases for the various oligo design routines"""
    def setUp(self):
        pass

    def test_deletion_to_oligo(self):
        """Test a deletion oligo"""
        deletion = Mutation("eq", "AAA=", 67)
        oligo = oligo_design.mut_to_oligo(deletion, ref_genome, 90)
        target = "ACGATATTAGTTGCAGTCACTAAGGTGCTAATTTTGCTATATTGCGGCGAAAACATTTAGTAATCCGCGAGAGGCATTCTGAAAGTTCCT"
        self.assertEqual(str(oligo), target)
        
    def test_insertion_to_olig_min(self):
        """Test an insertion oligo (minimal mut encoding)
        
        TODO: get notation straight
        """
        insertion = Mutation("eq", "=TAG", 107)
        oligo = oligo_design.mut_to_oligo(insertion, ref_genome, 90)
        target = "TTGCAAAGGCGAAAACATTTAGTAATCCGCGAGAGGCATTCTGAAtagAGTTCCTGAATGGTAATCCGCCAAACACTGTGTTCACTGTTC"
        target = "TTGCAAAGGCGAAAACATTTAGTAATCCGCGAGAGGCATTCTGAtagAAGTTCCTGAATGGTAATCCGCCAAACACTGTGTTCACTGTTC"
        self.assertEqual(str(oligo), target)
        
    def test_insertion_to_olig_extr(self):
        """Test an insertion oligo (including org seq)"""
        insertion = Mutation("eq", "A=ATAG", 107)
        oligo = oligo_design.mut_to_oligo(insertion, ref_genome, 90)
        target = "TGCAAAGGCGAAAACATTTAGTAATCCGCGAGAGGCATTCTGAatagAGTTCCTGAATGGTAATCCGCCAAACACTGTGTTCACTGTTCA"
        self.assertEqual(str(oligo), target)
        
    def test_mut_to_olig(self):
        """Test a mutation oligo (same length)"""
        mut = Mutation("eq", "AGA=GTG", 160)
        oligo = oligo_design.mut_to_oligo(mut, ref_genome, 90)
        target = "AATGGTAATCCGCCAAACACTGTGTTCACTGTTCACAAAAGTACgtgAAGTCATGTGACGAAAACGTATCGGTCTTGGTTAGTAGGTACA"
        self.assertEqual(str(oligo), target)
        
    def test_mut_to_olig_single(self):
        """Test a mutation oligo (arrow notation)"""
        mut = Mutation("arrow", "G->T", 73, mut_type="point_mutation")
        oligo = oligo_design.mut_to_oligo(mut, ref_genome, 90)
        target = "TTAGTTGCAGTCACTAAGGTGCTAATTTTGCTATATTGCAAAGGCtAAAACATTTAGTAATCCGCGAGAGGCATTCTGAAAGTTCCTGAA"
        self.assertEqual(str(oligo), target)
