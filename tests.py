#!/usr/bin/env python

import unittest

from mage_tool.tests import oligo_design_tests
from mage_tool.IO import Mutation
from mage_tool import oligo_design

ref_genome = """GACTAGGGTCCCCGCTCATAAACGATATTAGTTGCAGTCACTAAGGTGCTAATTTTGCTATATT
                GCAAAGGCGAAAACATTTAGTAATCCGCGAGAGGCATTCTGAAAGTTCCTGAATGGTAATCCGC
                CAAACACTGTGTTCACTGTTCACAAAAGTACAGAAAGTCATGTGACGAAAACGTATCGGTCTTG
                GTTAGTAGGTACAACCACGTAAGTAAAAACCTCCGTCGAATTTTTTATAATCTATTAGACA"""
ref_genome = "".join(ref_genome.split())

if __name__ == "__main__":
    suite = unittest.TestLoader().discover("./", pattern="*_tests.py")
    #suite = unittest.TestLoader().loadTestsFromTestCase(oligo_design_tests.TestOligoDesign)
    unittest.TextTestRunner(verbosity=2).run(suite)
