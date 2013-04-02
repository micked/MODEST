#!/usr/bin/env python

import unittest
import doctest

import ViennaRNA
import translation
import IO
import oligo_design

class TestMageTool(unittest.TestCase):
    """Test cases for the various oligo design routines"""
    def setUp(self):
        pass

    def test_nothing(self):
        self.assertTrue(True)

if __name__ == "__main__":
    #suite = unittest.TestLoader().discover("./", pattern="*_tests.py")
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMageTool)

    suite.addTests(doctest.DocTestSuite(ViennaRNA))
    suite.addTests(doctest.DocTestSuite(translation))
    suite.addTests(doctest.DocTestSuite(IO))
    suite.addTests(doctest.DocTestSuite(oligo_design))

    unittest.TextTestRunner(verbosity=2).run(suite)
