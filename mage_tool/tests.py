#!/usr/bin/env python

import unittest
import doctest
import cStringIO
import sys

from Bio import SeqIO

import ViennaRNA
import translation
import IO
import oligo_design

from IO import seqIO_to_genelist
from IO import create_config_tables
from IO import find_wobble_seq

genome = """LOCUS       FK_000001               1080 bp    DNA              BCT 11-JAN-2012
DEFINITION  Escherichia coli str. K-12 substr. MG1655 chromosome, complete
            genome.
ACCESSION   FK_001
VERSION     FK_001.1
DBLINK      Project:57779 BioProject:PRJNA57779
KEYWORDS    .
SOURCE      .
  ORGANISM  .
            .
COMMENT     Completely fake genome for testing purposes
FEATURES             Location/Qualifiers
     CDS             14..37
                     /gene="fakA"
                     /locus_tag="b0001"
                     /codon_start=1
                     /gene_synonym="fak001"
                     /transl_table=11
     CDS             complement(54..77)
                     /gene="fakD"
                     /locus_tag="b0004"
                     /codon_start=1
                     /gene_synonym="fak004"
                     /transl_table=11
     CDS             119..148
                     /gene="fakB"
                     /locus_tag="b0002"
                     /codon_start=1
                     /gene_synonym="fak002"
                     /transl_table=11
     CDS             155..193
                     /gene="fakC"
                     /locus_tag="b0003"
                     /codon_start=1
                     /gene_synonym="fak003"
                     /transl_table=11
ORIGIN
        1 agcttttcat tctgactgca acgggcaata tgtctctgtg tggattaaaa aaagagtgtc
       61 tgatagcagc ttctgaactg gttacctgcc gtgagtaaat taaaatttta ttgacttagg
      121 tcactaaata ctttaaccaa tataggcata gcgcacagac agataaaaat tacagagtac
      181 acaacatcca tgaaacgcat tagcaccacc attaccacca ccatcaccat taccacaggt
      241 aacggtgcgg gctgacgcgt acaggaaaca cagaaaaaag cccgcacctg acagtgcggg
      301 cttttttttt cgaccaaagg taacgaggta acaaccatgc gagtgttgaa gttcggcggt
      361 acatcagtgg caaatgcaga acgttttctg cgtgttgccg atattctgga aagcaatgcc
      421 aggcaggggc aggtggccac cgtcctctct gcccccgcca aaatcaccaa ccacctggtg
      481 gcgatgattg aaaaaaccat tagcggccag gatgctttac ccaatatcag cgatgccgaa
      541 cgtatttttg ccgaactttt gacgggactc gccgccgccc agccggggtt cccgctggcg
      601 caattgaaaa ctttcgtcga tcaggaattt gcccaaataa aacatgtcct gcatggcatt
      661 agtttgttgg ggcagtgccc ggatagcatc aacgctgcgc tgatttgccg tggcgagaaa
      721 atgtcgatcg ccattatggc cggcgtatta gaagcgcgcg gtcacaacgt tactgttatc
      781 gatccggtcg aaaaactgct ggcagtgggg cattacctcg aatctaccgt cgatattgct
      841 gagtccaccc gccgtattgc ggcaagccgc attccggctg atcacatggt gctgatggca
      901 ggtttcaccg ccggtaatga aaaaggcgaa ctggtggtgc ttggacgcaa cggttccgac
      961 tactctgctg cggtgctggc tgcctgttta cgcgccgatt gttgcgagat ttggacggac
     1021 gttgacgggg tctatacctg cgacccgcgt caggtgcccg atgcgaggtt gttgaagtcg
//
"""

config = {'Definition': 'Escherichia coli str. K-12 substr. MG1655',
'Locus': 'FK_000001',
'codon_usage': {'AAA': ['K', 0.73],'AAC': ['N', 0.53],'AAG': ['K', 0.27],
               'AAT': ['N', 0.47],'ACA': ['T', 0.13],'ACC': ['T', 0.47],
               'ACG': ['T', 0.24],'ACT': ['T', 0.16],'AGA': ['R', 0.02],
               'AGC': ['S', 0.33],'AGG': ['R', 0.03],'AGT': ['S', 0.14],
               'ATA': ['I', 0.07],'ATC': ['I', 0.35],'ATG': ['M', 1.0],
               'ATT': ['I', 0.58],'CAA': ['Q', 0.3],'CAC': ['H', 0.45],
               'CAG': ['Q', 0.7],'CAT': ['H', 0.55],'CCA': ['P', 0.14],
               'CCC': ['P', 0.13],'CCG': ['P', 0.55],'CCT': ['P', 0.17],
               'CGA': ['R', 0.07],'CGC': ['R', 0.44],'CGG': ['R', 0.07],
               'CGT': ['R', 0.36],'CTA': ['L', 0.05],'CTC': ['L', 0.1],
               'CTG': ['L', 0.46],'CTT': ['L', 0.12],'GAA': ['E', 0.7],
               'GAC': ['D', 0.35],'GAG': ['E', 0.3],'GAT': ['D', 0.65],
               'GCA': ['A', 0.21],'GCC': ['A', 0.31],'GCG': ['A', 0.38],
               'GCT': ['A', 0.11],'GGA': ['G', 0.13],'GGC': ['G', 0.46],
               'GGG': ['G', 0.12],'GGT': ['G', 0.29],'GTA': ['V', 0.17],
               'GTC': ['V', 0.18],'GTG': ['V', 0.4],'GTT': ['V', 0.25],
               'TAA': ['*', 0.64],'TAC': ['Y', 0.47],'TAG': ['*', 0.0],
               'TAT': ['Y', 0.53],'TCA': ['S', 0.15],'TCC': ['S', 0.11],
               'TCG': ['S', 0.16],'TCT': ['S', 0.11],'TGA': ['*', 0.36],
               'TGC': ['C', 0.58],'TGG': ['W', 1.0],'TGT': ['C', 0.42],
               'TTA': ['L', 0.15],'TTC': ['F', 0.43],'TTG': ['L', 0.12],
               'TTT': ['F', 0.57]},
'replication': {'ori': [50, 60],
               'ter': [700, 750],
               'ter extended': [600, 800]},
'start_codons': ['ATG', 'GTG', 'TTG', 'ATT', 'CTG']}


def testest():
    selfgenome = SeqIO.read(cStringIO.StringIO(genome), "genbank")
    selfconfig = create_config_tables(config)
    selfgenes = seqIO_to_genelist(selfgenome, config)

    exit(0)



class TestMageTool(unittest.TestCase):
    """Test cases for the various oligo design routines"""
    def setUp(self):
        self.genome = SeqIO.read(cStringIO.StringIO(genome), "genbank")
        self.config = create_config_tables(config)
        self.genes = seqIO_to_genelist(self.genome, self.config)

    def test_find_wobble(self):
        #No wobble
        self.assertEqual(find_wobble_seq(self.genome, 10, 13, self.config["codon_table"], self.config["dgn_table"]), None)
        #Right-side wobble
        self.assertEqual(find_wobble_seq(self.genome, 10, 19, self.config["codon_table"], self.config["dgn_table"]), "NNNGAYTGY")
        #Left-side wobble
        self.assertEqual(find_wobble_seq(self.genome, 30, 40, self.config["codon_table"], self.config["dgn_table"]), "YGTNWSNNNN")
        #Both-side wobble
        self.assertEqual(find_wobble_seq(self.genome, 10, 40, self.config["codon_table"], self.config["dgn_table"]), "NNNGAYTGYAAYGGNCARTAYGTNWSNNNN")
        #Inside wobble
        self.assertEqual(find_wobble_seq(self.genome, 13, 16, self.config["codon_table"], self.config["dgn_table"]), "GAY")
        #Right-side wobble (complement gene)
        self.assertEqual(find_wobble_seq(self.genome, 50, 60, self.config["codon_table"], self.config["dgn_table"]), "NNNNARNGTY")
        #Left side wobble (complement gene)
        self.assertEqual(find_wobble_seq(self.genome, 74, 80, self.config["codon_table"], self.config["dgn_table"]), "RAANNN")
        #TODO: Test wobble exceeding genome

    def test_sequences(self):
        self.assertEqual(str(self.genes["fakA"].cds), "GACTGCAACGGGCAATATGTCTCT")
        self.assertEqual(str(self.genes["fakD"].cds), "TTCAGAAGCTGCTATCAGACACTC")

    def test_leaders(self):
        #Standard leader
        self.assertEqual(str(self.genes["fakB"].leader), "ACCTGCCGTGAGTAAATTAAAATTTTATTGACTTA")
        #Extends beyond 0:
        self.assertEqual(str(self.genes["fakA"].leader), "CGATGCGAGGTTGTTGAAGTCGAGCTTTTCATTCT")
        #TODO: reverse leaders

    def test_wobbles(self):
        #No wobble
        self.assertEqual(self.genes["fakB"].leader_wobble, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
        self.assertEqual(self.genes["fakA"].leader_wobble, "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
        #Wobble on one side
        self.assertEqual(self.genes["fakC"].leader_wobble, "GNCAYTRRATHYTNTRRCCNATHTRRGCNNNNNNN")
        #print str(self.genes["fakC"].leader_wobble)
        #TODO:
        # test wobble on right side, both sides, and everywhere
        # test same on -1 strand


if __name__ == "__main__":
    if "t" in sys.argv:
        testest()

    #suite = unittest.TestLoader().discover("./", pattern="*_tests.py")
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMageTool)

    suite.addTests(doctest.DocTestSuite(ViennaRNA))
    suite.addTests(doctest.DocTestSuite(translation))
    suite.addTests(doctest.DocTestSuite(IO))
    suite.addTests(doctest.DocTestSuite(oligo_design))

    unittest.TextTestRunner(verbosity=2).run(suite)
