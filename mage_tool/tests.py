#!/usr/bin/env python2

from __future__ import print_function

import sys
import doctest
import unittest

try: import cStringIO as strIO
except ImportError: from io import StringIO as strIO

from Bio import SeqIO

import ViennaRNA
import translation
import IO
import oligo_design
import helpers
import manual

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
     CDS             200..229
                     /gene="fakE"
                     /locus_tag="b0005"
                     /codon_start=1
                     /gene_synonym="fak005"
                     /transl_table=11
     CDS             260..289
                     /gene="fakE"
                     /locus_tag="b0006"
                     /codon_start=1
                     /gene_synonym="fak006"
                     /transl_table=11
     CDS             complement(1053..1073)
                     /gene="fakZ"
                     /locus_tag="b0010"
                     /codon_start=1
                     /gene_synonym="fak010"
                     /transl_table=11
     CDS             complement(321..359)
                     /gene="fakF"
                     /locus_tag="b0007"
                     /codon_start=1
                     /gene_synonym="fak007"
                     /transl_table=11
     CDS             complement(381..419)
                     /gene="fakG"
                     /locus_tag="b0008"
                     /codon_start=1
                     /gene_synonym="fak008"
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
               'TAA': ['$', 0.64],'TAC': ['Y', 0.47],'TAG': ['$', 0.0],
               'TAT': ['Y', 0.53],'TCA': ['S', 0.15],'TCC': ['S', 0.11],
               'TCG': ['S', 0.16],'TCT': ['S', 0.11],'TGA': ['$', 0.36],
               'TGC': ['C', 0.58],'TGG': ['W', 1.0],'TGT': ['C', 0.42],
               'TTA': ['L', 0.15],'TTC': ['F', 0.43],'TTG': ['L', 0.12],
               'TTT': ['F', 0.57]},
'replication': {'ori': [50, 60],
               'ter': [700, 750],
               'ter extended': [600, 800]},
'start_codons': ['ATG', 'GTG', 'TTG', 'ATT', 'CTG'],
'operons': {'3781': {'start': 119,
                     'genes': ['b0002', 'b0003'],
                     'end': 193,
                     'strand': 1},
            '3782': {'start': 320,
                     'genes': ['b0006', 'b0007'],
                     'end': 419,
                     'strand': -1}}}

def testest():
    class Self: pass
    self = Self()
    self.genome = SeqIO.read(strIO.StringIO(genome), "genbank")
    self.config = create_config_tables(config)
    self.genes = seqIO_to_genelist(self.genome, config)


    #seq = self.genes["fakA"].leader
    #args = (self.config["codon_table"], self.config["dgn_table"])
    #No wobble
    #print(find_wobble_seq(self.genome, seq, 10, 13, *args).get_wobble_str())
    #self.assertEqual(find_wobble_seq(self.genome, gene, 10, 13, *args), gene)

    exit(0)



class TestMageTool(unittest.TestCase):
    """Test cases for the various oligo design routines"""
    def setUp(self):
        self.genome = SeqIO.read(strIO.StringIO(genome), "genbank")
        self.config = create_config_tables(config)
        self.genes = seqIO_to_genelist(self.genome, self.config, promoter_len=20)

    def test_seqIO_to_genelist(self):
        #Default instantiation, with everything.
        self.assertItemsEqual(["fakA", "fakB", "fakC", "fakD", "fakZ", "fakE",
                               "b0001", "b0002", "b0003", "b0004", "b0010",
                               "b0006", "b0005", "b0007", "b0008", "fakF",
                               "fakG"], self.genes)

        #Instantiation with a genelist and genome
        genes = seqIO_to_genelist(self.genome, self.config, ("fakA"), True)
        self.assertItemsEqual(genes, ["fakA", "genome"])

        #Instantiation with genes and locus_tags
        genes = seqIO_to_genelist(self.genome, self.config,
                                  ("fakA", "b0001", "b0002"))
        self.assertItemsEqual(["fakA", "b0001", "b0002"], genes)

        #Duplicated gene, throw an exception
        with self.assertRaises(IO.ParserError):
            genes = seqIO_to_genelist(self.genome, self.config, ("fakA", "fakE"))

    def test_find_wobble(self):
        args = (self.config["codon_table"], self.config["dgn_table"])

        #No wobble
        seq = oligo_design.Sequence("TCT")
        w = find_wobble_seq(self.genome, seq, 10, 13, *args)
        self.assertEqual(w.get_wobble_str(), "NNN")

        #Right-side wobble
        seq = oligo_design.Sequence("TCTGACTGC")
        w = find_wobble_seq(self.genome, seq, 10, 19, *args)
        self.assertEqual(w.get_wobble_str(), "NNNGAYTGY")

        #Left-side wobble
        seq = oligo_design.Sequence("TGTCTCTGTG")
        w = find_wobble_seq(self.genome, seq, 30, 40, *args)
        self.assertEqual(w.get_wobble_str(), "YGTNWSNNNN")

        #Both-side wobble
        seq = oligo_design.Sequence("TCTGACTGCAACGGGCAATATGTCTCTGTG")
        w = find_wobble_seq(self.genome, seq, 10, 40, *args)
        self.assertEqual(w.get_wobble_str(), "NNNGAYTGYAAYGGNCARTAYGTNWSNNNN")

        #Inside wobble
        seq = oligo_design.Sequence("GAC")
        w = find_wobble_seq(self.genome, seq, 13, 16, *args)

        self.assertEqual(w.get_wobble_str(), "GAY")
        #Right-side wobble (complement gene)
        seq = oligo_design.Sequence("AAAGAGTGTC")
        w = find_wobble_seq(self.genome, seq, 50, 60, *args)
        self.assertEqual(w.get_wobble_str(), "NNNNARNGTY")

        #Left side wobble (complement gene)
        seq = oligo_design.Sequence("GAACTG")
        w = find_wobble_seq(self.genome, seq, 74, 80, *args)
        self.assertEqual(w.get_wobble_str(), "RAANNN")

        #TODO: Test wobble exceeding genome

    def test_sequences(self):
        self.assertEqual(str(self.genes["fakA"].cds), "GACTGCAACGGGCAATATGTCTCT")
        self.assertEqual(str(self.genes["fakD"].cds), "TTCAGAAGCTGCTATCAGACACTC")
        self.assertEqual(self.genes["fakA"].pos, 13)
        self.assertEqual(self.genes["fakD"].pos, 77)

    def test_leaders(self):
        #Standard leader
        self.assertEqual(str(self.genes["fakB"].leader), "ACCTGCCGTGAGTAAATTAAAATTTTATTGACTTA")
        #Extends beyond 0:
        self.assertEqual(str(self.genes["fakA"].leader), "CGATGCGAGGTTGTTGAAGTCGAGCTTTTCATTCT")
        #Reverse leader
        self.assertEqual(str(self.genes["fakD"].leader), "AATAAAATTTTAATTTACTCACGGCAGGTAACCAG")
        #Reverse leader extends beyond length
        self.assertEqual(str(self.genes["fakZ"].leader), "TTGCCCGTTGCAGTCAGAATGAAAAGCTCGACTTC")

    def test_wobbles(self):
        #No wobble
        self.assertEqual(self.genes["fakB"].leader.get_wobble_str(), "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
        self.assertEqual(self.genes["fakA"].leader.get_wobble_str(), "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN")
        #Wobble on one side
        self.assertEqual(self.genes["fakC"].leader.get_wobble_str(), "GNCAYTRRATHYTNTRRCCNATHTRRGCNNNNNNN")
        #print str(self.genes["fakC"].leader_wobble)
        #TODO:
        # test wobble on right side, both sides, and everywhere
        # test same on -1 strand
        # test wobbles extending beyond limits
        # test multiple wobble

    def test_promoter_seq(self):
        #Operons, +1 and -1
        self.assertEqual(str(self.genes["fakB"].promoter), "ATTAAAATTTTATTGACTTA")
        self.assertEqual(str(self.genes["fakC"].promoter), "ATTAAAATTTTATTGACTTA")
        self.assertEqual(self.genes["fakB"].promoter_pos, -20)
        self.assertEqual(self.genes["fakC"].promoter_pos, -56)
        self.assertEqual(str(self.genes["fakG"].promoter), "TGGCCACCTGCCCCTGCCTG")
        self.assertEqual(str(self.genes["fakF"].promoter), "TGGCCACCTGCCCCTGCCTG")
        self.assertEqual(self.genes["fakG"].promoter_pos, -20)
        self.assertEqual(self.genes["fakF"].promoter_pos, -80)
        #Not in operon, +1 and -1
        self.assertEqual(str(self.genes["fakA"].promoter), "GAAGTCGAGCTTTTCATTCT")
        self.assertEqual(self.genes["fakA"].promoter_pos, -20)
        self.assertEqual(str(self.genes["fakD"].promoter), "TACTCACGGCAGGTAACCAG")
        self.assertEqual(self.genes["fakD"].promoter_pos, -20)

    def test_do_mutation(self):
        mut1 = oligo_design.Mutation("TG", "GT", 3)
        mut1 = self.genes["fakA"].do_mutation(mut1)
        self.assertEqual(str(mut1), "[TG=gt].16")

        mut2 = oligo_design.Mutation("CA", "GG", 2)
        mut2 = self.genes["fakD"].do_mutation(mut2)
        self.assertEqual(str(mut2), "[TG=cc].73")


    def test_KO(self):
        mut1 = translation.translational_KO(self.genes["fakA"])
        self.assertEqual(str(mut1), "[CAACGG=atAata].18")
        mut2 = translation.translational_KO(self.genes["fakC"])
        self.assertEqual(str(mut2), "[GACA=tAgt].157")
        mut3 = translation.translational_KO(self.genes["fakC"],["TAG", "TAA", "TGA"], 4)
        self.assertEqual(str(mut3), "[GACAGATAAA=tAgtGATAAt].157")

    def test_a_lot_of_sequence_mutations(self):
        """Shotgun a lot of mutations and do some automated tests."""
        s = oligo_design.Sequence("ATAGCACATAAGACTAGACAT")

        for i in range(500):
            s = s.random_mutation(5)
            sseql = [(a, b) for a in range(s.org_len)
                            for b in range(len(s.seq[a]))
                            if s.mutations[a][b] != "d"]
            if sseql != s.seql:
                raise Exception("aoe")




if __name__ == "__main__":
    if "t" in sys.argv:
        testest()

    #suite = unittest.TestLoader().discover("./", pattern="*_tests.py")
    suite = unittest.TestLoader().loadTestsFromTestCase(TestMageTool)

    suite.addTests(doctest.DocTestSuite(ViennaRNA))
    suite.addTests(doctest.DocTestSuite(translation))
    suite.addTests(doctest.DocTestSuite(IO))
    suite.addTests(doctest.DocTestSuite(oligo_design))
    suite.addTests(doctest.DocTestSuite(helpers))
    suite.addTests(doctest.DocTestSuite(manual))

    unittest.TextTestRunner(verbosity=2).run(suite)
