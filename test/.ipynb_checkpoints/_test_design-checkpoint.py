import pytest
from pcr_marker_design import design as d
from pybedtools import BedTool, Interval
from pyfaidx import Fasta


class TestDesign:

    def test_getseqslicedict(self):
        """
        Get dict suitable for primer3
        """
        test_seq = "./test/test-data/targets.fasta"
        annfile = "./test/test-data/targets.snps.bed"
        designer = d.PrimerDesign(test_seq, annfile, "Designer Test")
        target = BedTool("k69_93535 1146 1147", from_string=True)
        max_size = 250
        target_dic = {'REF_OFFSET': 896,
 'SEQUENCE_EXCLUDED_REGION': [(244, 1), (439, 1)],
 'SEQUENCE_ID': 'Designer Test',
 'SEQUENCE_TARGET': (250, 1),
 'SEQUENCE_TEMPLATE': 'AAATAATGGAGAATAGATGGTTCAAGAATGGATTCGAGCCTGTGAAATATTACATTGAGAATGATAGGTTTCATAAGTGGTGTAGCTTAGACGAAGAGAATGCTAATGACAACGAGGAGGTAGAATCTGGAGATGAATCAGACTCTTCAGTTGCTTCCTGCCCTCCTACACTTAATGAAGGAAAGAAAAAAAGGACAGGGAAGCTTCATAGGCCTTTGAGTCTGAACGCATTTGACATAATTTCCTTTTCCAGAGGATTTGATCTTTCAGGTTTGTTTGAAGAAACGGGAGATGAAACAAGATTTGTGTCGGGTGAAACGATACCAAACATCATATCGAAATTGGAGGAGATTGCAAAAGTGGGTAGTTTCACGTTTAGGAAGAAGGATTGTAGGGTTAGTTTAGAAGGAACGCGAGAAGGAGTGAAGGGCCCTCTTACGATTGGAGCTGAGATATTTGAGCTTACGCCTAGTTTGGTTGTTGTTGAGCTTAAGAAGAAAG',
 'TARGET_ID': 'k69_93535:1147-1147'}
        assert designer.getseqslicedict(target, max_size) == target_dic

    def test_getVCFseqslicedict(self):
        """
        Get dict suitable for primer3 from vcf reader
        """
        test_seq = "./test/test-data/AcCHR1_test.fasta"
        vcffile = "./test/test-data/AcCHR1_test.vcf.gz"
        designer = d.VcfPrimerDesign(test_seq, vcffile, "TestCHR1")
        target = Interval('CHR1',3000,3001)
        max_size = 307
        target_dic = {'REF_OFFSET': 2693,
 'SEQUENCE_EXCLUDED_REGION': [(100, 1),
                              (179, 1),
                              (214, 1),
                              (233, 1),
                              (272, 1),
                              (300, 1),
                              (346, 1),
                              (468, 1),
                              (529, 1),
                              (613, 2)],
 'SEQUENCE_ID': 'CHR1:2693-3308',
 'SEQUENCE_TARGET': (307, 1),
 'SEQUENCE_TEMPLATE': 'GGTTGGTCTATTCATCATTGCTCCTAACGCATTCCTCATGGCAATCTGCATTGCTGCCTCAATTTCTTTAGAAGCTTCCAGAGTTGTTGAATTGGCAGCGGCAACTACAGTCGCAACTGTTCCTAGCTTTGCAGAACCATTCCCACTCAAGGAATTCACGGACTCTTTATGTGCCTTCAGAACCAACTGTGTCGCACTGGGTTTTAAAGGAAATAAATAAATATGGAATAAAACATTGATATTACAAATAAAGGGTGCTTCTAGCTGAGTAGTCCTCCGATAAAGCACACGCATACAAAGGAATGAGAGAGAGAGAGAGAGGCGCTACCACATATAAAAGGGACAGCAAACATTTTAACATGAGCAAATCAGTGACACTAGGTAGGTGTTAGCACAAAAATGAACCTTGTTTACATCTGTTCACCACATCCTAGAACATCTTAGACACACACTGCAATAACATATGAGGTGGAGCATGGCACAGTGATACTGCAACAGTAGGATTCCCTGTAACTCTAATGCAACTTTTCATGTACTCAGCCTCTCAAATGATATCGCATGACAAAGTAAAATTAGGTTTTTTTAACTTTTAAACAATAAAACATGAAATTGGAA',
 'TARGET_ID': 'CHR1:3000-3000'}
        assert designer.getseqslicedict(target, max_size) == target_dic
