import pytest
from pcr_marker_design import design as d
from pybedtools import BedTool
from pyfaidx import Fasta


class TestDesign:

    def test_getseqslicedict(self):
        """
        Get dict suitable for primer3
        """
        test_seq="./test/test-data/targets.fasta"
        annfile= "./test/test-data/targets.gff"
        designer = d.PrimerDesign(test_seq,annfile,"targets")
        target=BedTool("k69_93535 1146 1147",from_string=True)
        max_size=250
        target_dic={'SEQUENCE_EXCLUDED_REGION': [(244, 2), (439, 2)],
        'SEQUENCE_ID': 'targets',
        'SEQUENCE_TARGET': (250, 1),
        'SEQUENCE_TEMPLATE': 'AAATAATGGAGAATAGATGGTTCAAGAATGGATTCGAGCCTGTGAAATATTACATTGAGAATGATAGGTTTCATAAGTGGTGTAGCTTAGACGAAGAGAATGCTAATGACAACGAGGAGGTAGAATCTGGAGATGAATCAGACTCTTCAGTTGCTTCCTGCCCTCCTACACTTAATGAAGGAAAGAAAAAAAGGACAGGGAAGCTTCATAGGCCTTTGAGTCTGAACGCATTTGACATAATTTCCTTTTCCAGAGGATTTGATCTTTCAGGTTTGTTTGAAGAAACGGGAGATGAAACAAGATTTGTGTCGGGTGAAACGATACCAAACATCATATCGAAATTGGAGGAGATTGCAAAAGTGGGTAGTTTCACGTTTAGGAAGAAGGATTGTAGGGTTAGTTTAGAAGGAACGCGAGAAGGAGTGAAGGGCCCTCTTACGATTGGAGCTGAGATATTTGAGCTTACGCCTAGTTTGGTTGTTGTTGAGCTTAAGAAGAAAG',
        'TARGET_ID': 'k69_93535_896_1397'}
        assert designer.getseqslicedict(target,max_size) == target_dic

    def test_getVCFseqslicedict(self):
        """
        Get dict suitable for primer3 from vcf reader
            """
        test_seq="./test/test-data/AcCHR1_test.fasta"
        vcffile= "./test/test-data/AcCHR1_test.vcf.gz"
        designer = d.VcfPrimerDesign(test_seq,vcffile,"TestCHR1")
        target=BedTool('CHR1 3000 3001',from_string=True)
        max_size=100
        target_dic={'SEQUENCE_EXCLUDED_REGION': [(7, 1), (26, 1), (65, 1), (93, 1), (139, 1)],
        'SEQUENCE_ID': 'TestCHR1',
        'SEQUENCE_TARGET': (100, 1),
        'SEQUENCE_TEMPLATE': 'AGGAAATAAATAAATATGGAATAAAACATTGATATTACAAATAAAGGGTGCTTCTAGCTGAGTAGTCCTCCGATAAAGCACACGCATACAAAGGAATGAGAGAGAGAGAGAGAGGCGCTACCACATATAAAAGGGACAGCAAACATTTTAACATGAGCAAATCAGTGACACTAGGTAGGTGTTAGCACAAAAATGAACCTT',
        'TARGET_ID': 'CHR1_2900_3101'}
        assert designer.getseqslicedict(target,max_size) == target_dic
