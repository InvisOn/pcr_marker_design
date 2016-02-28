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
        test_seq="./test/test-data/Kiwifruit_pseudomolecule_Chr.fa"
        vcffile= "./test/test-data/Chr9_Myb210.vcf.gz"
        designer = d.VcfPrimerDesign(test_seq,vcffile,"MybTest")
        target=BedTool("Chr9 1390065 1390066",from_string=True)
        max_size=100
        target_dic={'SEQUENCE_EXCLUDED_REGION': [(82, 4),(92, 1),(113, 1),(126, 1),(127, 1),(139, 3),
        (150, 1),(164, 2),(192, 1)],
        'SEQUENCE_ID': 'MybTest',
        'SEQUENCE_TARGET': (100, 1),
        'SEQUENCE_TEMPLATE': 'TAAATTTATAAAAAATATATAAATATATATATATATAAATATATAAATATATATATATATAAATATATAAATATATATATATATATATATATATATATATATGTGAGAAAGTATTGATTAACCATAGCGGAGTAAACAACAGAGTTAGAAAAACTCTTTAAAAAAATGTGCTCTCTTTTGTGTGGGGAAAAAGCGCCACAG',
        'TARGET_ID': 'Chr9_1389965_1390166'}
        assert designer.getseqslicedict(target,max_size) == target_dic
