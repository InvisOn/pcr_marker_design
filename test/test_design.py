import pytest
from design import *

reference = "/output/genomic/plant/Actinidia/chinensis/CK51F3_01/Genome/Assembly/PS1/1.68.5/AllChromosomes/PS1.1.68.5.fasta"
annotations = "/output/genomic/plant/Actinidia/chinensis/Resequencing/Variants/PS1.1.68.5/52DiploidGenomes/Combined_diploidCK_basic_NS30_Q50_SAFR3_DP50_PAIR0.8_PS1.1.68.5_ann.vcf.gz"
description = "PS1.1.68.5"
targets = "./test/test-data/targets.bed"
    

def test_PrimerDesign_init():
    """
    pass
    """
    designer = PrimerDesign(reference_file=reference,
                            annotation_file=annotations,
                            description=description,
                            targets_file=targets,
                            output_dir=None,
                            amplicon_size_range=(100, 500),
                            primer_size_range=(18, 25),
                            p3_globals=None)
    
    assert designer.amplicon_size_range == (100, 500)
    assert designer.primer_size_range == (18, 25)
    # Needs other asserts
    
    
def test_gettargetdict():
    """
    pass
    """
    designer = PrimerDesign(reference_file=reference,
                            annotation_file=annotations,
                            description=description,
                            targets_file=targets,
                            output_dir=None,
                            amplicon_size_range=(100, 500),
                            primer_size_range=(18, 25),
                            p3_globals=None)
    for targ in d.targets.fetch():
        d = designer.gettargetdict([targ.chrom, targ.start, targ.stop])

# class TestDesign:

#     def test_getseqslicedict(self):
#         """
#         Get dict suitable for primer3
#         """
#         test_seq = "./test/test-data/targets.fasta"
#         annfile = "./test/test-data/targets.snps.bed"
#         designer = d.PrimerDesign(test_seq, annfile, "Designer Test")
#         target = BedTool("k69_93535 1146 1147", from_string=True)
#         max_size = 250
#         target_dic = {'REF_OFFSET': 896,
#  'SEQUENCE_EXCLUDED_REGION': [(244, 1), (439, 1)],
#  'SEQUENCE_ID': 'Designer Test',
#  'SEQUENCE_TARGET': (250, 1),
#  'SEQUENCE_TEMPLATE': 'AAATAATGGAGAATAGATGGTTCAAGAATGGATTCGAGCCTGTGAAATATTACATTGAGAATGATAGGTTTCATAAGTGGTGTAGCTTAGACGAAGAGAATGCTAATGACAACGAGGAGGTAGAATCTGGAGATGAATCAGACTCTTCAGTTGCTTCCTGCCCTCCTACACTTAATGAAGGAAAGAAAAAAAGGACAGGGAAGCTTCATAGGCCTTTGAGTCTGAACGCATTTGACATAATTTCCTTTTCCAGAGGATTTGATCTTTCAGGTTTGTTTGAAGAAACGGGAGATGAAACAAGATTTGTGTCGGGTGAAACGATACCAAACATCATATCGAAATTGGAGGAGATTGCAAAAGTGGGTAGTTTCACGTTTAGGAAGAAGGATTGTAGGGTTAGTTTAGAAGGAACGCGAGAAGGAGTGAAGGGCCCTCTTACGATTGGAGCTGAGATATTTGAGCTTACGCCTAGTTTGGTTGTTGTTGAGCTTAAGAAGAAAG',
#  'TARGET_ID': 'k69_93535:1147-1147'}
#         assert designer.getseqslicedict(target, max_size) == target_dic

#     def test_getVCFseqslicedict(self):
#         """
#         Get dict suitable for primer3 from vcf reader
#         """
#         test_seq = "./test/test-data/AcCHR1_test.fasta"
#         vcffile = "./test/test-data/AcCHR1_test.vcf.gz"
#         designer = d.VcfPrimerDesign(test_seq, vcffile, "TestCHR1")
#         target = Interval('CHR1',3000,3001)
#         max_size = 307
#         target_dic = {'REF_OFFSET': 2693,
#  'SEQUENCE_EXCLUDED_REGION': [(100, 1),
#                               (179, 1),
#                               (214, 1),
#                               (233, 1),
#                               (272, 1),
#                               (300, 1),
#                               (346, 1),
#                               (468, 1),
#                               (529, 1),
#                               (613, 2)],
#  'SEQUENCE_ID': 'CHR1:2693-3308',
#  'SEQUENCE_TARGET': (307, 1),
#  'SEQUENCE_TEMPLATE': 'GGTTGGTCTATTCATCATTGCTCCTAACGCATTCCTCATGGCAATCTGCATTGCTGCCTCAATTTCTTTAGAAGCTTCCAGAGTTGTTGAATTGGCAGCGGCAACTACAGTCGCAACTGTTCCTAGCTTTGCAGAACCATTCCCACTCAAGGAATTCACGGACTCTTTATGTGCCTTCAGAACCAACTGTGTCGCACTGGGTTTTAAAGGAAATAAATAAATATGGAATAAAACATTGATATTACAAATAAAGGGTGCTTCTAGCTGAGTAGTCCTCCGATAAAGCACACGCATACAAAGGAATGAGAGAGAGAGAGAGAGGCGCTACCACATATAAAAGGGACAGCAAACATTTTAACATGAGCAAATCAGTGACACTAGGTAGGTGTTAGCACAAAAATGAACCTTGTTTACATCTGTTCACCACATCCTAGAACATCTTAGACACACACTGCAATAACATATGAGGTGGAGCATGGCACAGTGATACTGCAACAGTAGGATTCCCTGTAACTCTAATGCAACTTTTCATGTACTCAGCCTCTCAAATGATATCGCATGACAAAGTAAAATTAGGTTTTTTTAACTTTTAAACAATAAAACATGAAATTGGAA',
#  'TARGET_ID': 'CHR1:3000-3000'}
#         assert designer.getseqslicedict(target, max_size) == target_dic
