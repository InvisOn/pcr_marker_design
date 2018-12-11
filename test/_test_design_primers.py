import os
from design_primers import parse_args, design_primers


def test_design_primers():
    """
    Test for function design_primers.
    Uses argument parser to ensure input is in expected form.
    This is effectively an full integration test.
    """
    # create arguments for parser
    test_directory = os.path.dirname(os.path.abspath(__file__))
    targets_fasta = os.path.join(test_directory, 'test-data/targets.fasta')
    targets_gff = os.path.join(test_directory, 'test-data/targets.gff')
    targets = os.path.join(test_directory, 'test-data/targets')
    args = parse_args(['-i', targets_fasta, '-g', targets_gff, '-T', targets])

    # run design_primers and convert result to string
    result = '\n'.join(design_primers(args))

    # expected result
    expected = """SNP_Target_ID Position Ref_base Variant_base Amplicon_bp PRIMER_LEFT_SEQUENCE PRIMER_RIGHT_SEQUENCE ref_melt_Tm var_melt_Tm Tm_difference
k69_93535:SAMTOOLS:SNP:1336 1336 G A 149 GAACGCGAGAAGGAGTGAAG GCAACCCAGGTTTCAACTCC 0 0 0
k69_93535:SAMTOOLS:SNP:1336 1336 G A 252 GTGTCGGGTGAAACGATACC GCAACCCAGGTTTCAACTCC 0 0 0
k69_93535:SAMTOOLS:SNP:1336 1336 G A 183 GAACGCGAGAAGGAGTGAAG GAAGGAACACCGCCATTAGG 0 0 0
k69_93535:SAMTOOLS:SNP:1336 1336 G A 184 GAACGCGAGAAGGAGTGAAG GGAAGGAACACCGCCATTAG 0 0 0
k69_93535:SAMTOOLS:SNP:1336 1336 G A 286 GTGTCGGGTGAAACGATACC GAAGGAACACCGCCATTAGG 0 0 0
k69_98089:SAMTOOLS:SNP:625 625 A G 227 GGAGAAGGTCGAGGTCAGC ACGGCCGAATATACATACAACG 0 0 0
k69_98089:SAMTOOLS:SNP:625 625 A G 294 GGGAGACCGATCAGTGTTGG AACGGCCGAATATACATACAACG 0 0 0
k69_98089:SAMTOOLS:SNP:625 625 A G 225 AGAAGGTCGAGGTCAGCG ACGGCCGAATATACATACAACG 0 0 0
k69_98089:SAMTOOLS:SNP:625 625 A G 292 GGGAGACCGATCAGTGTTGG CGGCCGAATATACATACAACGTC 0 0 0
k69_98089:SAMTOOLS:SNP:625 625 A G 228 GGAGAAGGTCGAGGTCAGC AACGGCCGAATATACATACAACG 0 0 0"""

    assert result == expected


# def test_design_primers_umelt():
#     """
#     Test for function design_primers with umelt functionality.
#     This test is SLOW to run.
#     Uses argument parser to ensure input is in expected form.
#     This is effectively an full integration test.
#     """
#     # create arguments for parser
#     test_directory = os.path.dirname(os.path.abspath(__file__))
#     targets_fasta = os.path.join(test_directory, 'test-data/targets.fasta')
#     targets_gff = os.path.join(test_directory, 'test-data/targets.gff')
#     targets = os.path.join(test_directory, 'test-data/targets')
#     args = parse_args(['-i', targets_fasta, '-g', targets_gff, '-T', targets, '-u'])
#
#     # run design_primers and convert result to string
#     result = '\n'.join(design_primers(args))
#
#     # expected result
#     expected = """SNP_Target_ID Position Ref_base Variant_base Amplicon_bp PRIMER_LEFT_SEQUENCE PRIMER_RIGHT_SEQUENCE ref_melt_Tm var_melt_Tm Tm_difference
# k69_93535:SAMTOOLS:SNP:1336 1336 G A 149 GAACGCGAGAAGGAGTGAAG GCAACCCAGGTTTCAACTCC 88.7334273625 88.6833568406 0.0500705218618
# k69_93535:SAMTOOLS:SNP:1336 1336 G A 252 GTGTCGGGTGAAACGATACC GCAACCCAGGTTTCAACTCC 89.2341325811 89.1840620592 0.0500705218618
# k69_93535:SAMTOOLS:SNP:1336 1336 G A 183 GAACGCGAGAAGGAGTGAAG GAAGGAACACCGCCATTAGG 89.2341325811 89.0839210155 0.150211565585
# k69_93535:SAMTOOLS:SNP:1336 1336 G A 184 GAACGCGAGAAGGAGTGAAG GGAAGGAACACCGCCATTAG 89.2341325811 89.1339915374 0.100141043724
# k69_93535:SAMTOOLS:SNP:1336 1336 G A 286 GTGTCGGGTGAAACGATACC GAAGGAACACCGCCATTAGG 89.3843441467 89.284203103 0.100141043724
# k69_98089:SAMTOOLS:SNP:625 625 A G 227 GGAGAAGGTCGAGGTCAGC ACGGCCGAATATACATACAACG 85.7792665726 86.2799717913 0.500705218618
# k69_98089:SAMTOOLS:SNP:625 625 A G 294 GGGAGACCGATCAGTGTTGG AACGGCCGAATATACATACAACG 85.7792665726 86.2799717913 0.500705218618
# k69_98089:SAMTOOLS:SNP:625 625 A G 225 AGAAGGTCGAGGTCAGCG ACGGCCGAATATACATACAACG 85.7792665726 86.2799717913 0.500705218618
# k69_98089:SAMTOOLS:SNP:625 625 A G 292 GGGAGACCGATCAGTGTTGG CGGCCGAATATACATACAACGTC 85.5789844852 86.2799717913 0.700987306065
# k69_98089:SAMTOOLS:SNP:625 625 A G 228 GGAGAAGGTCGAGGTCAGC AACGGCCGAATATACATACAACG 85.7792665726 86.2799717913 0.500705218618"""
#
#     assert result == expected