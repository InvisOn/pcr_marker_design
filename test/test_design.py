import pytest
from pcr_marker_design import design as d
from pyfaidx import Fasta


class TestDesign:
    def test_getseqslice(self):
        test_seq="./test/test-data/targets.fasta"
        annfile= "./test/test-data/targets.gff"
        k69_93535_250_290='GACAAAGAGAAAATCCTCAAATCCGGCCTCGTCAACCACA'
        designer = d.PrimerDesign(test_seq,annfile)
        seqslice=designer.getseqslice('k69_93535',250,290)
        assert seqslice == k69_93535_250_290
