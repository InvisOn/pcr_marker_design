# Gabriel Besombes January 2019

# Tests for the different elements of design_bash.py.
# These can all be executed using the "python -m pytest" command in a bash shell

import pytest
import pybedtools as pb
import subprocess as sp
import pandas as pd
import pyfaidx as pf
import vcf
from ast import literal_eval
from Bio.Seq import Seq


reference = "/output/genomic/plant/Actinidia/chinensis/CK51F3_01/Genome/Assembly/PS1/1.68.5/AllChromosomes/PS1.1.68.5.fasta"
annotations = "/output/genomic/plant/Actinidia/chinensis/Resequencing/Variants/PS1.1.68.5/52DiploidGenomes/Combined_diploidCK_basic_NS30_Q50_SAFR3_DP50_PAIR0.8_PS1.1.68.5_ann.vcf.gz"
description = "PS1.1.68.5"
targets = "./test/test-data/targets.bed"


def test_design_bash():
    """
    pass
    """
    targets_file = pb.BedTool(targets)
    out = open("test/test_general_design_bash.out.csv", mode="+w")
    err = open("test/test_general_design_bash.err", mode="+w")
    
    # Test program with different sizes of targets/amplicons
    l = []
    for target in targets_file:
        l.append(sp.call(['./design_bash.py',
                          '-f', reference,
                          '-a', annotations,
                          '-d', description,
                          '-amplicon_size_range',
                          str(int(target.stop-target.start-100)),
                          str(int(target.stop-target.start+600)),
                          '-T', pb.BedTool("{} {} {}".format(target.chrom,
                                                             target.start,
                                                             target.stop),
                                           from_string=True).fn],
                         stdout=out,
                         stderr=err))

    out.close()
    err.close()
    # Assert that all finished normaly
    # if not, you can go see in test/test_general_design_bash.err
    # for detailed error
    assert(not any(l))
    
    # Assert that the primers don't overlap with features on the vcf
    out_csv = pd.read_csv("test/test_general_design_bash.out.csv")
    l = []
    for i,row in out_csv.iterrows():
        if row[0] != "REF_OFFSET":
            for j in range(0, int(row["PRIMER_LEFT_NUM_RETURNED"])):
                l.append([row["CHROMOSOME"],
                         literal_eval(row["PRIMER_LEFT_" + str(j)])[0] + int(row["REF_OFFSET"]),
                         literal_eval(row["PRIMER_LEFT_" + str(j)])[0] + int(row["REF_OFFSET"]) + literal_eval(row["PRIMER_LEFT_" + str(j)])[1],
                         row["PRIMER_LEFT_" + str(j) + "_SEQUENCE"]])
            for j in range(0, int(row["PRIMER_RIGHT_NUM_RETURNED"])):
                l.append([row["CHROMOSOME"],
                         literal_eval(row["PRIMER_RIGHT_" + str(j)])[0] + int(row["REF_OFFSET"]) - literal_eval(row["PRIMER_RIGHT_" + str(j)])[1] + 1,
                         literal_eval(row["PRIMER_RIGHT_" + str(j)])[0] + int(row["REF_OFFSET"]) + 1,
                         row["PRIMER_RIGHT_" + str(j) + "_SEQUENCE"]])
    anno = vcf.Reader(filename=annotations)
    for primer in l:
        annots = [X for X in anno.fetch(primer[0], primer[1], primer[2])]
        assert(len(annots) == 0)
    
    # Assert that the sequences and coordinates match
    ref = pf.Fasta(reference)
    for primer in l:
        seq = Seq(ref[primer[0]][primer[1]:primer[2]].seq)
        assert(seq == primer[3] or seq.reverse_complement() == primer[3])