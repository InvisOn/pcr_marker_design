# Gabriel Besombes January 2019

# Tests for the different elements of design.py.
# These can all be executed using the "python -m pytest" command in a bash shell

import pytest
import pyfaidx as pf
from Bio.Seq import Seq
from design import *

reference = "/output/genomic/plant/Actinidia/chinensis/CK51F3_01/Genome/Assembly/PS1/1.68.5/AllChromosomes/PS1.1.68.5.fasta"
annotations = "/output/genomic/plant/Actinidia/chinensis/Resequencing/Variants/PS1.1.68.5/52DiploidGenomes/Combined_diploidCK_basic_NS30_Q50_SAFR3_DP50_PAIR0.8_PS1.1.68.5_ann.vcf.gz"
description = "PS1.1.68.5"
targets = "test_targets.bed"
# targets = "./test/test-data/targets.bed"
    

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
    
    assert(designer.amplicon_size_range == (100, 500))
    assert(designer.primer_size_range == (18, 25))
    # Needs other asserts
    
    
def test_gettargetdict():
    """
    pass
    """
    i = 0
    targets_file = pb.BedTool(targets)
    for target in targets_file:
        designer = PrimerDesign(reference_file=reference,
                                annotation_file=annotations,
                                description=description,
                                targets_file=pb.BedTool("{} {} {}".format(target.chrom,
                                                                          target.start,
                                                                          target.stop),
                                                        from_string=True).fn,
                                output_dir=None,
                                amplicon_size_range=(int(target.stop-target.start-100),
                                                     int(target.stop-target.start+600)),
                                primer_size_range=(18, 25),
                                p3_globals=None)
        d = designer.gettargetdict([target.chrom, target.start, target.stop])
        start = max(0,
                    target.stop + designer.primer_size_range[0] - designer.amplicon_size_range[1])
        stop = min(len(designer.reference[target.chrom]),
                   target.start - designer.primer_size_range[0] + designer.amplicon_size_range[1])
        assert(start == d["REF_OFFSET"])
        assert(designer.reference[target.chrom][start:stop].seq == d["SEQUENCE_TEMPLATE"])
        for excl in d["SEQUENCE_EXCLUDED_REGION"]:
            s = start + excl[0]
            e = start + excl[0] + excl[1]
            annots = [X for X in designer.annotations.fetch(target.chrom, s, e)]
            assert(len(annots) == 1 and
                   excl[1] <= annots[0].stop - annots[0].start)
        if i == 100:
            break
        i += 1
            

def test_gettargetdict():
    """
    pass
    """
    pass

def test_run_p3():
    """
    pass
    """
    pass

def test_cleanoutput():
    """
    pass
    """
    pass

def test_auto():
    """
    pass
    """
    pass

def test_gettargets():
    """
    pass
    """
    pass
