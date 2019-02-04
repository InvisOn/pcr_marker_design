#! /usr/bin/python

# Just a small script to use the PrimerDesign class from the shell
# and set all the parameters for recurring use.


# IMPORTS

from design import *


# PARAMETERS

p3_globals = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 55.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 20.0,
        'PRIMER_MAX_GC': 80.0,
        'PRIMER_MAX_POLY_X': 100,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_NUM_RETURN' : 3
}

primer_size_range = (18, 25)
amplicon_size_range = (200, 300)

reference_file = "/output/genomic/plant/Actinidia/chinensis/CK51F3_01/Genome/Assembly/PS1/1.68.5/AllChromosomes/PS1.1.68.5.fasta"
description = "PS1.1.68.5"
annotation_file = "/output/genomic/plant/Actinidia/chinensis/Resequencing/Variants/PS1.1.68.5/52DiploidGenomes/Combined_diploidCK_basic_NS30_Q50_SAFR3_DP50_PAIR0.8_PS1.1.68.5_ann.vcf.gz"
targets_file = "/workspace/cfljam/PrimerDesignJobs/PVR_Marker_design/BedFilesCK/CHR1.targets.bed"
output_dir = "tmp"

region = {"CHR1":[[0, 5000]]}





# EXECUTION

# For more detail on the different possiblities go to Example/Example.ipynb

d = PrimerDesign(reference_file,
                 description, 
                 annotation_file,
                 targets_file, 
                 output_dir,
                 amplicon_size_range,
                 primer_size_range,
                 p3_globals)

d.run_p3(region)
d.cleanoutput()
