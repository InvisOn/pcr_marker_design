#!/usr/bin/env python3

# Gabriel Besombes January 2019

# This file includes a parsing function to use the PrimerDesign class from design.py
# in a bash shell.

import argparse
import pybedtools as pb
import pandas as pd
import re
import ast
import sys
from design import *

def design_primers(args):
    """
    Design primers using the PrimerDesign class from design.py
    passing it the argumets parsed in main().
    """
    
    if args.region or args.region_file:
        region = {}
    else:
        region = None
    if args.region:
        for l in [X.replace(":", "-").split("-") for X in args.region.split(",")]:
            if l[0] in region:
                region[l[0]].append([int(l[1]), int(l[2])])
            else:
                region[l[0]] = [[int(l[1]), int(l[2])]]
    if args.region_file:
        bed_reg = pb.BedTool(args.region_file)
        for reg in bed_reg:
            if reg.chrom in region:
                region[reg.chrom].append([int(reg.start), int(reg.stop)])
            else:
                region[reg.chrom] = [[int(reg.start), int(reg.stop)]]
    
    
#     for arg in vars(args):
#         print("{}={}".format(arg, getattr(args, arg)))
    
    # Set the default primer3 parameters
    p3_globals = {
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 0,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
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
        'PRIMER_PRODUCT_SIZE_RANGE': [200,300],
        'PRIMER_NUM_RETURN' : 2
    }
    
    # Modify p3_globals if specified
    if args.p3_globals:
        
        # Get dictionnary from string and update p3_globals
        p3_globals.update(ast.literal_eval(args.p3_globals))
        
    
    # Create primer design
    designer = PrimerDesign(reference_file=args.reference,
                            description=args.description, 
                            annotation_file=args.annotation, 
                            targets_file=args.targets, 
                            output_dir=args.output_directory,
                            amplicon_size_range=args.amplicon_size_range,
                            primer_size_range=args.primer_size_range,
                            p3_globals=args.p3_globals)
    
    # If the -i flag was specifyed use the stream from stdin as targets
    if args.interactive:
        for line in sys.stdin:
            target = pb.BedTool(line, from_string=True)
            designer.targets = target.fn
            designer.run_p3(region=region)
            pb.cleanup()
        if args.output_directory != "stdout":
            designer.cleanoutput()
    else:
        if not args.targets:
            designer.gettargets(region=region, write="targets.bed")
            designer.targets = "targets.bed"
        designer.run_p3(region=region)
        if args.output_directory != "stdout":
            designer.cleanoutput()     
            

def main():
    """
    Main function that will be automaticaly executed.
    It parses the arguments and runs the design_primers function
    """
    parser = argparse.ArgumentParser(description="""Design primers""")
    parser.add_argument("-f", "--fasta-ref",
                        help="FASTA FILE", dest="reference", type=str, required=True)
    parser.add_argument("-a", "--annotations",
                        help="VCF ANNOTATION FILE", dest="annotation", type=str, default=None)
    parser.add_argument("-d", "--description",
                        help="description", dest="description", type=str, required=True)
    parser.add_argument("-T", "--targets-file",
                        help="TARGETS BED FILE", dest="targets", type=str, default=None)
    parser.add_argument("-o", "--output_directory",
                        help="output path", dest="output_directory", type=str, default="stdout")
    parser.add_argument("-p3_globals", help="dictionnary to specify primer3 parameters",
                        dest='p3_globals', type=str, default=None)
    parser.add_argument("-amplicon_size_range", help="amplicon size range", dest="amplicon_size_range",
                        nargs=2, type=int, default=None)
    parser.add_argument("-primer_size_range", help="primer size range", dest="primer_size_range",
                        nargs=2, type=int, default=None)
    parser.add_argument("-r", "--regions",
                        help="chr:from-to,chr:-to,chr:from-", dest="region", type=str, default=None)
    parser.add_argument("-R", "--regions-file",
                        help="BED FILE OF REGIONS", dest="region_file", type=str, default=None)
    parser.add_argument("-interactive", help="stream targets from stdin", dest='interactive',
                        action='store_true', required=False)
    args = parser.parse_args()
    design_primers(args)

# If this file is not imported but executed directly, execute main().
if __name__ == '__main__':
    main()