#!/usr/bin/env python3



import argparse
import pybedtools as pb
import pandas as pd
import re
import ast
import sys
from design import *

def design_primers(args):
    """
    Design primers.
    """
    for arg in vars(args):
        print("{}={}".format(arg, getattr(args, arg)))
    
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
        
    # Get dictionnary from string if needed
    if args.region:
        args.region = ast.literal_eval(args.region)
        
    
    
    
    designer = PrimerDesign(reference_file=args.reference,
                            annotation_file=args.annotation, 
                            description=args.description, 
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
            designer.run_p3(region=args.region)
            pb.cleanup()
        if args.output_directory != "stdout":
            designer.cleanoutput()
    else:
        designer.run_p3(region=args.region)
        if args.output_directory != "stdout":
            designer.cleanoutput()
    
#     # Create a VcfPrimerDesign object
#     designer = d.VcfPrimerDesign(args.fastaref, args.vcffile, args.desc)
    
#     # Open the bed target file
#     targets = pb.BedTool((args.bedfile_path + args.bedfile))
    
    
#     # Set the default primer3 parameters
#     p3_globals = {
#         'PRIMER_OPT_SIZE': 20,
#         'PRIMER_PICK_INTERNAL_OLIGO': 0,
#         'PRIMER_INTERNAL_MAX_SELF_END': 8,
#         'PRIMER_MIN_SIZE': 18,
#         'PRIMER_MAX_SIZE': 25,
#         'PRIMER_OPT_TM': 60.0,
#         'PRIMER_MIN_TM': 55.0,
#         'PRIMER_MAX_TM': 63.0,
#         'PRIMER_MIN_GC': 20.0,
#         'PRIMER_MAX_GC': 80.0,
#         'PRIMER_MAX_POLY_X': 100,
#         'PRIMER_INTERNAL_MAX_POLY_X': 100,
#         'PRIMER_SALT_MONOVALENT': 50.0,
#         'PRIMER_DNA_CONC': 50.0,
#         'PRIMER_MAX_NS_ACCEPTED': 0,
#         'PRIMER_MAX_SELF_ANY': 12,
#         'PRIMER_MAX_SELF_END': 8,
#         'PRIMER_PAIR_MAX_COMPL_ANY': 12,
#         'PRIMER_PAIR_MAX_COMPL_END': 8,
#         'PRIMER_PRODUCT_SIZE_RANGE': [200,300],
#         'PRIMER_NUM_RETURN' : 2
#     }
    
#     # Modify p3_globals if specified
#     if args.p3_globals:
        
#         # Get dictionnary from string
#         custom_p3_globals = ast.literal_eval(args.p3_globals)
        
#         for key in custom_p3_globals:
#             p3_globals[key] = custom_p3_globals[key]
    
#     # If output_path is left as default print everything to standard output
#     if args.output_path == "stdout":
#         print('\t'.join([
#             'CHR',
#             'START',
#             'END',
#             'SEQUENCE_START',
#             'SEQUENCE_END',
#             'TARGET_START',
#             'TARGET_END',
#             'PRIMER_LEFT',
#             'PRIMER_LEFT_SEQUENCE',
#             'PRIMER_RIGHT',
#             'PRIMER_RIGHT_SEQUENCE'
#         ]))
#         for X in targets:
#             for Y in P3.run_P3(global_dict=p3_globals,target_dict = designer.getseqslicedict(X,p3_globals['PRIMER_PRODUCT_SIZE_RANGE'][1])):
#                 reg = re.split(':|-',Y['AMPLICON_REGION'])
#                 seq = re.split(':|-',Y['SEQUENCE_ID'])
#                 targ = re.split(':|-',Y['TARGET_ID'])
#                 print('\t'.join([
#                     reg[0],
#                     reg[1],
#                     reg[2],
#                     seq[1],
#                     seq[2],
#                     targ[1],
#                     targ[2],
#                     str(Y['PRIMER_LEFT']),
#                     Y['PRIMER_LEFT_SEQUENCE'],
#                     str(Y['PRIMER_RIGHT']),
#                     Y['PRIMER_RIGHT_SEQUENCE']
#                 ]))
                
#     else:
#         # Create a pandas DataFrame that will hold the primers data
#         df = pd.DataFrame(columns=[
#             'CHR',
#             'START',
#             'END',
#             'SEQUENCE_START',
#             'SEQUENCE_END',
#             'TARGET_START',
#             'TARGET_END',
#             'PRIMER_LEFT',
#             'PRIMER_LEFT_SEQUENCE',
#             'PRIMER_RIGHT',
#             'PRIMER_RIGHT_SEQUENCE'
#         ])

#         # Populate the data frame with the run_P3 output from the targets
#         for X in targets:
#             for Y in P3.run_P3(global_dict=p3_globals,target_dict = designer.getseqslicedict(X,300)):
#                 reg = re.split(':|-',Y['AMPLICON_REGION'])
#                 seq = re.split(':|-',Y['SEQUENCE_ID'])
#                 targ = re.split(':|-',Y['TARGET_ID'])
#                 df = df.append({
#                     'CHR':reg[0],
#                     'START':reg[1],
#                     'END':reg[2],
#                     'SEQUENCE_START':seq[1],
#                     'SEQUENCE_END':seq[2],
#                     'TARGET_START':targ[1],
#                     'TARGET_END':targ[2],
#                     'PRIMER_LEFT':Y['PRIMER_LEFT'],
#                     'PRIMER_LEFT_SEQUENCE':Y['PRIMER_LEFT_SEQUENCE'],
#                     'PRIMER_RIGHT':Y['PRIMER_RIGHT'],
#                     'PRIMER_RIGHT_SEQUENCE':Y['PRIMER_RIGHT_SEQUENCE']
#                 }, ignore_index=True)
        
#         df.to_csv("{X.output_path}primers_{X.bedfile}.csv".format(X=args))
            
            

def main():
    """
    Main function that will be automaticaly executed.
    It parses the arguments and runs the design_primers function
    """
    parser = argparse.ArgumentParser(description="""Design primers""")
    parser.add_argument("-reference", help="fasta reference file", dest="reference", type=str, required=True)
    parser.add_argument("-annotation", help="vcf annotation file", dest="annotation", type=str, required=True)
    parser.add_argument("-description", help="description", dest="description", type=str, required=True)
    parser.add_argument("-targets", help="target bed file", dest="targets", type=str, default=None)
    parser.add_argument("-output_directory", help="output path", dest="output_directory", type=str, default="stdout")
    parser.add_argument("-p3_globals", help="dictionnary to specify primer3 parameters",
                        dest='p3_globals', type=str, default=None)
    parser.add_argument("-amplicon_size_range", help="amplicon size range", dest="amplicon_size_range",
                        nargs=2, type=int, default=None)
    parser.add_argument("-primer_size_range", help="primer size range", dest="primer_size_range",
                        nargs=2, type=int, default=None)
    parser.add_argument("-region", help="working region as '{\"CHR1\":[[start, end], ...], ...}'",
                        dest="region", type=str, default=None)
    parser.add_argument("-interactive", help="stream targets from stdin", dest='interactive',
                        action='store_true', required=False)
    args = parser.parse_args()
    design_primers(args)

if __name__ == '__main__':
    main()