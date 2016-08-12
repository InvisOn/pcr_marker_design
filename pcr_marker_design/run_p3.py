#!/usr/bin/python
##run primer3 by passing Python dictionary

#run_P3.py

p3_globals={
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_MAX_SELF_END': 8,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
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
        'PRIMER_PRODUCT_SIZE_RANGE': [[60,200]],
    }

import primer3

##call P3 with dict of args, returns dict, no exception handling
def run_P3(target_dict,global_dict):
    P3_dict=primer3.bindings.designPrimers(target_dict,global_dict)
    ##return iterable list
    primer_list=[dict(PRIMER_RIGHT_SEQUENCE=P3_dict.get('PRIMER_RIGHT_'+ str(X) + '_SEQUENCE'),\
    PRIMER_LEFT=P3_dict.get('PRIMER_LEFT_'+ str(X) ),\
    PRIMER_RIGHT=P3_dict.get('PRIMER_RIGHT_'+ str(X) ),\
    PRIMER_LEFT_SEQUENCE=P3_dict.get('PRIMER_LEFT_'+ str(X) + '_SEQUENCE')) for X in range(0,int(P3_dict.get('PRIMER_RIGHT_NUM_RETURNED'))-1)]
    return(primer_list)
