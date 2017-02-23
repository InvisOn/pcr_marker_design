#!/usr/bin/python

import primer3

# run primer3 by passing Python dictionary

# run_P3.py

p3_globals = {
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
        'PRIMER_NUM_RETURN': 5,
        'PRIMER_INTERNAL_MAX_POLY_X': 100,
        'PRIMER_SALT_MONOVALENT': 50.0,
        'PRIMER_DNA_CONC': 50.0,
        'PRIMER_MAX_NS_ACCEPTED': 0,
        'PRIMER_MAX_SELF_ANY': 12,
        'PRIMER_MAX_SELF_END': 8,
        'PRIMER_PAIR_MAX_COMPL_ANY': 12,
        'PRIMER_PAIR_MAX_COMPL_END': 8,
        'PRIMER_PRODUCT_SIZE_RANGE': [[60, 200]],
    }


# call P3 with dict of args, returns dict, no exception handling


def run_P3(target_dict, global_dict):
    P3_dict = primer3.bindings.designPrimers(target_dict, global_dict)
    # return iterable list
    my_offset=target_dict.get('REF_OFFSET',0)
    my_seq_id=target_dict.get('SEQUENCE_ID')
    primer_list=[]
    for i in range(0, int(P3_dict.get('PRIMER_RIGHT_NUM_RETURNED')) - 1):
        primer_dict=dict(TARGET_ID=target_dict.get('TARGET_ID'),
                         SEQUENCE_ID=my_seq_id)
        primer_dict['PRIMER_LEFT_SEQUENCE']=P3_dict.get('PRIMER_LEFT_' + str(i) + '_SEQUENCE')
        primer_dict['PRIMER_RIGHT_SEQUENCE'] = P3_dict.get('PRIMER_RIGHT_' + str(i) + '_SEQUENCE')
        pr_left=P3_dict.get('PRIMER_LEFT_' + str(i))
        primer_dict['PRIMER_LEFT']=(pr_left[0] + my_offset,pr_left[1])
        pr_right = P3_dict.get('PRIMER_RIGHT_' + str(i))
        primer_dict['PRIMER_RIGHT'] = (pr_right[0] + my_offset, pr_right[1])
        primer_dict['AMPLICON_REGION']= my_seq_id.split(':')[0] + ":" +\
                                        str(pr_left[0] + my_offset + 1) + "-" + \
                                        str(pr_right[0] + my_offset +1)

        primer_list.append(primer_dict)

    return primer_list

