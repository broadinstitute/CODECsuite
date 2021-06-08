#!/usr/bin/env python

import argparse
import logging
import os
import sys
import pandas as pd

logger = logging.getLogger("{}".format(__file__))

def rc_base(b):
    rc_hash = {'T': 'A', 'A': 'T', 'C': 'G', 'G': 'C'}
    return rc_hash[b]

def get_arguments():

    parser = argparse.ArgumentParser(prog="foo", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mutfam", type=str, help="mutfam files")
    parser.add_argument("met", type=str, help="met files")
    args = parser.parse_args()
    return args

def absolute_path(path):
    """convert relative path to absolute path"""
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(os.getcwd(), path)

def process(opts) :
    spname = opts.mutfam.replace( "_mutant_families.txt", "")
    mutfam = pd.read_csv(opts.mutfam, sep='\t', low_memory = False)
    mutfam = mutfam[mutfam['snv_base_qual'] > 0]
    if not mutfam.empty:
        mutfam['refbase'] = mutfam.apply(lambda x : rc_base(x['ref']) if x['ref'] == 'T' or x['ref'] == 'G' else x['ref'], 1)
        mutfam['mutbase'] = mutfam.apply(lambda x : rc_base(x['alt']) if x['ref'] == 'T' or x['ref'] == 'G' else x['alt'], 1)
        mutfam = mutfam.groupby(['refbase', 'mutbase']).read_name.count()
    met = pd.read_csv(opts.met, sep='\t', low_memory = False, nrows=1)
    den = {'A' : 0, 'C': 0}
    den['A'] = met['n_A_eval'] + met['n_T_eval']
    den['C'] = met['n_C_eval'] + met['n_G_eval']


    print("sample_id", "ref_allele", "alt_allele", "total_base_errors", "total_bases_eval", sep='\t')
    for ref,alt in zip(['A'] * 3 + ['C'] * 3, ['C', 'G', 'T', 'A', 'G', 'T']):
        if (ref, alt) in mutfam.index:
            print(spname, ref, alt, mutfam[ref, alt], int(den[ref]), sep='\t')
        else:
            print(spname, ref, alt, 0, int(den[ref]), sep='\t')

if __name__ == '__main__':
    sys.exit(process(get_arguments()))

