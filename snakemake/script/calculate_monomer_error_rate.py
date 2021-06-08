#!/usr/bin/env python
import pandas as pd
import argparse
import numpy as np


def parse_cl_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tnd", help="trinucleotide count result", required=True)
    parser.add_argument("--mutfam", help="mutant family result", required=True)
    return parser.parse_args()


def rc_base(b):
    rc_hash = {'T': 'A', 'A': 'T', 'C': 'G', 'G': 'C'}
    return rc_hash[b]


def process(opts) :
    spname = opts.mutfam.replace( "_mutant_families.txt", "")
    assert (spname == opts.tnd.replace( "_trinucleotide_context.txt", ""))
    tnd = pd.read_csv(opts.tnd, sep='\t', low_memory=False)
    tnd['middle_base'] = tnd.context.apply(lambda x: x[1])
    tnd = tnd.groupby('middle_base')['base_count'].agg(np.sum)
    mutfam = pd.read_csv(opts.mutfam, sep='\t', low_memory = False)
    mutfam = mutfam[mutfam['pass_filter'] == 1]
    if not mutfam.empty:
        mutfam['refbase'] = mutfam.apply(lambda x : rc_base(x['ref_allele']) if x['ref_allele'] == 'T' or x['ref_allele'] == 'G' else x['ref_allele'], 1)
        mutfam['mutbase'] = mutfam.apply(lambda x : rc_base(x['alt_allele']) if x['ref_allele'] == 'T' or x['ref_allele'] == 'G' else x['alt_allele'], 1)
        mutfam = mutfam.groupby(['refbase', 'mutbase']).pass_filter.count()

    print("sample_id", "ref_allele", "alt_allele", "total_base_errors", "total_bases_eval", sep=',')
    for ref,alt in zip(['A'] * 3 + ['C'] * 3, ['C', 'G', 'T', 'A', 'G', 'T']):
        if (ref, alt) in mutfam.index:
            print(spname, ref, alt, mutfam[ref, alt], tnd.loc[ref], sep=',')
        else:
            print(spname, ref, alt, 0, tnd.loc[ref], sep=',')


if __name__ == "__main__":
    process(parse_cl_args())
