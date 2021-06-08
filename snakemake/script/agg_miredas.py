#!/usr/bin/env python
import argparse
import logging
import sys
import pandas as pd

logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="aggregate miredas reult(s)", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-o", "--out_prefix", type=str, help="output prefix", required=True)
    parser.add_argument("-e", "--error_metrics", nargs = "+", type=str, help="error_metrics files", required=True)
    parser.add_argument("-t", "--trinuc_context", nargs = "+", type=str, help="tri-nucleotide context files", required=True)
    parser.add_argument("-m", "--mutant_families", nargs = "+", type=str, help="mutant_families files", required=True)
    #parser.add_argument("--si_intermol_bamout", type=str, help="output bam for singleinsert intermolcular reads", required=True)
    args = parser.parse_args()
    return args

def concat_df(inputs):
    totaldf = pd.DataFrame()
    for ii in inputs:
        df = pd.read_csv(ii, sep='\t', low_memory=False)
        if totaldf.empty:
            totaldf = df
        else:
            totaldf = pd.concat([totaldf, df])
    return totaldf

def process(opts):
    ret = concat_df(opts.error_metrics)
    ret = ret.groupby('sample_id').agg({'total_base_errors': sum, 'total_bases_eval': sum}).reset_index()
    ret.to_csv(f"{opts.out_prefix}_error_metrics.txt", sep='\t', index=False)

    ret = concat_df(opts.trinuc_context)
    ret = ret.groupby('context').agg({'site_count': sum, 'base_count': sum}).reset_index()
    ret.to_csv(f"{opts.out_prefix}_trinucleotide_context.txt", sep='\t', index=False)

    ret = concat_df(opts.mutant_families)
    ret.to_csv(f"{opts.out_prefix}_mutant_families.txt", sep='\t', index=False)

if __name__ == '__main__':
    sys.exit(process(get_arguments()))