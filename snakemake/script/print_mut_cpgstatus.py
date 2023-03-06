#!/usr/bin/env python

import argparse
import logging
import os
import sys
import pysam
import pandas as pd

logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="print whether mutant families in CpG context or not", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mutfam", type=str, help="mutant families")
    parser.add_argument("ref", type=str, help="reference file")
    args = parser.parse_args()
    return args

def absolute_path(path):
    """convert relative path to absolute path"""
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(os.getcwd(), path)

def process(opt):
    mf = pd.read_csv(opt.mutfam, sep='\t', low_memory=False)
    mf = mf[mf['pass_filter'] == 1].reset_index()
    ncol = mf.shape[1]
    mf['CpG'] = 'NA'
    fasta = pysam.FastaFile(opt.ref)
    for idx, row in mf.iterrows():
        if row['ref_allele'] == 'C':
            nextb = fasta.fetch(str(row['contig']), row['position'], row['position'] + 1)
            if nextb == 'G' or nextb == 'g':
                mf.iloc[idx, ncol] = '1'
            else:
                mf.iloc[idx, ncol] = '0'
        if row['ref_allele'] == 'G':
            nextb = fasta.fetch(str(row['contig']), row['position']-2, row['position'] - 1)
            if nextb == 'C' or nextb == 'c':
                mf.iloc[idx, ncol] = '1'
            else:
                mf.iloc[idx, ncol] = '0'
    print(mf)

if __name__ == '__main__':
    sys.exit(process(get_arguments()))

