#!/usr/bin/env python
import argparse
import logging
import sys
import pandas as pd

logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="convert mutect file to maf file", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mutect", type=str, help="mutect output")
    parser.add_argument("-t", "--tlod", type=int, default=10, help="t lod cutoff. 10 default")
    args = parser.parse_args()
    return args


def process(opts):

    header= ['Hugo_Symbol', 'Chromosome', 'Start_position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'probed', 'tumor_f',
             't_ref_count', 't_alt_count', 'n_ref_count', 'n_alt_count']
    print('\t'.join(header))
    df = pd.read_csv(opts.mutect, sep='\t', low_memory=False, skiprows=1)
    for idx, row in df.iterrows():
        if row['judgement'] == "KEEP" and row['t_lod_fstar'] >= opts.tlod:
            print('NA',
                  row['contig'],
                  row['position'],
                  row['ref_allele'],
                  row['alt_allele'],
                  1,
                  row['tumor_f'],
                  row['t_ref_count'],
                  row['t_alt_count'],
                  row['n_ref_count'],
                  row['n_alt_count'],
                  sep='\t')


if __name__ == '__main__':
    sys.exit(process(get_arguments()))