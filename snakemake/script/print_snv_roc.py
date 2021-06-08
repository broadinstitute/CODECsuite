#!/usr/bin/env python

import argparse
import logging
import sys
import glob
import os
import pandas as pd

logger = logging.getLogger("{}".format(__file__))
def get_arguments():

    parser = argparse.ArgumentParser(prog="print whether mutant families in CpG context or not", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcfeval_folder", type=str, help="vcfeval output folder")
    args = parser.parse_args()
    return args

def process(opts):
    outs = glob.glob(os.path.join(opts.vcfeval_folder, "*/*/snp_roc.tsv.gz"))
    total_df = pd.DataFrame()
    for o in outs:
        fields = o.split('/')
        df = pd.read_csv(o, sep='\t', skiprows=6)
        df['cutoff'] = fields[-2]
        df['vaf'] =  fields[-3][2:]
        if total_df.empty:
            total_df = df
        else:
            total_df = pd.concat([total_df, df])
        total_df.sort_values(by=['vaf', 'cutoff'], inplace=True)
        total_df = total_df[['vaf','cutoff', 'precision', 'sensitivity', 'false_positives', 'false_negatives', 'true_positives_baseline', 'f_measure']]
    print (total_df)
    total_df.to_csv("summary.tsv", sep='\t', index=False)


if __name__ == '__main__':
    sys.exit(process(get_arguments()))

