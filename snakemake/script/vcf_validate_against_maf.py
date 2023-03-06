#!/usr/bin/env python

import argparse
import logging
import sys
from cyvcf2 import VCF, Writer
import pandas as pd

logger = logging.getLogger("{}".format(__file__))
def get_arguments():

    parser = argparse.ArgumentParser(prog="change the variant genotype", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", type=str, help="vcf input")
    parser.add_argument("out", type=str, help="vcf output")
    parser.add_argument("maf", type=str, help="maf file for truth")
    parser.add_argument("--whitelist", type=str, action='append', help="whitelist varainst in chr:pos format")
    args = parser.parse_args()
    return args

def process(opts):
    inputvcf = VCF(opts.vcf)
    writer = Writer(opts.out, inputvcf, "wz")
    maf = pd.read_csv(opts.maf, sep='\t', low_memory=False)
    maf = maf.astype({'Chromosome' : 'str'})
    for record in inputvcf:
        assert(len(record.gt_types) == 1)
        if not maf[(maf['Chromosome'] == record.CHROM) & (maf['Start_position'] == record.start + 1) & (maf['Tumor_Seq_Allele2'].isin(record.ALT))].empty:
            writer.write_record(record)
        else:
            found = False
            if opts.whitelist:
                for var in opts.whitelist:
                    chr, pos = var.split(':')
                    if record.CHROM == chr and record.start + 1 == int(pos):
                       found = True
                       break
                if found:
                    writer.write_record(record)
                else:
                    print(record.CHROM, record.start, record.start + 1, sep='\t')
            else:
                print(record.CHROM, record.start, record.start + 1, sep='\t')


if __name__ == '__main__':
    sys.exit(process(get_arguments()))

