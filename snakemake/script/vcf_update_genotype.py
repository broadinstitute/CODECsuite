#!/usr/bin/env python

import argparse
import logging
import sys
from cyvcf2 import VCF, Writer

logger = logging.getLogger("{}".format(__file__))
def get_arguments():

    parser = argparse.ArgumentParser(prog="change the variant genotype", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("vcf", type=str, help="vcf input")
    parser.add_argument("out", type=str, help="vcf output")
    args = parser.parse_args()
    return args

def process(opts):
    inputvcf = VCF(opts.vcf)
    writer = Writer(opts.out, inputvcf, "wz")
    for record in inputvcf:
        assert(len(record.gt_types) == 1)
        if record.gt_types[0] == 0:
            record.genotypes = [[0,1,False]]
            writer.write_record(record)
        else:
            writer.write_record(record)
    writer.close()

if __name__ == '__main__':
    sys.exit(process(get_arguments()))

