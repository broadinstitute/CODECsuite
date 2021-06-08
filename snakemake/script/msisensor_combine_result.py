#!/usr/bin/env python
import argparse
import logging
import sys
import numpy as np
import os
import pprint
from collections import defaultdict

logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="combine msisensor results", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("table", type=str, help="msisensor all table")
    parser.add_argument("dist", type=str, help="msisensor dist file")
    args = parser.parse_args()
    return args


def process(opts):
    collector = defaultdict(lambda : defaultdict( lambda : defaultdict(int)))
    pp = pprint.PrettyPrinter(indent=4)
    hp_agg = defaultdict(lambda : defaultdict(int))
    with open(opts.dist, 'r') as fh:
        for nrow, line in enumerate(fh):
            if nrow % 2 == 0:
                cols = line.split()
                chrom = cols[0]
                pos = int(cols[1])
            else:
                duo = line.split(": ")
                vec = [int(x) for x in duo[1].split()]
                nonzeroind = np.nonzero(vec)[0]
                for ind in nonzeroind:
                    collector[chrom][pos][ind + 1] = vec[ind]
    with open(opts.table, 'r') as fh:
        for line in fh:
            if line.startswith("chromosome"):
                continue
            cols = line.strip().split('\t')
            # no germline variatns
            if cols[-1] == '0':
                repeat_times = cols[3]
                chrom = cols[0]
                pos = int(cols[1])
                repeat = cols[4]
                unit_cnt = collector[chrom][pos]
                for k,v in unit_cnt.items():
                    key = repeat + "," + repeat_times
                    hp_agg[key][k] += v
                    if repeat == 'A' and repeat_times == '10':
                        if k == 5:
                            print(line)

    #pp.pprint(hp_agg)



if __name__ == '__main__':
    sys.exit(process(get_arguments()))
