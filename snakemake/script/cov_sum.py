#!/usr/bin/env python

import argparse
import logging
import os
import sys

logger = logging.getLogger("{}".format(__file__))

def process():
    cov_count = {}
    for line in sys.stdin:
        if line.startswith('REF'):
            continue
        cols = line.strip().split("\t")
        cov = int(cols[2])
        if cov not in cov_count:
            cov_count[cov] = 1
        else:
            cov_count[cov] += 1

    for k in sorted(cov_count.keys()):
        print(k, cov_count[k], sep='\t')

if __name__ == '__main__':
    sys.exit(process())

