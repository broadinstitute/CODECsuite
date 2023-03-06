#!/usr/bin/env python

import argparse
import logging
import os
import sys
import pysam

def get_arguments():

    parser = argparse.ArgumentParser(prog="", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bam", type=str, help="input bam file")
    parser.add_argument("out", type=str, help="output bam file")
    args = parser.parse_args()
    return args

def process(options):
    inbam = pysam.AlignmentFile(options.bam, "rb")
    outbam = pysam.AlignmentFile(options.out, "wb", template=inbam)
    for read in inbam:
        if read.is_reverse:
            read.query_qualities = read.query_qualities[::-1]
        outbam.write(read)



if __name__ == '__main__':
    sys.exit(process(get_arguments()))

