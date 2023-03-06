import argparse
import sys
from math import ceil

def read_bed(bedfile):
    """ Creates generator from bed file or interval_list """
    interval_list = bedfile.endswith("interval_list")
    with open(bedfile, "r") as bed:
        for line in bed:
            if line.startswith("@"):
                continue
            line = line.strip()
            chrom, start, stop = line.split()[0:3]
            start, stop = int(start), int(stop)
            if interval_list:
                start -= 1
            yield chrom, start, stop

def get_arguments():

    parser = argparse.ArgumentParser(prog="get midpoint from interval and print like a bed file", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bed", type=str, help="interval list or bed file")
    parser.add_argument("out", type=str, help="output for bed file")
    args = parser.parse_args()
    return args

def process(options):
    outh = open(options.out, 'w')
    for chrom, start, stop in read_bed(options.bed):
        mid = ceil((start+stop)/2)
        print(chrom, mid-1, mid, sep='\t', file=outh)

if __name__ == '__main__':
    sys.exit(process(get_arguments()))
