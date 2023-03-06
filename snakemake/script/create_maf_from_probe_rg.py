#!/usr/bin/env python
import argparse
import logging
import sys
import pysam

logger = logging.getLogger("{}".format(__file__))
def read_bed(bedfile):
    """ Creates generator from bed file or interval_list """
    logger.info("Reading region file...")
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

    parser = argparse.ArgumentParser(prog="Parse Fastp json result(s)", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--bed", type=str, help="json file output by fastp", required=False)
    parser.add_argument("--ref", type=str, help="json file output by fastp", required=False)
    args = parser.parse_args()
    return args

def process(opts):
    header= ['Hugo_Symbol', 'Chromosome', 'Start_position', 'Reference_Allele', 'Tumor_Seq_Allele2', 'probed']
    reffa = pysam.FastaFile(opts.ref)
    print('\t'.join(header))
    for chrom, s, e in read_bed(opts.bed):
        halfw = int(e) - int(s)
        start = int(s) + halfw
        REF = reffa.fetch(chrom, start , start + 1)
        ALT = REF
        line = ['NA', chrom, str(start + 1), REF, ALT, '1']
        print('\t'.join(line))

if __name__ == '__main__':
    sys.exit(process(get_arguments()))
