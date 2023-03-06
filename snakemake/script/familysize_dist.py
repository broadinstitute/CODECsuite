#!/usr/bin/env python

import argparse
import logging
import os
import sys
import pysam
import pandas as pd
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
import seaborn as sns
from scipy import stats

logger = logging.getLogger("{}".format(__file__))

def read_pair_generator(bam, region_string=None):
    """
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.
    """
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(until_eof=True):
        if read.is_secondary or read.is_supplementary:
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]

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

def is_overlap(read1, read2):
    if read1.is_unmapped or read2.is_unmapped:
        return False
    if read1.reference_name != read2.reference_name:
        return False
    if read1.reference_start < read2.reference_end and read2.reference_start < read1.reference_end:
        return True

def overlap_len(read1, read2):
    if not is_overlap(read1,read2):
        return 0
    else:
        return min(read1.reference_end, read2.reference_end) - max(read1.reference_start, read2.reference_start)

def overlap_span_ratio(read1, read2):
    ol = overlap_len(read1, read2)
    if ol == 0:
        return 0;
    else:
        span = max(read1.reference_end, read2.reference_end) - min(read1.reference_start, read2.reference_start)
        return ol/span

def get_arguments():

    parser = argparse.ArgumentParser(prog="", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("bam", type=str, help="input bam file")
    parser.add_argument("bed", type=str, help="target bed")
    parser.add_argument("--im_dist_cutoff", default=5000, help="intermolecu reads min distance to be considered", type=int, required=False)
    args = parser.parse_args()
    return args

def absolute_path(path):
    """convert relative path to absolute path"""
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(os.getcwd(), path)

def bin_counts(vec, cut):
    cat = pd.cut(np.array(vec), cut)
    print(cat.value_counts())

def cdf_plot(vec, bins, figure, cumulative, xlab):
    #values, base = np.histogram(vec, bins=bins)
    #cumulative = np.cumsum(values)
    #plt.plot(base[:-1], cumulative, c='blue')
    fig, ax = plt.subplots()
    ax.hist(vec, bins, density=True, cumulative=cumulative, histtype='step')
    ax.set_axisbelow(True)
    ax.yaxis.grid(color='gray', linestyle='dashed')
    ax.xaxis.grid(color='gray', linestyle='dashed')
    #loc = plticker.MultipleLocator(base=(max(vec) - min(vec)) / 10) # this locator puts ticks at regular intervals
    #ax.xaxis.set_major_locator(loc)
    #ax.yaxis.set_major_locator(loc)
    ax.set_xlabel(xlab)
    dir = " > " if cumulative == 1 else " < "
    ax.set_ylabel("cumulative fraction of reads" + dir + xlab)
    ax.set_title("number of reads " + str(len(vec)))
    fig.savefig(figure)
    plt.close()


class CdsStat:
    def __init__(self, sample_id):
        self.sample_id = sample_id
        self.total_frag =  0
        self.num_intermol = 0 # intermolecular
        self.num_overlap = 0
        self.read_len_diff = []
        self.aln_len_diff = []
        self.ol_ratios = []

    def plot_read_len_diff(self):
        cdf_plot(self.read_len_diff, 100, "read_len_diff.png", 1, "read len diff")

    def plot_ol_ratios(self):
        cdf_plot(self.ol_ratios, 100, "ol_ratios.png", -1, "overlap ratio")

    def __str__(self):
        return f"{self.sample_id}\t{self.total_frag}\t{self.num_intermol}\t{self.num_overlap}\t{self.num_intermol/self.total_frag}\t{self.num_overlap/self.total_frag}"



def process(options):
    samfile = pysam.AlignmentFile(options.bam, "rb")
    fid_size_dict = {}
    for chrom, start, end in read_bed(options.bed):
        for read in samfile.fetch(chrom, start, end):
            if read.is_unmapped or read.mate_is_unmapped:
                continue
            if read.is_read1 and read.has_tag('UG'):
                famid = read.get_tag('UG')
                intermol = False
                if not read.is_proper_pair or read.template_length == 0 or read.tempalte_length > options.im_dist_cutoff:
                    intermol = True
                if famid not in fid_size_dict:
                   fid_size_dict[famid] = [intermol, 1]
                else:
                    assert(intermol == fid_size_dict[famid][0])
                    fid_size_dict[famid][1] += 1
    keys = zip(*fid_size_dict.values)
    d = {'family_id' : list(fid_size_dict.keys()), 'is_intermol' : keys[0], 'size': keys[1]}


if __name__ == '__main__':
    sys.exit(process(get_arguments()))

