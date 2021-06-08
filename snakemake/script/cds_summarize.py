#!/usr/bin/env python
import argparse
import logging
import sys
import json
import os
import pandas as pd
import numpy as np
import pysam
from collections import defaultdict
#from Bio import pairwise2

logger = logging.getLogger("{}".format(__file__))

#def check_tandem_adpt(seq):
    #linker = "AGATCGGAAGAGCTTCATCATTAGATCCATTAATGTTACACTTCAACTCTTCACCCACATCAGATTAGTACCAGCTTCGAGGATCAACACGTCAGAGTCTAGCTGGTGATAGGAAGTGTAGGTAACATAGACGAAGTTATCAACAATGTGTAACTGACTTAACGCTCTTCCGATCT"
    #res = pairwise2.align.localms(seq, linker, 1, -4, -6,-2)

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

def is_overlapped(read1, read2):
    if read1.is_unmapped or read2.is_unmapped:
        return False
    if read1.reference_name != read2.reference_name:
        return False
    if read1.reference_start < read2.reference_end and read2.reference_start < read1.reference_end:
        return True

def is_complete_overlapped_excluding_sclips(read1, read2):
    if not is_overlapped(read1, read2):
        return False
    if read1.reference_start != read2.reference_start:
        return False
    if read1.reference_end != read2.reference_end:
        return False
    return True

def overlap_len(read1, read2):
    if not is_overlapped(read1,read2):
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

    parser = argparse.ArgumentParser(prog="Parse Fastp json result(s)", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fastp", type=str, help="json file output by fastp", required=False)
    parser.add_argument("--trim_log", type=str, default = "", help="trim linker log file", required=False)
    parser.add_argument("--highconf_bam", type=str, default = "", help="high confident CDS reads", required=False)
    parser.add_argument("--lowconf_bam", type=str, default = "", help="low confident CDS reads", required=False)
    parser.add_argument("--si_hiconf_bam", type=str, default = "", help="single insert highconf bam", required=False)
    parser.add_argument("--si_lowconf_bam", type=str, default = "", help="single insert lowconf bam", required=False)
    parser.add_argument("--trim_one_bam", type=str, default = "", help="Where linker has been trimmed from only one end", required=False)
    parser.add_argument("--untrim_both_bam", type=str, default = "", help="Where linker has not been trimmed from both ends", required=False)
    parser.add_argument("--hs_metrics", type=str, help="hs metrics output", default="")
    parser.add_argument("--sample_id", type=str, help="sample id", required=True)
    parser.add_argument("--cds_intermol_bamout", type=str, help="output bam for cds intermolcular reads", default = "", required=False)
    parser.add_argument("--demux_log", type=str, help="demux log file", default = "", required=False)
    #parser.add_argument("--si_intermol_bamout", type=str, help="output bam for singleinsert intermolcular reads", required=True)
    args = parser.parse_args()
    return args

class CdsMetrics:
    """
    ## uf = unfiltered
    ## pf = passed filter
    """
    def __init__(self, sid):
        self.sample_id = sid
        self.n_raw_frag = 0
        self.n_unmapped = 0
        self.pct_raw_uf_q20 = 0
        self.pct_raw_uf_q30 = 0
        self.n_raw_pf_frag = 0
        self.pct_aligned_frag = 0
        self.on_target_rate = 0
        self.mean_bait_cov = 0
        self.mean_target_cov = 0

#CDS specific
        self.n_high_conf = 0
        self.n_adp_dimer_frag = 0
        self.n_double_ligation = 0
        self.n_low_conf = 0
        self.n_intermol = 0 # intermolecular
        self.n_single_hiconf = 0
        self.n_single_lowconf = 0
        self.n_insuf_trim = 0
        self.n_close_proxim = 0
        # self.read_len_diff = []
        # self.aln_len_diff = []
        # self.ol_ratios = defaultdict(list)

    def n_raw_uf_frag(self):
        return self.n_raw_frag
    def n_categorized(self):
        return self.n_high_conf + self.n_low_conf + self.n_adp_dimer_frag + self.n_double_ligation + \
               self.n_intermol + self.n_single_hiconf + self.n_single_lowconf + self.n_insuf_trim + self.n_unmapped + self.n_close_proxim
    def __str__(self):
        header = ["sample_id",
                  ### CDS specific
                  "pct_cds_highconf",
                  "pct_cds_lowconf",
                  "pct_single_lowconf",
                  "pct_single_highconf",
                  "pct_double_ligation",
                  "pct_insuf_trim",
                  "pct_adp_dimer",
                  "pct_intermol",
                  "pct_unmapped",
                  "pct_close_proxim",
                  "pct_categorized",
                  "n_cds_highconf",
                  "n_cds_lowconf",
                  "n_single_lowconf",
                  "n_single_highconf",
                  "n_double_ligation",
                  "n_insuf_trim",
                  "n_adp_dimer",
                  "n_intermol",
                  "n_unmapped",
                  "n_close_proxim",
                  "n_categorized",
                  "n_total",
                  #####
                  "pct_aligned_frag",
                  "on_target_rate",
                  "mean_bait_cov",
                  "mean_target_cov",
                  "n_pass_fastp_filter",
                  "pct_pass_fastp_filter",
                  "pct_raw_uf_q20",
                  "pct_raw_uf_q30",
                  ]

        header_str = "\t".join(header)
        return f"{header_str}\n" \
               f"{self.sample_id}\t" \
               f"{self.n_high_conf / self.n_raw_uf_frag()}\t" \
               f"{self.n_low_conf / self.n_raw_uf_frag()}\t" \
               f"{self.n_single_lowconf / self.n_raw_uf_frag()}\t" \
               f"{self.n_single_hiconf / self.n_raw_uf_frag()}\t" \
               f"{self.n_double_ligation / self.n_raw_uf_frag()}\t" \
               f"{self.n_insuf_trim / self.n_raw_uf_frag()}\t" \
               f"{self.n_adp_dimer_frag / self.n_raw_uf_frag()}\t" \
               f"{self.n_intermol / self.n_raw_uf_frag()}\t" \
               f"{self.n_unmapped / self.n_raw_uf_frag()}\t" \
               f"{self.n_close_proxim / self.n_raw_uf_frag()}\t" \
               f"{self.n_categorized() / self.n_raw_uf_frag()}\t" \
               f"{self.n_high_conf}\t" \
               f"{self.n_low_conf}\t" \
               f"{self.n_single_lowconf}\t" \
               f"{self.n_single_hiconf}\t" \
               f"{self.n_double_ligation}\t" \
               f"{self.n_insuf_trim}\t" \
               f"{self.n_adp_dimer_frag}\t" \
               f"{self.n_intermol}\t" \
               f"{self.n_unmapped}\t" \
               f"{self.n_close_proxim}\t" \
               f"{self.n_categorized()}\t" \
               f"{self.n_raw_frag}\t" \
               f"{self.pct_aligned_frag}\t" \
               f"{self.on_target_rate}\t" \
               f"{self.mean_bait_cov}\t" \
               f"{self.mean_target_cov}\t" \
               f"{self.n_raw_pf_frag}\t" \
               f"{self.n_raw_pf_frag / self.n_raw_frag}\t" \
               f"{self.pct_raw_uf_q20}\t" \
               f"{self.pct_raw_uf_q30}"

def parse_linker_trim_log(log_file, cdsm, adap_v2):
    with open(log_file, 'r') as f:
        for line in f:
            k, v = line.split(":")
            if k == "TOTAL" or k == "TOTOL":
                cdsm.n_raw_frag = int(v)
            if k == "LOST_BOTH":
                cdsm.n_adp_dimer_frag = int(v)
            elif not adap_v2 and k == "DOUBLE_LIGATION":
                cdsm.n_double_ligation += int(v)
            elif adap_v2 and k== "LOST_READ1":
                cdsm.n_double_ligation += int(v)
            elif adap_v2 and k== "LOST_READ2":
                cdsm.n_double_ligation += int(v)


def alignment_analysis(bam, cdsm, trim_type, intermol_bam = None, im_dist_cutoff = 5_000, adap_v2 = False):
    assert(trim_type in ['HighConf', "LowConf", "TrimOne", "UntrimBoth", "SingleHiconf", "SingleLowconf"])
    samfile = pysam.AlignmentFile(bam, "rb")
    total_frag = 0
    for read1, read2 in read_pair_generator(samfile):
        total_frag += 1
        if read1.is_unmapped or read2.is_unmapped:
            cdsm.n_unmapped += 1
            continue

        if read1.reference_name != read2.reference_name or abs(read1.tlen) > im_dist_cutoff:
            cdsm.n_intermol += 1
            if intermol_bam:
                intermol_bam.write(read1)
                intermol_bam.write(read2)
            continue

        if trim_type == "HighConf":
            if adap_v2:
                if is_overlapped(read1, read2):
                    cdsm.n_high_conf += 1
                else:
                    cdsm.n_close_proxim += 1
            else:
                if is_complete_overlapped_excluding_sclips(read1, read2):
                    cdsm.n_high_conf += 1
        elif trim_type == "LowConf":
            if is_complete_overlapped_excluding_sclips(read1, read2):
                cdsm.n_low_conf += 1
        elif trim_type == "SingleHiconf":
            if (not read1.has_tag('tm') or read1.get_tag('tm') != 4 ) and \
                    (not read2.has_tag('tm') or read2.get_tag('tm') != 4 ) and \
                    is_complete_overlapped_excluding_sclips(read1, read2):
                    cdsm.n_single_hiconf += 1
        elif trim_type == "SingleLowconf":
            if is_complete_overlapped_excluding_sclips(read1, read2):
                cdsm.n_single_lowconf += 1
        elif trim_type == "TrimOne" or trim_type == "UntrimBoth":
            if not adap_v2 and is_overlapped(read1, read2):
                cdsm.n_insuf_trim += 1

    return total_frag

def process(opts):
    cdsm = CdsMetrics(opts.sample_id)
    adap_v2 = False
    if opts.fastp:
        with open(opts.fastp, 'r') as f:
            sample_dict = json.load(f)
            cdsm.n_raw_frag  = int(sample_dict['summary']["before_filtering"]["total_reads"])/2
            cdsm.n_raw_pf_frag = int(sample_dict['summary']['after_filtering']['total_reads'])/2
            cdsm.pct_raw_uf_q20 = sample_dict['summary']['before_filtering']['q20_rate']
            cdsm.pct_raw_uf_q30 = sample_dict['summary']['before_filtering']['q30_rate']

    if opts.cds_intermol_bamout:
        bam = pysam.AlignmentFile(opts.highconf_bam, "rb")
        cds_intermol_writer = pysam.AlignmentFile(opts.cds_intermol_bamout + ".tmp.bam", "wb", template=bam)
        bam.close()

    if opts.lowconf_bam:
        num_frag_processed = alignment_analysis(opts.lowconf_bam, cdsm, trim_type = 'LowConf')
        adap_v2 = True if num_frag_processed == 0 else False
    else:
        adap_v2 = True
    if opts.highconf_bam:
        if opts.cds_intermol_bamout:
            alignment_analysis(opts.highconf_bam, cdsm, trim_type = 'HighConf', intermol_bam=cds_intermol_writer, adap_v2=adap_v2)
        else:
            alignment_analysis(opts.highconf_bam, cdsm, trim_type='HighConf', adap_v2=adap_v2)
    if opts.trim_one_bam:
        alignment_analysis(opts.trim_one_bam, cdsm, trim_type = 'TrimOne')
    if opts.untrim_both_bam:
        alignment_analysis(opts.untrim_both_bam, cdsm, trim_type = 'UntrimBoth')
    if opts.si_hiconf_bam:
        alignment_analysis(opts.si_hiconf_bam, cdsm, trim_type = 'SingleHiconf')
    if opts.si_lowconf_bam:
        alignment_analysis(opts.si_lowconf_bam, cdsm, trim_type = 'SingleLowconf')

    if opts.trim_log:
        parse_linker_trim_log(opts.trim_log, cdsm, adap_v2)

    if opts.hs_metrics:
        hs_metrics_df = pd.read_csv(opts.hs_metrics, skiprows=6, nrows=1, sep='\t', low_memory=False)
        cdsm.pct_aligned_frag = hs_metrics_df['PCT_PF_UQ_READS_ALIGNED'][0]
        cdsm.on_target_rate = hs_metrics_df['PCT_SELECTED_BASES'][0]
        cdsm.mean_bait_cov = hs_metrics_df['MEAN_BAIT_COVERAGE'][0]
        cdsm.mean_target_cov = hs_metrics_df['MEAN_TARGET_COVERAGE'][0]
    print(cdsm)
    if opts.cds_intermol_bamout:
        cds_intermol_writer.close()
        os.system(f"samtools sort -n {opts.cds_intermol_bamout}.tmp.bam -o {opts.cds_intermol_bamout} && rm {opts.cds_intermol_bamout}.tmp.bam")


if __name__ == '__main__':
    sys.exit(process(get_arguments()))
