#!/usr/bin/env python

import pandas as pd
import sys
import pysam
import subprocess
import argparse
import os

def get_arguments():

    parser = argparse.ArgumentParser(prog="grep false positive reads and print from gbu bam. Additionaly can print fastq files",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fq1", type=str, help="raw fastq1 file", default="")
    parser.add_argument("--fq2", type=str, help="raw fastq2 file", default="")
    parser.add_argument("bam", type=str, help="groupby_umi_bam")
    parser.add_argument('df', type=str, help="error family file")
    parser.add_argument('prefix', type=str, help="prefix for output")

    args = parser.parse_args()
    return args

def process(options):
    prefix_id = options.prefix

    df = pd.read_csv(options.df, sep="\t", low_memory=False)
    new = df['family_id'].str.split('/', n = 1, expand = True)
    df['family_num']  = new[0]
    df['strand'] = new[1]
    df = df.astype({'family_num' : 'int64'})
    df = df.sort_values(by=['family_num', 'strand'])
    df = df[df['pass_filter'] == 1]
    gbubam = pysam.AlignmentFile(options.bam, "rb")
    tmpbam = pysam.AlignmentFile(prefix_id + ".tmp.bam", "wb", template=gbubam)
    f = open(prefix_id + ".run", 'w')
    rowit = df.iterrows()
    cur_row = next(rowit)[1]
    in_match = False
    counter = 0
    for aln in gbubam.fetch(until_eof=True):
        mitag = aln.get_tag("MI")
        if mitag == cur_row['family_id']:
            if not in_match:
                print(cur_row['family_id'])
                in_match = True
                counter = 0
            if options.fq1:
                if aln.is_read1:
                    counter += 1
                    cmd = ["zegrep", "-A 3", "-m 1", aln.query_name, options.fq1, ">>", prefix_id + ".1.fastq #" + mitag]
                    print(" ".join(cmd), file=f)
            if options.fq2:
                if aln.is_read2:
                    cmd = ["zegrep", "-A 3", "-m 1", aln.query_name, options.fq2, ">>", prefix_id + ".2.fastq #" + mitag]
                    print(" ".join(cmd), file=f)
            tmpbam.write(aln)
        if in_match and mitag != cur_row['family_id']:
            print("matched", counter, "times")
            in_match = False
            try:
                cur_row = next(rowit)[1]
            except StopIteration:
                break
    tmpbam.close()
    cmd = "samtools sort {0} -o {1} && samtools index {1} && rm {0}".format(prefix_id + ".tmp.bam", prefix_id + ".bam")
    subprocess.check_call(cmd, shell=True)
if __name__ == '__main__':
    sys.exit(process(get_arguments()))
