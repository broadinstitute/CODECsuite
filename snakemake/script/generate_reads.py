#!/usr/bin/env python

import argparse
import logging
import os
import sys
import pandas as pd
from Bio.Seq import reverse_complement

logger = logging.getLogger("{}".format(__file__))

"""
Simulate CDS reads from the xlxs file provided by Jin
"""
def get_arguments():

    parser = argparse.ArgumentParser(prog="foo", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("xlsx", type=str, help="input xlsx file")
    parser.add_argument("--readlen", type=int, default=300, help="read len")
    args = parser.parse_args()
    return args

def absolute_path(path):
    """convert relative path to absolute path"""
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(os.getcwd(), path)

def print_fastq(name, seq, filename):
    filename.write("@"+name)
    filename.write("\n")
    filename.write(seq)
    filename.write("\n")
    filename.write("+\n")
    filename.write("I" * len(seq))
    filename.write("\n")

def generate_read(template, readlen):
    read1 = None
    read2 = None
    if len(template) > readlen:
        ss = len(template) - readlen
        read1 = template[:readlen]
        read2 = reverse_complement(template[ss:])
    else:
        read1 = template
        read2 = reverse_complement(template)
    return (read1,read2)

def process(options):
    xl = pd.ExcelFile(options.xlsx)
    fastq1 = open("read1.fastq", "w")
    fastq2 = open("read2.fastq", "w")
    rl = options.readlen

    for sid, sheetname in enumerate(xl.sheet_names):
        sheet = xl.parse(sheetname)
        for cid, col in enumerate(sheet.columns):
            for rid, row in enumerate(sheet[col].tolist()):
                if rid != 0: continue
                row = row.strip(" '")
                n1 = "test_" + str(sid) + "_" + str(cid) + "_" + str(rid)
                n2 = "test_" + str(sid) + "_" + str(cid) + "_" + str(rid)
                read1, read2 = generate_read(row, rl)
                print_fastq(n1, read1, fastq1)
                print_fastq(n2, read2, fastq2)

                #row_rc = reverse_complement(row)
                #n1 = "test_rc" + str(sid) + "_" + str(cid) + "_" + str(rid)
                #n2 = "test_rc" + str(sid) + "_" + str(cid) + "_" + str(rid)
                #read1, read2 = generate_read(row_rc, rl)
                #print_fastq(n1, read1, fastq1)
                #print_fastq(n2, read2, fastq2)

if __name__ == '__main__':
    sys.exit(process(get_arguments()))
