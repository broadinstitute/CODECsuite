#!/usr/bin/env python

import argparse
import logging
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.ticker import MaxNLocator
import seaborn as sns
sns.set(font_scale=2)
logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="plot error rates stratified by family size", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("accu", type=str, help="input accu file")
    args = parser.parse_args()
    return args

def absolute_path(path):
    """convert relative path to absolute path"""
    if os.path.isabs(path):
        return path
    else:
        return os.path.join(os.getcwd(), path)

def process(options):
    fsize_accu = {}
    with open(options.accu, 'r') as fh:
        for line in fh:
            cols = line.strip().split('\t')
            if cols[1] == "chrY" or cols[1] == "chrX":
                continue
            fields = cols[0].split('_')
            fsize = int(fields[-1])
            if fsize not in fsize_accu:
                fsize_accu[fsize] = [int(cols[3]), int(cols[4]), 0]
            else:
                fsize_accu[fsize][0] += int(cols[3])
                fsize_accu[fsize][1] += int(cols[4])
                fsize_accu[fsize][2] += 1

    accus = [0] * 10
    counts = [0] * 10
    for k,v in fsize_accu.items():
        if k <= 10:
            print (k, v[0], v[1], v[1]/v[0], v[2])
            accus[k-1] = v[1]/v[0] * 1e5
            counts[k-1] = v[2] * 1e-6
    fig, ax = plt.subplots(figsize=(11.7, 8.27))
    ax.bar(np.arange(10) + 1, accus)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel("Family size")
    plt.ylabel("Error per 100k")
    fig.savefig('error_rate_by_family_size.png')

    fig1, ax1 = plt.subplots(figsize=(11.7, 8.27))
    ax1.bar(np.arange(10) + 1, counts)
    ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
    #ax1.ticklabel_format(axis='y', useMathText=True)
    plt.xlabel("Family size")
    plt.ylabel("Million Fragment")
    fig1.savefig('num_frag_family_size.png')

if __name__ == '__main__':
    sys.exit(process(get_arguments()))

