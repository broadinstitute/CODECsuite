#!/usr/bin/env python
import argparse
import logging
import pandas as pd
import sys

logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="aggregate miredas reult(s)", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("log", type=str, nargs="+", help="trim adapter logs")
    parser.add_argument("out", type=str, help="output name")
    args = parser.parse_args()
    return args

def process(opts):
    tot = pd.DataFrame()
    for f in opts.log:
        df = pd.read_csv(f, sep=":", names=['cat', "count"])
        if tot.empty:
            tot = df
        else:
            tot['count'] = tot['count'].add(df['count'])
    tot['cat'] = tot['cat'] + ":"
    tot.to_csv(opts.out, sep=" ", index=False, header=False)

if __name__ == '__main__':
    sys.exit(process(get_arguments()))
