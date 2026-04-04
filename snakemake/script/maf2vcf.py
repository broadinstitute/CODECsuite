#!/usr/bin/env python

import argparse
import logging
import os
import sys
import pandas as pd
import subprocess as sp

logger = logging.getLogger("{}".format(__file__))

def get_arguments():

    parser = argparse.ArgumentParser(prog="convert maf to vcf file", formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("maf", type=str, help="input maf file")
    parser.add_argument("-r", "--ref_fasta", required=True, help="reference fasta")
    parser.add_argument("-p", "--path_to_maf2vcf", default = "/usr/bin/maf2vcf.pl" , type=str, help="path to the conversion script")
    parser.add_argument("-n", "--t_alt", type=int, default=1, help="min t_alt (1)")
    parser.add_argument("-O", "--out_type", choices=['vcf', 'maf'], default="vcf", help="output type")
    parser.add_argument("-o", "--out_dir", default="./", help="output dir")
    parser.add_argument("-s", "--sample_name", default="", help="sample name")
    parser.add_argument("-d", "--max_del_len", type=int, default=1e9, help="maximum deletion length allowed")
    parser.add_argument("-i", "--min_indel_len", type=int, default=1, help="min indel lngth allowed")
    args = parser.parse_args()
    return args


def process(options):
    maf = pd.read_csv(options.maf, comment="#", sep='\t', low_memory=False)
    print("input variant count:", maf.shape[0])
    if 't_alt_count' in maf.columns:
        maf = maf[maf['t_alt_count'] >= options.t_alt]
    if options.max_del_len < 1e9:
        maf = maf[maf['Reference_Allele'].str.len() <= options.max_del_len + 1]
    if options.min_indel_len > 1:
        maf = maf[abs(maf['Reference_Allele'].str.len() - maf['Tumor_Seq_Allele2'].str.len()) >= options.min_indel_len]
    print("output variant count:",  maf.shape[0])
    if not options.sample_name:
        sname = os.path.basename(options.maf)
        sname = sname[:-4]
    else:
        sname = options.sample_name
    outmaf = os.path.join(options.out_dir, sname+".filtered.maf")
    outtsv = os.path.join(options.out_dir, sname+".filtered.pairs.tsv")
    #maf = maf[["Chromosome", "Start_Position", "End_Position",
    #              "Variant_Type", "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele2", "t_alt_count", "t_ref_count", "n_alt_count", "n_ref_count"]]
    maf = maf[["Chromosome", "Start_Position", "Variant_Type", "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele2", "t_alt_count"]]
    maf.to_csv(outmaf, sep='\t', index=False)
    if options.out_type == "vcf":
        cmd=f"perl {options.path_to_maf2vcf} --input-maf {outmaf} --output-dir {options.out_dir} --ref-fasta {options.ref_fasta}"
        print(f"converting to vcf: {cmd}")
        sp.check_output(cmd, shell=True)
        sp.check_output(f"rm {outmaf}", shell=True)
        sp.check_output(f"rm {outtsv}", shell=True)


if __name__ == '__main__':
    sys.exit(process(get_arguments()))
