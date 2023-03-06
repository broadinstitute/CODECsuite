#!/usr/bin/env python

import argparse
import sys
import os
import pysam
import random
from random import sample

def parse_cl_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--in_bam", help="input BAM", required=True)
    parser.add_argument("--sample_id", help="sample_id for series of bams", required=True)
    parser.add_argument("--outdir", help="output directory", default="./")
    parser.add_argument("--min_family_size", type = int, help="minimum family size (A+B) for downsampling", default=20)
    parser.add_argument("--max_family_size", type = int, help="max_family_size for output a bam", default=20)
    parser.add_argument("--min_strand_specific_family_size", type = int, default = 10)
    parser.add_argument("--seed", help="seed for random sampling", default=7)
    parser.add_argument("--is_codec", help="if is CODEC library", default=False, action='store_true')

    return parser.parse_args()

class PairEnd:
    def __init__(self, aln, is_cds):
        self.name = aln.query_name
        self.reads = []
        self.fid = None
        self.strand = None
        self.push(aln, is_cds)

    def push(self, aln, is_cds):
        assert(aln.query_name == self.name)

        try:
            mitag = aln.get_tag("MI")
        except KeyError:
            sys.stderr.write(aln.query_name + " has no MI tag\n")
        if is_cds:
            fid = mitag
            strand = None
        else:
            if "/" in mitag:
                fid, strand = mitag.split("/")
            else:
                fid = mitag
                strand = None
        if self.fid:
            assert(self.fid == fid)
        else:
            self.fid = fid
        if strand and self.strand:
            assert(strand == self.strand)
        elif strand:
            self.strand = strand

        self.reads.append(aln)


class Duplex:
    def __init__(self, pairend, is_cds):
        self.A_reads = []
        self.B_reads = []
        self.cds_reads = []
        self.fid = None
        self.push(pairend, is_cds)
        self.is_cds = is_cds

    def sizeA(self):
        if self.is_cds:
            return len(self.cds_reads)
        else:
            return len(self.A_reads)

    def sizeB(self):
        if self.is_cds:
            return len(self.cds_reads)
        else:
            return len(self.B_reads)

    def size(self):
        if self.is_cds:
            return len(self.cds_reads)
        else:
            return self.sizeA() + self.sizeB()

    def push(self, pairend, is_cds):
        if self.fid:
            assert(self.fid == pairend.fid)
        else:
            self.fid = pairend.fid
        if is_cds:
            self.cds_reads.append(pairend)
        else:
            if pairend.strand == 'A':
                self.A_reads.append(pairend)
            elif pairend.strand == 'B':
                self.B_reads.append(pairend)
            else:
                raise ValueError(pairend.name + " MI tag malformed\n")


def sub_sample(duplex, target_size, is_cds):
    assert(duplex.size() >= target_size)

    ret = []
    if is_cds:
        draw = sample(list(range(duplex.size())), target_size)
        ret = [duplex.cds_reads[ii] for ii in draw]
        return ret

    if target_size == 1:
        draw = sample(list(range(duplex.size())), target_size)
        ret.append(duplex.A_reads[draw[0]] if draw[0] < duplex.sizeA() else duplex.B_reads[draw[0] - duplex.sizeA()])
        return ret

    if duplex.sizeA() < duplex.sizeB():
        if duplex.sizeA() < target_size / 2:
            ret = duplex.A_reads
        else:
            idx = sample(list(range(duplex.sizeA())), int(target_size / 2))
            ret = [duplex.A_reads[ii] for ii in idx]
        rest_idx = sample(list(range(duplex.sizeB())), target_size - len(ret))
        ret = ret + [duplex.B_reads[ii] for ii in rest_idx]
    else:
        if duplex.sizeB() < target_size / 2:
            ret = duplex.B_reads
        else:
            idx = sample(list(range(duplex.sizeB())), int(target_size / 2))
            ret = [duplex.B_reads[ii] for ii in idx]
        rest_idx = sample(list(range(duplex.sizeA())), target_size - len(ret))
        ret = ret + [duplex.A_reads[idx] for idx in rest_idx]

    return ret


class DuplexFamilyBamsWriter:
    def __init__(self, in_bam, max_fs, min_fs, min_sp_fs, outdir, sid, is_codec):
        self.max_family_size = max_fs
        self.min_family_size_to_downsample = min_fs
        self.min_strand_specific_size_to_downsample = min_sp_fs
        self.sample_id = sid
        self.out_bams = []
        self.is_codec = is_codec
        for ii in range(max_fs):
            fname = sid + "_" + str(ii + 1) + ".bam"
            fpath = os.path.join(outdir, fname)
            self.out_bams.append(pysam.AlignmentFile(fpath, "wb", template=in_bam))

    def write_duplex(self, duplex):
        if duplex.size() >= self.min_family_size_to_downsample and  \
            duplex.sizeA() >= self.min_strand_specific_size_to_downsample and \
                duplex.sizeB() >= self.min_strand_specific_size_to_downsample:
            for fs in range(self.max_family_size):
                if duplex.size() > fs:
                    fragments = sub_sample(duplex, fs + 1, self.is_codec)
                    for frag in fragments:
                        for record in frag.reads:
                            self.out_bams[fs].write(record)

    def close_all(self):
        for bam in self.out_bams:
            bam.close()


def process(opts):

    random.seed(opts.seed)
    in_bam = pysam.AlignmentFile(opts.in_bam, "rb")
    dstack = [] #duplex stack
    fstack = [] #fragment stack
    writer = DuplexFamilyBamsWriter(in_bam, opts.max_family_size, opts.min_family_size, opts.min_strand_specific_family_size,
                                    opts.outdir, opts.sample_id, opts.is_codec)

    for aln in in_bam.fetch(until_eof=True):
        if fstack:
            if fstack[-1].name == aln.query_name:
                fstack[-1].push(aln, is_cds=opts.is_codec)
                continue
            else:
                pe = fstack.pop()
                # do something
                if dstack:
                    if dstack[-1].fid == pe.fid:
                        dstack[-1].push(pe, is_cds=opts.is_codec)
                    else:
                        dpx = dstack.pop()
                        writer.write_duplex(dpx)
                        dpx = Duplex(pe, is_cds=opts.is_codec)
                        dstack.append(dpx)
                else:
                    dpx = Duplex(pe, is_cds=opts.is_codec)
                    dstack.append(dpx)

        pe = PairEnd(aln, opts.is_codec)
        fstack.append(pe)

    # process the last read
    if fstack:
        pe = fstack.pop()
        # do something
        if dstack:
            if dstack[-1].fid == pe.fid:
                dstack[-1].push(pe, is_cds = opts.is_codec)
            else:
                dpx = dstack.pop()
                writer.write_duplex(dpx)
                dpx = Duplex(pe, is_cds=opts.is_codec)
                dstack.append(dpx)

    if dstack:
        dpx = dstack.pop()
        writer.write_duplex(dpx)

    in_bam.close()
    writer.close_all()


if __name__ == "__main__":
    process(parse_cl_args())
