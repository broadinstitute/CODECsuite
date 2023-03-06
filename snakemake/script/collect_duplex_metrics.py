import logging

import click
import numpy as np
import pysam
import math
from collections import defaultdict
import sys
from random import sample, seed

from dpx.bam_iterator import iterate_read_pairs, get_unclipped_fragment_ends, is_proper_orientation, read_generator
from dpx.downsampler import Downsampler, get_overlap
from dpx.intervals import create_interval_dict_from_bed


seed(10)
logger = logging.getLogger("CollectDuplexMetrics")

def get_family_tag(alignment):
    return alignment.get_tag("MI") if alignment.has_tag("MI") \
        else alignment.query_name.split(":")[1] \
        if ":" in alignment.query_name else "null"


class ReadPair:
    def __init__(self, read1, read2):
        self.read1 = read1
        self.read2 = read2
        readrepr = None
        if read1:
            readrepr = read1
        else:
            readrepr = read2
        if readrepr.has_tag("MI"):
            self.family_tag = readrepr.get_tag("MI")
            self.family = int(self.family_tag.split("/")[0])
            self.strand = self.family_tag[-1]
        # elif readrepr.has_tag("UG"):
        #     self.family = readrepr.get_tag("UG")
        #     self.strand = "B" if readrepr.is_reverse else "A"
        else:
            self.family = None
            self.strand = None
        self._set_coordinates()

    def __str__(self):
        attributes = ["coordinate_id", "family", "strand", "target", "overlap"]
        attr_dict = {attr: getattr(self, attr) for attr in attributes}
        return str(attr_dict)

    def get_tag(self, tag):
        return self.read1.get_tag(tag)

    def _set_coordinates(self):
        start, end = get_unclipped_fragment_ends(
            read1 = self.read1, read2 = self.read2
        )
        self.chrom = self.read1.reference_name
        self.start = start
        self.end = end
        self.insert_size = self.end - self.start

    def are_ends_overlapped(self):
        return self.read1 and self.read2 and \
               self.read1.reference_name == self.read2.reference_name and \
               self.read1.reference_start < self.read2.reference_end and \
               self.read1.reference_end > self.read2.reference_start

    def get_overlap(self, interval_dict):
        """ Get target and overlap with given interval dict """
        r1s = self.read1.reference_start
        r1e = self.read1.reference_end
        r2s = self.read2.reference_start
        r2e = self.read2.reference_end
        chrom = str(self.chrom)
        r1res = None
        r2res = None
        if self.read1.reference_name in interval_dict:
            r1res = interval_dict[self.read1.reference_name].search(r1s, r1e)
        if self.read2.reference_name in interval_dict:
            r2res = interval_dict[self.read2.reference_name].search(r2s, r2e)
        common_iv = None
        if r1res:
            for iv in r1res:
                if iv in r2res:
                    common_iv = iv
                    break
        if common_iv:
            start, end = self.start, self.end
            tstart, tend = common_iv.start, common_iv.end
            overlap = min(end, tend) - max(start, tstart)
            target = f"{chrom}:{tstart}-{tend}"
            return target, overlap
        if bool(r1res) != bool(r2res):
            if r1res:
                start = self.read1.reference_start
                end = self.read1.reference_end
                tstart, tend = r1res[0].start, r1res[0].end
                overlap = min(end, tend) - max(start, tstart)
                target = f"{self.read1.reference_name}:{tstart}-{tend}"
            if r2res:
                start = self.read2.reference_start
                end = self.read2.reference_end
                tstart, tend = r2res[0].start, r2res[0].end
                overlap = min(end, tend) - max(start, tstart)
                target = f"{self.read2.reference_name}:{tstart}-{tend}"
            return target, overlap
        return None, 0

    @property
    def coordinate_id(self):
        """ Region string of fragment """
        return f"{self.chrom}:{self.start}-{self.end}"

def generate_chunks(
    bam, records, interval_dict, insert_sizes, min_overlap=1, per_target=False
):
    """ Create chunks every *records* read pairs where all reads sharing same
    start/stop positions are kept together. If per_target = True, rather than
    grouping by start/stop position of read, it will group reads by
    start/stop position of the probe.

    TODO This logic should be cleaned up...
    """
    min_insert, max_insert = insert_sizes
    current_id = ""
    current_target = ""
    current_bucket = []
    res = []

    logger.info("Begin iterating over BAM file")

    for read1, read2 in iterate_read_pairs(bam):
        if not is_proper_orientation(read1 = read1, read2 = read2):
            continue
        read_pair = ReadPair(read1, read2)
        # Exclude if less than min_overlap with probe
        iter_id = read_pair.family
        insert_size = read_pair.insert_size
        if not (min_insert <= insert_size <= max_insert):
            continue
        if per_target:
            target, overlap = read_pair.get_overlap(interval_dict)
            if overlap < min_overlap:
                continue
            if not target:
                continue
            if target != current_target:
                if current_bucket:
                    yield current_bucket
                current_bucket = [read_pair]
                current_target = target
                current_id = iter_id
            else:
                if len(current_bucket) >= records:
                    # Keep adding to bucket if same start/stop
                    if iter_id == current_id:
                        current_bucket.append(read_pair)
                    else:
                        yield current_bucket
                        current_bucket = [read_pair]
                        current_id = iter_id
                else:
                    current_bucket.append(read_pair)
        else:
            if iter_id != current_id:
                if len(res) + len(current_bucket) <= records:
                    res.extend(current_bucket)
                    current_bucket = [read_pair]
                else:
                    yield res
                    res = current_bucket
                    res.append(read_pair)
                    current_bucket = []
                current_id = iter_id
            else:
                current_bucket.append(read_pair)
    # If there is anything left, yield it
    if res or current_bucket:
        res.extend(current_bucket)
        yield res
    logger.info(f"Processed {iterate_read_pairs.counts:,} read pairs.")

def generate_chunks_cds(
        bam, batch_size, interval_dict, insert_sizes, min_overlap=1, per_target=False):

    """
    TDDO: should support both umi tools and fgbio groupbyumi output
    """
    min_insert, max_insert = insert_sizes
    current_target = ""
    current_bucket = []
    read_dict = defaultdict(lambda: [None, None])

    logger.info("Begin iterating over BAM file")
    res = []
    for read1, read2 in iterate_read_pairs(bam):
        if not is_proper_orientation(read1 = read1, read2 = read2):
            continue

        # Exclude if less than min_overlap with probe
        target, overlap = get_overlap(read, interval_dict)
        insert_size = abs(read.tlen)
        if overlap < min_overlap:
            continue
        if not (min_insert <= insert_size <= max_insert):
            continue
        if not target:
            continue
        ## yield readpair use read_dict
        qname = read.query_name
        read_pair = None
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                read_pair = ReadPair(read, read_dict[qname][1])
            else:
                read_pair = ReadPair(read_dict[qname][0], read)
            del read_dict[qname]

        if per_target:
            if target != current_target:
                if current_target:
                    # for unpaired reads
                    for qname, preads in read_dict.items():
                        itm_read_pair = ReadPair(preads[0], preads[1])
                        current_bucket.append(itm_read_pair)
                    read_dict.clear()
                    yield current_bucket
                if read_pair:
                    current_bucket = [read_pair]
                current_target = target
            else:
                if read_pair:
                    current_bucket.append(read_pair)
        else:
            if read_pair:
                current_bucket.append(read_pair)
            if current_target and target != current_target:
                if len(current_bucket) > batch_size:
                    for qname, preads in read_dict.items():
                        itm_read_pair = ReadPair(preads[0], preads[1])
                        current_bucket.append(itm_read_pair)
                    read_dict.clear()
                    yield current_bucket
                    del current_bucket[:]
            current_target = target

    if current_bucket or read_dict:
        for qname, preads in read_dict.items():
            read_pair = ReadPair(preads[0], preads[1])
            current_bucket.append(read_pair)
        yield current_bucket

def create_row(counts, info, cols):
    """ Create row from downsampling counts dictionary """
    res = info
    for col in cols:
        if col == "ds_fraction_duplexes":
            res.append(
                counts.get("ds_duplexes", 0) / counts.get("ds_families", 1)
            )
        else:
            res.append(counts.get(col, 0))
    return res


def write_metrics_to_file(downsampler, per_target, output):
    """ Write downsampler results to file mimicking Fgbio output """
    logger.info(f"Writing results to {output}")
    with open(output, "w") as outfile:
        data_cols = [
            "read_pairs",
            "cs_families",
            "ss_families",
            "ds_families",
            "ds_duplexes",
            "ds_fraction_duplexes",
        ]
        info = ["fraction"]
        if per_target:
            info = ["target", "fraction"]
        headers = info + data_cols
        outfile.write("\t".join(headers) + "\n")
        if per_target:
            for target, downsample in downsampler.counts.items():
                for downsampling, counts in downsample.items():
                    row = create_row(counts, [target, downsampling], data_cols)
                    outfile.write("\t".join(map(str, row)) + "\n")

        else:
            for downsampling, counts in downsampler.counts.items():
                row = create_row(counts, [downsampling], data_cols)
                outfile.write("\t".join(map(str, row)) + "\n")
    return None


class FloatList(click.ParamType):
    name = "floatlist"

    def convert(self, value, param, ctx):
        try:
            return list(map(float, value.split(",")))
        except ValueError:
            self.fail(
                f"{value} is not a comma separated float list", param, ctx,
            )


@click.command(context_settings=dict(help_option_names=["-h", "--help"]))
@click.option(
    "-b",
    "--bam_file",
    type=click.Path(exists=True),
    help="BAM file from Fgbio GroupReadsByUmi",
    required=True,
)
@click.option(
    "-o",
    "--output_file",
    type=click.Path(writable=True),
    help="File to write output",
    required=True,
)
@click.option(
    "-l",
    "--interval_list",
    type=click.Path(exists=True),
    help="Interval list or BED file with baited regions",
)
@click.option(
    "-m",
    "--min_overlap",
    type=click.INT,
    default=1,
    help="Minimum overlap with bait regions to accept read pair fragment",
    show_default=True,
)
@click.option(
    "-d",
    "--downsampling",
    help="Comma separated downsampling [default: np.logspace(-4, 0, 30)]",
    default=",".join(map(str, np.logspace(-4, 0, 30))),
    type=FloatList(),
)
@click.option(
    "-i",
    "--insert_sizes",
    type=FloatList(),
    help="Comma separated insert size range [default: [0, inf]",
    default="0,inf",
)
@click.option(
    "--min_reads",
    type=click.INT,
    default=2,
    help="Minimum number of reads for both strand to consider a duplex. min(aD, bD) >= this",
    show_default=True,
)
@click.option(
    "--min_max_strand_reads",
    type=click.INT,
    default=0,
    help="max(aD, bD) >= this",
    show_default=True,
)
@click.option(
    "-p",
    "--per_target",
    is_flag=True,
    help="Output metrics on a per-target level",
)
@click.option(
    "-c",
    "--is_cds",
    is_flag=True,
    help="Is CDS library",
)

@click.option(
    "-r",
    "--reduce_inv",
    is_flag=True,
    help="reduce target interval to just the mid-point",
)

def main(
    bam_file,
    output_file,
    interval_list,
    min_overlap,
    downsampling,
    insert_sizes,
    min_reads,
    min_max_strand_reads,
    per_target,
    is_cds,
    reduce_inv,
):
    """ Tool mimicking Fgbio's CollectDuplexSeqMetrics but with additional features.
    With this tool, you can output metrics on a per-target level and set probe overlap
    requirements. """

    if not interval_list:
        min_overlap = 0
        interval_dict = {}
    else:
        logger.info("Creating interval tree from region file")
        interval_dict = create_interval_dict_from_bed(interval_list, reduce_inv)
        logger.info("Finished creating interval tree")
    bam = pysam.AlignmentFile(bam_file, "rb")
    probabilities = downsampling
    chunk_size = 100_000

    downsampler = Downsampler(
        probabilities=probabilities,
        min_min_strand_reads=min_reads,
        min_max_strand_reads=min_max_strand_reads,
        per_target=per_target,
        interval_dict=interval_dict,
        is_cds = is_cds
    )
    if is_cds:
        for chunk in generate_chunks_cds(
            bam,
            chunk_size,
            insert_sizes=insert_sizes,
            min_overlap=min_overlap,
            interval_dict=interval_dict,
            per_target=per_target
        ):
            downsampler.run_downsamplings(chunk)
    else:
        for chunk in generate_chunks(
            bam,
            chunk_size,
            insert_sizes=insert_sizes,
            min_overlap=min_overlap,
            interval_dict=interval_dict,
            per_target=per_target
        ):
            downsampler.run_downsamplings(chunk)
    write_metrics_to_file(downsampler, per_target, output_file)


if __name__ == "__main__":
    main()
