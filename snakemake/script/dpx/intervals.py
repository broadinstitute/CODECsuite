from quicksect import IntervalTree
from collections import defaultdict
from math import ceil
def read_bed(bedfile):
    """ Creates generator from bed file or interval_list """
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


def create_interval_dict_from_bed(bedfile, midpoint=False):
    """
    Used for marking on/off target fragments by creating interval
    trees for each chromosome in bedfile.
    """
    interval_dict = defaultdict(IntervalTree)
    if bedfile:
        for chrom, start, stop in read_bed(bedfile):
            if midpoint:
                mid = ceil((start+stop)/2)
                start = mid - 1
                stop = mid
            interval_dict[str(chrom)].add(start, stop)
    return interval_dict
