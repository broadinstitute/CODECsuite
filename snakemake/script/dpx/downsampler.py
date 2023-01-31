from collections import defaultdict
from collections import Counter

import numpy as np
import sys


np.random.seed(10)

def get_overlap(read, interval_dict):
    rchrom = str(read.reference_name)
    rstart = read.reference_start
    rend = read.reference_end
    if rchrom in interval_dict:
        olap = interval_dict[rchrom].search(rstart, rend)
        if olap:
            tstart, tend = olap[0].start, olap[0].end
            overlap = min(rend, tend) - max(rstart, tstart)
            target = f"{rchrom}:{tstart}-{tend}"
            return target, overlap
    return None, 0

class Downsampler:
    """ Class to perform downsampling on a list of family IDs
    from Fgbio GroupReadsByUmi """

    def __init__(self, probabilities, min_min_strand_reads, min_max_strand_reads, per_target, interval_dict, is_cds):
        self.probabilities = probabilities
        self.min_min_strand_reads = min_min_strand_reads
        self.min_max_strand_reads = min_max_strand_reads
        self.kept_families = defaultdict(list)
        self.counts = (
            defaultdict(Counter)
            if not per_target
            else defaultdict(lambda: defaultdict(Counter))
        )
        self.per_target = per_target
        self.interval_dict = interval_dict
        self.is_cds = is_cds

    def downsample(self, read_pairs, probability):
        """ Downsamples list of read pairs at a given probability """
        duplexes = defaultdict(lambda: defaultdict(lambda: 0))
        summary_counts = defaultdict(lambda: 0)
        min_min_strand_reads = self.min_min_strand_reads
        min_max_strand_reads = self.min_max_strand_reads
        kept_reads = []
        previous_coordinate = set()
        for read_pair in read_pairs:
            strand = read_pair.strand
            family = read_pair.family
            if not family:
                continue
            coordinate_id = read_pair.coordinate_id
            if np.random.random() <= probability:
                summary_counts["read_pairs"] += 1
                if coordinate_id not in previous_coordinate:
                    summary_counts["cs_families"] += 1
                    previous_coordinate.add(coordinate_id)
                duplexes[family][strand] += 1
                if self.is_cds and read_pair.are_ends_overlapped():
                    if strand == "A":
                        duplexes[family]["B"] += 1
                    elif strand == "B":
                        duplexes[family]["A"] += 1
                #if self.per_target:
                kept_reads.append(read_pair)
        for family, count in duplexes.items():
            ss_families = (count["A"] > 0) + (count["B"] > 0)
            ds_families = int(ss_families > 0)
            summary_counts["ss_families"] += ss_families
            summary_counts["ds_families"] += ds_families
            if min(count["A"], count["B"]) >= min_min_strand_reads and \
                max(count["A"], count["B"]) >= min_max_strand_reads:
                summary_counts["ds_duplexes"] += 1
        return summary_counts, kept_reads

    def run_downsamplings(self, reads, serial_sampling=True):
        """ When serial sampling is true, we use reads from the sampling of the
        next highest probability."""
        probs = np.sort(self.probabilities)[::-1]
        adj_probs = probs
        if serial_sampling:
            # If serial sampling, we need to adjust probability as kept_reads
            # is smaller after each sampling.
            adj_probs = adj_probs / np.insert(adj_probs, 0, 1)[:-1]

        for actual, prob in zip(probs, adj_probs):
            summary_counts, kept_reads = self.downsample(reads, prob)
            summary_counts = Counter(summary_counts)
            if serial_sampling:
                reads = kept_reads
            if self.per_target:
                if kept_reads:
                    if kept_reads[0].read1 and kept_reads[0].read2:
                        target, _ = kept_reads[0].get_overlap(self.interval_dict)
                    elif kept_reads[0].read1:
                        target, _ = get_overlap(kept_reads[0].read1, self.interval_dict)
                    elif kept_reads[0].read2:
                        target, _ = get_overlap(kept_reads[0].read2, self.interval_dict)
                    if not target:
                        print(kept_reads[0].read1, "\n", kept_reads[0].read2)
                    assert(target)
                    self.counts[target][actual] = (
                        self.counts[target][actual] + summary_counts
                    )
            else:
                self.counts[actual] = self.counts[actual] + summary_counts
