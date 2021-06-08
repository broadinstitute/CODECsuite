# CODECsuite
CODECsuite is the software to process the CODEC data. It has 4 core functions: demultiplexing, adapter trimming, consensus calling and computing accuracy, which are written in c++14. We also use [snakemake](./snakemake) to manage the end-to-end workflow. 

## Installation
Tested on Red Hat 7 and Ubuntu 18.04

prerequisite for C++ based programs. For snakemake workflow check out [here](./snakemake)
1. gcc 5.2 or above with c++14 support
2. cmake 3.18 or above
3. elfutils:
https://sourceware.org/elfutils/. After download and install elfutils, set the environmental variable `ELFUTILS_ROOT` to the directory. 

First, create a build directory which will holds the installion files and final executables.

`mkdir build ` 

Next, build the program with cmake. 

`cd build && cmake .. && make`

After this, you should be able to see an executable named `codec` in the build folder you just created. 

```
Usage: codec demux [options]

-p/--library_param,                    Sample, barcode mapping in CSV format. Header must be "SampleName,IndexBarcode1,IndexBarcode2", required
-1/--q1,                               Input read1, required
-2/--q2,                               Input read2, required
-b/--index_begin,                      The read position where the index begins (Default: 3)
-l/--index_len,                        Index length (Default: 18)
-e/--max_ed,                           Maximum edit distance allowed as a match (Default: 3)
-o/--outprefix,                        Output path, e.g., /tmp/test
-r/--ref,                              Reference genome fasta file, for judging index hopping
-i/--include_non_pf,                   Include non-pass filter reads
-v/--verbose,                          Print verbose information
-c/--count_pf,                         Just count number of pass filter pairs. Do not do anything else
```

```
Usage: codec trim [options]

-1/--R1,                               R1.fastq. [required]
-2/--R2,                               R1.fastq. [required]
-o/--prefix,                           Output prefix. [/tmp/cds.]
-n/--nbase_skip,                       Do not trim in the first n base of the read. [n=0]
-i/--rgid,                             Read group id when output bam [A]
-s/--rgsm,                             Read group sample when output bam [A]
-T/--min_bq_from_back,                 After adpters being removed, trim bases with bq less than this value from the back [3].
-u/--r1_umi_len,                       num of umi bases to be trimmed from 5'end of read1 [0].
-U/--r2_umi_len,                       num of umi bases to be trimmed from 5'end of read2 [0].
-F/--trim_full_linker_len,             Minimum trimming length if linker is found in the middle of a read [30].
-P/--trim_prefix_linker_len,           Minimum trimming length if linker is found at the 3'end of a read [3]
-t/--post_trim3,                       Num of bases to be trimmed at the 3'end of a read after adapter removed [0]
-f/--post_trim5,                       Num of bases to be trimmed at the 5'end of a read after adapter removed [0]
-A/--match                             Score for a sequence match [1]
-B/--mismatch                          Penalty for a mismatch [4]
-G/--gap                               Penalty for open or extend a gap [5]
-d/--debug,                            1: detail pairwise alignment plot, 2: Qscore plot, default no debug[0]
```

```
Usage: codec consensus [options]
-b/--bam,                              Input bam [required]
-o/--outbam,                           Output unmapped bamfile [required].
-m/--mapq,                             Min mapping quality [10].
-q/--baseq,                            Min base quality for consensus. Otherwise, masked with N [20].
-l/--load_supplementary,               Include supplementary alignment [false].
-t/--trim_overhang,                    When perform paired-end consensus, if true then only do consensus of the overlapped region [false].
-C/--clip3,                            trim the 3'end soft clipping [false].
-p/--pair_min_overlap,                 When using selector, the minimum overlap between the two ends of the pair [1]. -1 means two segs complete overlap excluding the soft clip part [-1]. 0 for allowing non-overlapping pair
-d/--dirtmp,                           Temporary dir for sorted bam [/tmp]
-T/--thread,                           Number of threads for sort [1]
-i/--allow_nonoverlapping_pair,        Allow output of non-overlaping pairs, usually caused by intermolecular ligation. This will simply print the original reads.  [false]
```

```
Usage: codec accuracy [options]
General Options:
-v/--verbose,                          [default 0]
-b/--bam,                              input bam
-L/--bed,                              targeted region
-m/--mapq,                             min mapping quality [10].
-S/--load_supplementary,               include supplementary alignment [false].
-S/--load_secondary,                   include secondary alignment [false].
-u/--load_unpair,                      include unpaired alignment [false].
-V/--vcfs,                             comma separated VCF file(s) for germline variants or whitelist variants[null].
-M/--maf,                              MAF file for somatic variants [null].
-s/--sample,                           sample from the VCF file [null].
-r/--reference,                        reference sequence in fasta format [null].
-a/--accuracy_stat,                    output reporting accuracy for each alignment [accuracy_stat.txt].
-e/--error_prof_out,                   Error profile output in plain txt format [error_prof_out.txt].
-k/--known_var_out,                    Output for known var. [known_var_out.txt].
--detail_qscore_prof,                  Output finer scale qscore cutoffs, error rates profile. The default is only q0, q30 [false].
--read_level_stat,                     Output read level stat.

Filtering Options:
-q/--bqual_min,                        Skip bases with baseQ smaller than this when calculating error rate [0].
-n/--max_edit_filter,                  Skip a read if its NM tag is larger than this value [INT_MAX].
-x/--max_nonNedit_filter,              Skip a read if the number of non-N bases edits is larger than this value [INT_MAX].
-d/--fragend_dist_filter,              Consider a variant if its distance to the fragment end is at least this value [0].
-p/--pair_min_overlap,                 When using selector, the minimum overlap between the two ends of the pair. -1 for complete overlap, 0 no overlap required [0].
-O/--overlap_only,                     Count only overlapped region of a read pair. Default [false].
```
