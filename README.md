# CODECsuite
CODECsuite is the software to process the CODEC data. It has 3 core functions: demultiplexing, adapter trimming, and single fragment mutation calling (SFC), which are written in c++14. 

## Installation
Tested on Red Hat 7 and Ubuntu 18.04

prerequisite for C++ based programs. For snakemake workflow check out [here](./snakemake)
1. git
2. gcc 5.2 or above with c++14 support
3. cmake 3.18.3 or above

First, recursive clone the repo and create a build directory which will holds the installion files and final executables.

`git clone --recursive git@github.com:broadinstitute/CODECsuite.git && cd CODECsuite && mkdir build`

Next, build the program with cmake.

`cd build && cmake .. && make`

After this, you should be able to see an executable named `codec` in the build folder you just created.

## Demultiplexing
CODECsuite is expected to work with raw lane-level fastq.gz. This can be obtained from illumina [bcl2fastq](https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html).
The first step is demultiplexing and it requires a sample sheet in csv format for each lane which looks the following.

| SampleName | IndexBarcode1 | IndexBarcode2 |
|------------|---------------|---------------|
| sample_A   |GAGCCTACTCAGTCAACG|GTGTCGAACACTTGACGG|
| sample_B   |CTTGAACGGACTGTCCAC|CACCGAGCGTTAGACTAC|

`codex demux -1 reads.r1.fastq.gz -2 reads.r2.fastq.gz -p sample_sheet.csv -o demux_outprefix `

Given the toy sample_sheet.csv and code this command will generate 
```
demux_outprefix.sample_A.1.fastq.gz, demux_outprefix.sample_A.2.fastq.gz
demux_outprefix.sample_B.1.fastq.gz, demux_outprefix.sample_B.2.fastq.gz
```

## Adapter trimming
After demultiplexing CODEC reads still contain in-situ sample barcode and adapter sequences. The next step is to trim 
these out since they could interfere alignment

`codec trim -1 demux_outprefix.sample_A.1.fastq.gz -2 demux_outprefix.sample_A.2.fastq.gz -o trim_outprefix -u 3 -U 3 -f 2 -t 2 -s sample_A` 

This tells the CODECsuite that first 3bp of a read is the UMI and to trim off the next two 2bp. 
The output files of the adapter trimming step looks like
```
trim_outprefix.sample_A.trim.bam
trim_outprefix.sample_A.trim.log
```
By default, single-end byproducts are also output to the `trim.bam`. To split the output use `-S/--split_bam_output`.

The bam file is standard uBam (unmapped bam) with additional tags
```
RX: UMI sequence from R1 and R2, concatenated by a hyphen
QX: UMI quality scores
bc: Index barcode sequence
s5: 5' adapter sequence (same as Index barcode)
q5: 5' adapter quality scores
s3: 3' adapter sequence (same as Index barcode of the mate)
q3: 3' adapter quality scores
sl: the rest of 3' adapter  sequence
ql: the rest of 3' adapter quality scores
```

After adapter trimming. The codec reads can be mapped by standard NGS tools such as BWA. For our end-to-end pipeline 
please see [snakemake](./snakemake).

## Single fragment caller (SFC) and mutation rate computation

After GATK best practice (alignment, markduplicate, indel realignment)  for example. Of note, BQSR should NOT be run for 
CODEC data since CODEC has a different quality score distribution. I do not recommend BQSR in general since the modern
Illumina sequencers' quality scores having been improved and BQSR almost doubles the bam size.

Now, we can run SFC to call mutations. SFC is designed to call somatic mutations. For the best results, we need to have
a bed file which contains the high confident regions (e.g., [this](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.1/GRCh37/LowComplexity/GRCh37_notinAllTandemRepeatsandHomopolymers_slop5.bed.gz)) in the reference genome and a germline bam for masking the germline 
variants. If there is no germline bam, it is recommend to have a germline vcf file. The population based vcf (e.g., [dbsnp](https://data.broadinstitute.org/snowman/hg19/Homo_sapiens_assembly19.dbsnp.vcf)) is almost
surely recommended since it can account for contamination and low sequencing depth of germline bam. There is a set of 
fragment level and base level filters to improve the precision of the mutation calls at the cost of data loss and potential
loss of real mutations. Depends on the applications, we have presets of parameters: `-p/--preset`
```
stringent: setting where high precision calling is needed. Situations like calling background mutation rate in white 
blood cells

lenient: setting where certain sensitivity is prefer. situations like calling cancer mutations in tumor biopsy or high 
tumor fraction ctDNA samples. 

null: as little as filtering possible. For advanced users who want to filtering on there own ends.  
```

It is highly recommended that a user try and play with different parameters and figure out the ones that are best for
themselves. Trimming the end `-d 12` and use Q30 cutoff `-q 30` is always recommended. The others are ad hoc. However, 
the most effective parameter is probably `-Q/--min_passQ_frac`: the fraction of Q30+Q30 bases in the overlap region. The fraction essentially
measures the cluster quality, which is important for single fragment calling. Some examples of running SFC
```
codec call -b input.mark_duplicated.bam -L highconfidentregions.bed  -r hg19.fa -n germline.bam -p lenient -o output

```

The output of the SFC are
```
output.mutation_metrics.txt: includes SNV_rate, INDEL_rate and etc. 
output.variatns_called.txt: mutations from single fragments
output.context_count.txt: trinucleotide context and dinucleotide context counts
```

## Other notes
1. For CODEC-MSI please refer to [msi](./msi). And by default CMAKE will not build CODEC-MSI. Please uncomment the last two
lines if you indeed want to build CODEC-MSI

2. The Snakemake is hard-coded to de-multiplex 4 lanes simultaneously (e.g. for NovaSeq 6000). If you need to de-multiplex
less #lanes (e.g. for NovaSeq SP), comment out entire rules for DemuxL3 and DemuxL4. If you have more than 4 lanes (e.g. HiSeq X)
either do 4 lane at a times or add more rules yourself. 

3. The Snakemake pipeline setup file `qsub_wrapper.py` is specific to [UGE](https://en.wikipedia.org/wiki/Univa_Grid_Engine).
You may need to change settings for your computing environment. 