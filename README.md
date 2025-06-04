UPDATES
* (05/19/25) Version 1.1.5 introduces new scripts for converting CODEC single fragment duplex variants output to MAF or VCF format.
  [See details here](#reformatting-codec-sfc-variant-output-into-maf-or-vcf-format)

*  (01/08/25) Version 1.1.4 introduces a new script `codec filter` designed to filter consensus BAM files. It retains only the reads and bases relevant for variant calling. Fragments (read-pairs) that do not pass fragment-level filtering are excluded from the output BAM. Bases that fail the filters are assigned a minimum base quality score (Q2), ensuring they are ignored by most coverage analysis and variant calling tools.
   It can be run as the following:
   
   ``` codec filter -b mol_consensus.sorbybyname.bam -o duplex_only.bam -r reference.fa -q 30 -m 60 -Q 0.7 -B 0.5 -N 0.05 ...```
> [!NOTE]
> The input BAM file must be sorted by read name, and the output BAM will also be query-name sorted. For consistent filtering results, it is recommended to use the same parameters as those in codec call.


# CODECsuite
The CODEC analysis pipeline, CODECsuite, comprises five key steps: demultiplexing, adapter trimming, alignment, duplicate collapsing, and single-fragment mutation calling. Duplicate collapsing and alignments are performed using the third-party tools Fgbio and BWA, respectively. After removing byproducts and applying fragment-level filtering, mutations were identified exclusively from duplexes in the overlapped regions, where bases from each read align and match. Bases within these regions underwent stringent filtering based on criteria such as base quality, proximity to fragment ends, overlap with germline mutations, and other factors. Notably, a single read pair is sufficient to form a duplex, as each read represents one strand. Refer to the [paper](https://www.nature.com/articles/s41588-023-01376-0) for more details. 

## Installation
Tested on Red Hat 7 and Ubuntu 18.04

prerequisite for C++ based programs. For snakemake workflow check out [here](./snakemake)
1. git
2. tested with gcc 5.2 and 7.3 with c++14 support
3. cmake 3.18.3 or above

First, recursive clone the repo and create a build directory which will holds the installion files and final executables.

`git clone --recursive git@github.com:broadinstitute/CODECsuite.git && cd CODECsuite && mkdir build`

Next, build the program with cmake.

`cd build && cmake .. && make`

After this, you should be able to see an executable named `codec` in the build folder you just created.

## Demultiplexing
CODECsuite is expected to work with raw lane-level fastq.gz. This can be obtained from illumina [bcl2fastq](https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html).
The first step is demultiplexing and it requires a sample sheet in csv format for each lane which looks the following.
Currently, we have used 12 barcodes. For good cluster generation, we recommend to have at least 4 sample barcodes per
sequencing lane. 

| SampleName | IndexBarcode1 | IndexBarcode2 |
|------------|---------------|---------------|
|Sample01|CTTGAACGGACTGTCCAC|CACCGAGCGTTAGACTAC|
|Sample02|GAGCCTACTCAGTCAACG|GTGTCGAACACTTGACGG|
|Sample03|AGCTTGTAAGGCAGGTTA|ACTGATCTTCAGCTGACT|
|Sample04|TCAAGCGTCTTACATGGT|TGAATCTGAGGCACTGTA|
|Sample05|CTGGTCCAAGAACGTCTG|CTCTGAACGATCGAGCTC|
|Sample06|GATCCAGTTCTGTCGAGC|GAGGTGCATGCACCTTAG|
|Sample07|ACCTATAGGTGCAACGAA|ACTAACTTCCATTGCACT|
|Sample08|TGAAGGTCCACTGTATCT|TGACCTGGATGGATAGGA|
|Sample09|CACTGCTTCGAGACGAAG|CTCCAGTTACTGAGACGG|
|Sample10|GTGATACCTCGATGCTCC|GAGGTCCAGTCTCTGTCC|
|Sample11|ACTCAGAGAACTCATGGA|ACTACAGGTGGATCCAAT|
|Sample12|TGAGCTGAGTTCGTACTT|TGATGTACCAACGATGTA|

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

> [!NOTE]
> We recommend using SMaHT duplex reference genome which is basically HG38 without decoy sequences. See reasons here: https://smaht-dac.github.io/pipelines-docs/DOCS/REFERENCE_FILES/Genome_Builds/1_Build_GRCh38.html

## Single fragment caller (SFC) and mutation rate computation

After GATK best practice (alignment, markduplicate, indel realignment)  for example. Of note, BQSR should NOT be run for 
CODEC data since CODEC has a different quality score distribution. I do not recommend BQSR in general since the modern
Illumina sequencers' quality scores having been improved and BQSR almost doubles the bam size.

Now, we can run SFC to call mutations. SFC is designed to call somatic mutations. For the best results, we need to have
a bed file which contains the high confident regions (e.g., [this](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.3/GRCh38@all/LowComplexity/GRCh38_notinAllTandemRepeatsandHomopolymers_slop5.bed.gz)) in the reference genome and a germline bam for masking the germline 
variants. If there is no germline bam, it is recommend to have a germline vcf file. The population based vcf (e.g., [gnomad vcf](https://storage.googleapis.com/gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz)) is almost
surely recommended since it can account for contamination and low sequencing depth of germline bam. However, to avoid over-filtering true somatic mutations, a minimum allele frequency threshold is recommend for the population vcf (e.g., 0.01%)

There is a set of 
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
codec call -b input.mark_duplicated.bam -L highconfidentregions.bed  -r hg38.fa -n germline.bam -p lenient -o output

```

The output of the SFC are
```
output.mutation_metrics.txt: includes SNV_rate, INDEL_rate and etc. 
output.variatns_called.txt: mutations from single fragments
output.context_count.txt: trinucleotide context and dinucleotide context counts
```

> [!NOTE]
> All CODEC related resources can be found at https://console.cloud.google.com/storage/browser/codec_cloud_resources. Including the population based vcf: https://storage.googleapis.com/codec_cloud_resources/alfa_all.freq.breakmulti.hg38.af0001.vcf.gz

## Reformatting CODEC SFC variant output into MAF or VCF format 
Scripts `codec2maf` and `maf2vcf.py` can be found in folder `snakemake/script/`.  `maf2vcf.py` depends on `maf2vcf.pl` script from the perl package [mskcc/vcf2maf](https://github.com/mskcc/vcf2maf/tree/main)
1. CODEC txt file To MAF: `codec2maf -i output.variatns_called.txt -o output.variatns_called.maf`
2. MAF To VCF: `maf2vcf.py output.variatns_called.maf -r hg38.fa -o outdir -p /usr/bin/maf2vcf.pl`


## Other notes
1. For CODEC-MSI please refer to [msi](./msi). And by default CMAKE will not build CODEC-MSI. Please uncomment the last two
lines if you indeed want to build CODEC-MSI

2. The Snakemake is hard-coded to de-multiplex 4 lanes simultaneously (e.g. for NovaSeq 6000). If you need to de-multiplex
less #lanes (e.g. for NovaSeq SP), comment out entire rules for DemuxL3 and DemuxL4. If you have more than 4 lanes (e.g. HiSeq X)
either do 4 lane at a times or add more rules yourself. 

3. The Snakemake pipeline setup file `qsub_wrapper.py` is specific to [UGE](https://en.wikipedia.org/wiki/Univa_Grid_Engine).
You may need to change settings for your computing environment. 
