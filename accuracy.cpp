//
// Created by Ruolin Liu on 3/3/20.
//
//
// Created by Ruolin Liu on 2/18/20.
//

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <chrono>
#include <ctime>

#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamWriter.h"
#ifdef BGZF_MAX_BLOCK_SIZE
#pragma push_macro("BGZF_MAX_BLOCK_SIZE")
#undef BGZF_MAX_BLOCK_SIZE
#define BGZF_MAX_BLOCK_SIZE_BAK
#endif

#ifdef BGZF_BLOCK_SIZE
#pragma push_macro("BGZF_BLOCK_SIZE")
#undef BGZF_BLOCK_SIZE
#define BGZF_BLOCK_SIZE_BAK
#endif

#include "InsertSeqFactory.h"
#include "ReadVCF.h"
#include "Variant.h"
#include "Alignment.h"
#include "StringUtils.h"
#include "MAF.h"
#include "TargetLayout.h"
#include "MutCounter.h"
#include "Algo.h"
#include "pileup.h"

//#define OPT_QSCORE_PROF   261
#define OPT_READ_LEVEL_STAT   262
#define OPT_CYCLE_LEVEL_STAT   263
//#define OPT_MIN_QPASS_RATE_T2G   264
//#define OPT_ACCU_BURDEN 265
#define OPT_MIN_QPASS_RATE_TT   266
#define OPT_MAX_GERM_INDEL_DIST  267
#define OPT_ALLOW_INDEL_NEAR_SNV  268
#define OPT_MIN_GERM_MAPQ  269
#define OPT_MIN_NEAREST_SNV  270
#define OPT_MIN_NEAREST_INDEL  271

using std::string;
using std::vector;
int MYINT_MAX = std::numeric_limits<int>::max();

struct AccuOptions {
 //input output
  string bam;
  string germline_bam;
//  string raw_bam;
  string reference;
  string bed_file;
  vector<string> vcfs;
  string maf_file;

  string mutation_metrics;
  string variants_out;
  string context_count;
  string known_var_out;
  string sample;
  string read_level_stat;
  string cycle_level_stat;
  string preset;


// filtering options
  bool load_supplementary = false;
  bool load_secondary = false;
  bool load_unpair = false;
  bool load_duplicate = false;
  bool load_proper_pair_only = false;

  int count_read = 0;
  int min_indel_len = 0;
  bool standard_ngs_filter = false;

  int mapq = 20;
  int bqual_min = 20;
  int fragend_dist_filter = 0;
  float min_passQ_frac_TT = 0;
  float min_passQ_frac = 0;
  bool filter_5endclip = false;
  int max_N_filter = MYINT_MAX;
  int max_snv_filter = MYINT_MAX;
  int min_fraglen = 30;
  int max_fraglen = MYINT_MAX;
  int verbose = 0;
  int clustered_mut_cutoff = MYINT_MAX;
  float max_frac_prim_AS = std::numeric_limits<float>::max();
  int germline_minalt = 1;
  int germline_minmapq = 20;
  int germline_mindepth = 5;
  int germline_var_maxdist = 5;
  float germline_cutoff_vaf = 1.0;
  bool allow_indel_near_snv = false;
  float max_pair_mismatch_frac = 1.0;
  string MID_tag = "";
  int min_indel_dist_from_readend = 3;
  int8_t IndelAnchorBaseQ = 53; // 33 as shift
  int MIN_NEAREST_SNV = 3;
  int MIN_NEAREST_INDEL = 10;

  //hidden parameter
  int indel_anchor_size= 1;
  int germline_minbq = 10;
  bool all_mutant_frags = false;
  bool detail_qscore_prof = false;
  int pair_min_overlap = 0;
  float family_agree_rate = 0.95;

//obsolete
//  double max_mismatch_frac = 1.0;
//  float min_passQ_frac_T2G = 0;
//  string trim_bam = "trimmed.bam";
};


static struct option  accuracy_long_options[] = {
    {"bam",                      required_argument,      0,        'b'},
    {"normal_bam",               required_argument,      0,        'n'},
//    {"raw_bam",                  required_argument,      0,        'l'},
    {"bed",                      required_argument,      0,        'L'},
    {"preset",                   required_argument,      0,        'p'},
    {"mutation_metrics",         required_argument,      0,        'a'},
    {"load_unpair",              no_argument,            0,        'u'},
    {"load_supplementary",       no_argument,            0,        'S'},
    {"load_secondary",           no_argument,            0,        '2'},
    {"load_duplicate",           no_argument,            0,        'D'},
    {"load_only_proper_pair",    no_argument,            0,        'P'},
    {"MID_tag",                  required_argument,      0,        'U'},
    {"mapq",                     required_argument ,     0,        'm'},
    {"vcfs",                     required_argument ,     0,        'V'},
    {"maf",                      required_argument ,     0,        'M'},
    {"output",                   required_argument,      0,        'o'},
    {"reference",                required_argument,      0,        'r'},
    {"variants_out",             required_argument,      0,        'e'},
    {"bqual_min",                required_argument,      0,        'q'},
    {"min_fraglen",              required_argument,      0,        'g'},
    {"max_fraglen",              required_argument,      0,        'G'},
    {"min_germdepth",            required_argument,      0,        'Y'},
    {"max_germindel_dist",       required_argument,      0,        OPT_MAX_GERM_INDEL_DIST},
    {"min_dist_to_nearest_SNV",  required_argument,      0,        OPT_MIN_NEAREST_SNV},
    {"min_dist_to_nearest_INDEL",required_argument,      0,        OPT_MIN_NEAREST_INDEL},
    {"max_frac_prim_AS",         required_argument,      0,         'B'},
//    {"max_mismatch_frac",        required_argument,      0,        'F'},
    {"min_germline_alt",         required_argument,      0,        'W'},
    {"min_germline_mapq",        required_argument,      0,        OPT_MIN_GERM_MAPQ},
    {"min_indel_anchor_baseq",   required_argument,      0,        'f'},
    {"min_indel_dist_readend",   required_argument,      0,        'E'},
    {"germline_cutoff_vaf",      required_argument,      0,        'i'},
    {"disable_5endclip_filtering",    no_argument,            0,        '5'},
    {"allow_indel_near_snv",     no_argument,            0,        OPT_ALLOW_INDEL_NEAR_SNV},
    {"min_indel_len",            required_argument,            0,        'I'},
    {"min_passQ_frac",           required_argument,      0,        'Q'},
    {"max_pair_mismatch_frac",               required_argument,      0,        'N'},
//    {"min_passQ_frac_T2G",       required_argument,      0,        OPT_MIN_QPASS_RATE_T2G},
    {"min_passQ_frac_TT",        required_argument,      0,        OPT_MIN_QPASS_RATE_TT},
    {"known_var_out",            required_argument,      0,        'k'},
    {"context_count",            required_argument,      0,        'C'},
//    {"pair_min_overlap",         required_argument,      0,        'p'},
    {"count_read",               required_argument,            0,        'R'},
    {"fragend_dist_filter",      required_argument,      0,        'd'},
    {"max_N_filter",             required_argument,      0,        'y'},
    {"max_snv_filter",           required_argument,      0,        'x'},
    {"standard_ngs_filter",      no_argument,            0,        's'},
    {"clustered_mut_cutoff",     required_argument,      0,        'c'},
    {"verbose",                  required_argument,      0,        'v'},
    {"all_mutant_frags",         no_argument,             0,        'A'},
//    {"detail_qscore_prof",       no_argument,            0,        OPT_QSCORE_PROF},
    {"read_level_stat",          required_argument,      0,        OPT_READ_LEVEL_STAT},
    {"cycle_level_stat",         required_argument,      0,        OPT_CYCLE_LEVEL_STAT},
//   {"accu_burden",         no_argument,      0,        OPT_ACCU_BURDEN},
    {0,0,0,0}
};
const char* accuracy_short_options = "b:a:m:v:S2uo:r:e:q:k:R:M:p:d:n:x:V:L:DC:N:5Q:g:G:I:c:B:Y:W:U:f:i:sPE:";

void accuracy_print_help()
{
  std::cerr<< "---------------------------------------------------\n";
  std::cerr<< "Usage: codec accuracy [options]\n";
  std::cerr<< "Suggested Options:\n";
  std::cerr<< "-p/--preset,                           options preset. choice of [stringent, lenient, null]\n";
  std::cerr<< "-b/--bam,                              input bam\n";
  std::cerr<< "-n/--normal_bam,                       input normal bam [optional]\n";
//  std::cerr<< "-l/--raw_bam,                          input raw bam [optional]\n";
  std::cerr<< "-L/--bed,                              targeted region\n";
  std::cerr<< "-o/--output,                           output prefix. -a, -e, -C overwrite this options [null].\n";
  std::cerr<< "-r/--reference,                        reference sequence in fasta format [null].\n";
  std::cerr<< "-V/--vcfs,                             comma separated VCF file(s) of blacklist variants such as germline variants or populations variants (e.g. dbSNP)[null].\n";

  std::cerr<< "\nFor standard NGS, use -u -R 1 -s:\n";
  std::cerr<< "-u/--load_unpair,                      include unpaired alignment [false].\n";
  std::cerr<< "-R/--count_read,                       0: collapse R1,R2; only overlapped region. 1: collapse R1,R2; include overhangs. 2: count independently for R1 and R2, including overhang [0].\n";
  std::cerr<< "-s/--standard_ngs_filter,              Filter for standard NGS (where you expect very little overlap between R1 and R2) [false].\n";

  std::cerr<< "\nOther Options:\n";
  std::cerr<< "-v/--verbose,                          [default 0]\n";
  std::cerr<< "-U/--MID_tag,                          molecular identifier tag name [MI] (MI is used in fgbio). This is only used when -D is on. During this mode, reads will\n";
  std::cerr<< "-M/--maf,                              MAF file for somatic variants [null].\n";
  std::cerr<< "-k/--known_var_out,                    Output variants which match the blacklist vcfs (from -V). [default null].\n";

  std::cerr<< "\nFiltering Options:\n";
  std::cerr<< "-S/--load_supplementary,               include supplementary alignment [false].\n";
  std::cerr<< "-P/--load_only_proper_pair,            improper pair is included by default [false].\n";
  std::cerr<< "-2/--load_secondary,                   include secondary alignment [false].\n";
  std::cerr<< "-D/--load_duplicate,                   include alignment marked as duplicates [false].\n";
  std::cerr<< "-q/--bqual_min,                        Skip bases if min(q1, q2) < this when calculating error rate. q1, q2 are baseQ from R1 and R2 respectively [20].\n";
  std::cerr<< "-m/--mapq,                             min mapping quality [20].\n";
  std::cerr<< "-Q/--min_passQ_frac,                   Filter out a read if the fraction of bases passing quality threshold (together with -q) is less than this number [0].\n";
  std::cerr<< "-x/--max_snv_filter,                   Skip a read if the number of mismatch bases is larger than this value [INT_MAX].\n";
  std::cerr<< "-y/--max_N_filter,                     Skip a read if its num of N bases is larger than this value [INT_MAX].\n";
  std::cerr<< "-d/--fragend_dist_filter,              Consider a variant if its distance to the fragment end is at least this value [0].\n";
//  std::cerr<< "-F/--max_mismatch_frac,                Filter out a read if its mismatch fraction is larger than this value  [1.0].\n";
  std::cerr<< "-N/--max_pair_mismatch_frac,           Filter out a read-pair if its R1 and R2 has mismatch fraction larger than this value  [1.0].\n";
  std::cerr<< "-W/--min_germline_alt,                 Minimum number of germline alt reads to be consider as a germline site  [1].\n";
  std::cerr<< "--min_germline_mapq,                   Minimum mapq of germline reads [20].\n";
  std::cerr<< "--max_germindel_dist,                  Filter out a INDEL if its distance to a germline INDEL is less than this number [5].\n";
  std::cerr<< "-i/--germline_cutoff_vaf,              Consider a variant is germline is the VAF is larger than this number [1.0].\n";
//  std::cerr<< "--min_passQ_frac_T2G,                  Filter out T>G SNV if the fraction of of pass baseq smaller than this value  [0].\n";
  std::cerr<< "--min_passQ_frac_TT,                   Filter out T>G SNV in the context of TT if the fraction of of pass baseq smaller than this value  [0].\n";
  std::cerr<< "-5/--filter_5endclip,                  Filtering out reads with 5'end soft clipping [False].\n";
  std::cerr<< "-I/--min_indel_len,                    any positive values indicate indels calls only [0].\n";
  std::cerr<< "-f/--min_indel_anchor_baseq,           minimum baseq for the anchoring bases of a indel [20].\n";
  std::cerr<< "-E/--min_indel_dist_readend            minimum distant of a INDEL to the end of a read [3].\n";
  std::cerr<< "--allow_indel_near_snv,                allow SNV to pass filter if a indel is found in the same read [false].\n";
  std::cerr<< "-g/--min_fraglen,                      Filter out a read if its fragment length is less than this value [30].\n";
  std::cerr<< "-B/--max_frac_prim_AS,                 Filter out a read if the AS of secondary alignment is within this fraction of the primary alignment [0.75].\n";
  std::cerr<< "-Y/--min_germdepth,                    Minimum depth in germline bam [5].\n";
  std::cerr<< "-G/--max_fraglen,                      Filter out a read if its fragment length is larger than this value [INT_MAX].\n";
  //std::cerr<< "-p/--pair_min_overlap,                 When using selector, the minimum overlap between the two ends of the pair. -1 for complete overlap, 0 no overlap required [0].\n";
  std::cerr<< "-c/--clustered_mut_cutoff,             Filter out a read if at least this number of mutations occur in a window of 30 bp near the read end (<15 bp). [INT_MAX].\n";
  std::cerr<< "--min_dist_to_nearest_SNV,             The indel of interest is at least x nt from another SNV. [x = 4].\n";
  std::cerr<< "--min_dist_to_nearest_INDEL,           The indel of interest is at least x nt from another INDEL. [x = 10].\n";

  std::cerr<< "\n Obsolete Options:\n";
  std::cerr<< "-a/--mutation_metrics,                 mutation metrics file [.mutation_metrics.txt].\n";
  std::cerr<< "-e/--variants_out,                     mutations in plain txt format [.variants_called.txt].\n";
  std::cerr<< "-C/--context_count,                    Output for trinucleotide and dinucleotide context context in the reference. [.context_count.txt].\n";
  //std::cerr<< "--qscore_prof,                         Output qscore prof. First column is qscore cutoff; second column is number of bases in denominator\n";
  std::cerr<< "--detail_qscore_prof,                  Output finer scale qscore cutoffs, error rates profile. The default is only q0, q30 [false]. \n";
  std::cerr<< "--read_level_stat,                     Output read level error metrics.\n";
  std::cerr<< "--cycle_level_stat,                    Output cycle level error metrics.\n";
  //std::cerr<< "-A/--all_mutant_frags,                 Output all mutant fragments even if not pass failters. Currently only works for known vars [false].\n";
}

int accuracy_parse_options(int argc, char* argv[], AccuOptions& opt) {
  int option_index;
  int next_option = 0;
  do {
    next_option = getopt_long(argc, argv, accuracy_short_options, accuracy_long_options, &option_index);
    switch (next_option) {
      case -1:break;
      case 'v':
        opt.verbose = atoi(optarg);
        break;
      case 'b':
        opt.bam = optarg;
        break;
      case 'n':
        opt.germline_bam = optarg;
        break;
//      case 'l':
//        opt.raw_bam = optarg;
//        break;
      case 'L':
        opt.bed_file = optarg;
        break;
      case 'a':
        opt.mutation_metrics = optarg;
        break;
      case 'm':
        opt.mapq = atoi(optarg);
        break;
      case 'q':
        opt.bqual_min = atoi(optarg);
        break;
      case 'Q':
        opt.min_passQ_frac = atof(optarg);
        break;
      case 'y':
        opt.max_N_filter = atoi(optarg);
        break;
      case 'x':
        opt.max_snv_filter = atoi(optarg);
        break;
      case 'B':
        opt.max_frac_prim_AS = atof(optarg);
        break;
      case 'Y':
        opt.germline_mindepth = atoi(optarg);
        break;
      case 'f':
        opt.IndelAnchorBaseQ = atoi(optarg) + 33;
        break;
      case 'E':
        opt.min_indel_dist_from_readend = atoi(optarg);
        break;
      case OPT_MAX_GERM_INDEL_DIST:
        opt.germline_var_maxdist = atoi(optarg);
        break;
      case OPT_MIN_NEAREST_SNV:
        opt.MIN_NEAREST_SNV = atoi(optarg);
        break;
      case OPT_MIN_NEAREST_INDEL:
        opt.MIN_NEAREST_INDEL = atoi(optarg);
        break;
      case 'c':
        opt.clustered_mut_cutoff = atoi(optarg);
        break;
      case 'd':
        opt.fragend_dist_filter = atoi(optarg);
        break;
//      case 'F':
//        opt.max_mismatch_frac = atof(optarg);
//        break;
      case 'N':
        opt.max_pair_mismatch_frac = atof(optarg);
        break;
//      case OPT_MIN_QPASS_RATE_T2G:
//        opt.min_passQ_frac_T2G = atof(optarg);
//        break;
      case OPT_MIN_QPASS_RATE_TT:
        opt.min_passQ_frac_TT = atof(optarg);
        break;
      case 'g':
        opt.min_fraglen = atoi(optarg);
        break;
      case 'G':
        opt.max_fraglen = atoi(optarg);
        break;
      case 'S':
        opt.load_supplementary = true;
        break;
      case 'U':
        opt.MID_tag = optarg;
        break;
      case '2':
        opt.load_secondary = true;
        break;
      case 'u':
        opt.load_unpair = true;
        break;
      case 'D':
        opt.load_duplicate = true;
        break;
      case 'P':
        opt.load_proper_pair_only = true;
        break;
      case '5':
        opt.filter_5endclip = true;
        break;
      case OPT_ALLOW_INDEL_NEAR_SNV:
        opt.allow_indel_near_snv = true;
        break;
      case 's':
        opt.standard_ngs_filter = true;
        break;
      case 'I':
        opt.min_indel_len = atoi(optarg);
        break;
      case 'R':
        opt.count_read = atoi(optarg);
        break;
      case 'o':
        opt.sample = optarg;
        break;
      case 'V':
        opt.vcfs = cpputil::split(optarg, ",");
        break;
      case 'M':
        opt.maf_file = optarg;
        break;
      case 'r':
        opt.reference = optarg;
        break;
      case 'e':
        opt.variants_out = optarg;
        break;
      case 'k':
        opt.known_var_out = optarg;
        break;
      case 'C':
        opt.context_count = optarg;
        break;
//      case OPT_QSCORE_PROF:
//        opt.detail_qscore_prof = true;
//        break;
      case OPT_READ_LEVEL_STAT:
        opt.read_level_stat = optarg;
        break;
      case OPT_CYCLE_LEVEL_STAT:
        opt.cycle_level_stat = optarg;
        break;
      case 'p':
        opt.preset = optarg;
        break;
      case 'W':
        opt.germline_minalt = atoi(optarg);
        break;
      case 'i':
        opt.germline_cutoff_vaf = atof(optarg);
        break;
      default:accuracy_print_help();
        return 1;
    }
  } while (next_option != -1);

  return 0;
}


void CycleBaseCount(const cpputil::Segments& seg,
                    const SeqLib::GenomicRegion* const gr,
                    cpputil::ErrorStat& es) {
  if (seg.size() == 1) {
    std::pair<int,int> range;
    if (gr) {
      range = cpputil::GetBamOverlapQStartAndQStop(seg.front(), *gr);
    } else {
      range.first = seg.front().AlignmentPosition();
      range.second = seg.front().AlignmentEndPosition();
    }
    if (seg.front().FirstFlag()) {
      cpputil::AddMatchedBasesToCycleCount(seg.front(), es.R1_q0_cov, es.R1_q30_cov, range.first, range.second);
    }
    else {
      cpputil::AddMatchedBasesToCycleCount(seg.front(), es.R2_q0_cov, es.R2_q30_cov, range.first, range.second);
    }
  } else {
    std::pair<int, int> overlap_front, overlap_back;
    if (gr) {
      overlap_front = cpputil::GetBamOverlapQStartAndQStop(seg.front(), *gr);
      overlap_back = cpputil::GetBamOverlapQStartAndQStop(seg.back(), *gr);
    } else {
      overlap_front.first = seg.front().AlignmentPosition();
      overlap_front.second = seg.front().AlignmentEndPosition();
      overlap_back.first = seg.back().AlignmentPosition();
      overlap_back.second = seg.back().AlignmentEndPosition();
    }
    if (seg.front().FirstFlag()) {
      cpputil::AddMatchedBasesToCycleCount(seg.front(), es.R1_q0_cov, es.R1_q30_cov, overlap_front.first, overlap_front.second);
      cpputil::AddMatchedBasesToCycleCount(seg.back(), es.R2_q0_cov, es.R2_q30_cov, overlap_back.first, overlap_back.second);
    } else {
      cpputil::AddMatchedBasesToCycleCount(seg.front(), es.R2_q0_cov, es.R2_q30_cov, overlap_front.first, overlap_front.second);
      cpputil::AddMatchedBasesToCycleCount(seg.back(), es.R1_q0_cov, es.R1_q30_cov, overlap_back.first, overlap_back.second);
    }
  }
}


int n_false_mut(vector<bool> v) {
  return count(v.begin(), v.end(), false);
}
int n_true_mut(vector<bool> v) {
  return count(v.begin(), v.end(), true);
}

void ErrorRateDriver(vector<cpputil::Segments>& frag,
                    const int pair_nmismatch,
                    const float nqpass,
                    int olen,
                    const SeqLib::BamHeader& bamheader,
                    const SeqLib::RefGenome& ref,
                    const SeqLib::BWAWrapper& bwa,
                    const SeqLib::GenomicRegion* const gr,
                    const std::set<int> blacklist,
                    const AccuOptions& opt,
                    const cpputil::BCFReader& bcf_reader,
                    const cpputil::BCFReader& bcf_reader2,
                    const cpputil::MAFReader& mafr,
                    cpputil::PileHandler& pileup,
                    cpputil::PileHandler& tumorpileup,
                    std::ofstream& ferr,
                    std::ofstream& known,
                    std::ofstream& readlevel,
                    cpputil::ErrorStat& errorstat) {

  for (auto& seg : frag) { // a fragment may have multiple segs due to supplementary alignments
    std::vector<std::string> orig_qualities(2);
    std::vector<std::string> orig_seqs(2);
    //EOF filter
    if (seg.size() == 2) {
      orig_seqs = {seg[0].Sequence(), seg[1].Sequence()};
      orig_qualities = {seg[0].Qualities(), seg[1].Qualities()};
      ++errorstat.n_pass_filter_pairs;
      cpputil::TrimPairFromFragEnd(seg.front(), seg.back(), opt.fragend_dist_filter);
    } else if(seg.size() == 1) {
      if (seg[0].FirstFlag()) {
        orig_seqs[0] = seg[0].Sequence();
        orig_qualities[0] = seg[0].Qualities();
      } else {
        orig_seqs[1] = seg[0].Sequence();
        orig_qualities[1] = seg[0].Qualities();
      }
      ++errorstat.n_pass_filter_singles;
      cpputil::TrimSingleFromFragEnd(seg.front(), opt.fragend_dist_filter);
    }

    //for readlevel output
    //int r1_q0_den = 0, r2_q0_den = 0;
    //int r1_q0_nerror = 0, r2_q0_nerror = 0;
    int r1_nerror = 0, r2_nerror = 0;
    string chrname = seg.front().ChrName(bamheader);
    // get denominator
    std::pair<int, int> q0den(0, 0);
    auto den = cpputil::CountDenom(seg, gr, ref, chrname, blacklist, errorstat, opt.bqual_min, q0den, true, opt.count_read, false);
    int r1_den = den.first;
    int r2_den = den.second;

    std::string aux_prefix = std::to_string(pair_nmismatch) + "\t" + std::to_string(nqpass) +
        "\t" + std::to_string(olen) + "\t" + std::to_string(abs(seg[0].InsertSize()));
    if (opt.count_read) {
      errorstat.neval += r1_den + r2_den;
    } else {
      errorstat.neval += r1_den;
    }
    errorstat.qcut_neval[0].first += q0den.first;
    errorstat.qcut_neval[0].second += q0den.second;
    if (opt.bqual_min > 0) {
      errorstat.qcut_neval[opt.bqual_min].first += r1_den;
      errorstat.qcut_neval[opt.bqual_min].second += r2_den;
    }
    bool fail_alignment_filter = false;

    if (!opt.cycle_level_stat.empty()) {
      CycleBaseCount(seg, gr, errorstat);
    }
    // first pass gets all variants in ROI
    std::map<cpputil::Variant, vector<cpputil::Variant>> var_vars;
    int rstart = 0;
    int rend = std::numeric_limits<int>::max();
    if (opt.count_read > 0) {
      if (gr) {
        rstart = gr->pos1;
        rend = gr->pos2;
      }
    } else {
      std::tie(rstart,rend) = cpputil::GetPairOverlapRStartAndRStop(seg.front(), seg.back());
      if (gr) {
        rstart = std::max(rstart, gr->pos1);
        rend = std::min(rend, gr->pos2);
      }
    }

    std::vector<cpputil::Variant> r1_vars, r2_vars;
    for (auto &s : seg) {
      auto vars = cpputil::GetVar(s, bamheader, ref);
      if (s.FirstFlag()) r1_vars = vars;
      else r2_vars = vars;
      for (auto &var : vars) {
        if (var.contig_start >= rstart && var.EndPos() <= rend)
          var_vars[var].push_back(var);
      }
    }

    std::vector<std::vector<cpputil::Variant>> refined_vars;
    refined_vars.reserve(var_vars.size());
    //Consolidate SNPs to break doublets and etc if not all bases pass bq filter.
    for (const auto& it: var_vars) {
      if (it.first.isIndel()) {
        refined_vars.push_back(it.second);
        continue;
      }
      if(it.second.size() == 1) {
        if (!it.first.isMNV() or it.second[0].var_qual >= opt.bqual_min) {
          refined_vars.push_back(it.second);
        } else {
          auto vars = var_atomize(it.second[0]);
          for (auto& var : vars) {
            refined_vars.push_back({var});
          }
        }
      } else {
        auto vars = cpputil::ConsolidateSNVpair(it.second[0], it.second[1], opt.bqual_min);
        refined_vars.insert(refined_vars.end(), vars.begin(), vars.end());
      }
    }
    // sort by var type first so that indels are in front
//    sort(refined_vars.begin(), refined_vars.end(), [](const auto& lhs, const auto& rhs) {
//      return lhs[0].Type() < rhs[0].Type();
//    });
    for (const auto& readpair_var: refined_vars) {

      auto var = cpputil::squash_vars(readpair_var);
      //int known_indel_len = 0;

      // alignment filter
      int XS;
      if (opt.max_frac_prim_AS < 1.0 and not seg[0].GetIntTag("XS", XS)) {
        int primary_score = 0, sec_as=0;
        //alignment filter
        std::string seq;
        if (seg.size() == 2) {
          seq = cpputil::MergePairSeq(seg, orig_seqs, false);
        } else {
          seq = seg[0].Sequence();
        }
        mem_alnreg_v ar;
        ar = mem_align1(bwa.GetMemOpt(), bwa.GetIndex()->bwt, bwa.GetIndex()->bns, bwa.GetIndex()->pac,
                        seq.length(), seq.data());
        if (ar.n >= 100) {
          fail_alignment_filter = true;
          free(ar.a);
          continue;
        }
        for (size_t idx = 0; idx < ar.n; ++idx) {
          if (ar.a[idx].secondary < 0) {
            primary_score = ar.a[idx].score;
            break;
          }
        }
        size_t idx = 0;
        for (;idx < ar.n; ++idx) {
          if (ar.a[idx].secondary >= 0) {
            sec_as = std::max(sec_as, ar.a[idx].score);
            if (ar.a[idx].score > primary_score * opt.max_frac_prim_AS) {
              break;
            }
          }
        }
        //std::cerr << seg[0].Qname() << "\t" << opt.max_frac_prim_AS << "\t" << primary_score <<"\t" << sec_as <<"\t" << ar.n << "\n";
        if (idx < ar.n) {
          fail_alignment_filter = true;
          free(ar.a);
          continue;
        }
        free(ar.a);
      }
      //string aux_output = aux_prefix + "\t" + std::to_string(primary_score) + "\t"  + std::to_string(sec_as);

      //pass alignment filters
      int nerr = 0;
//      std::vector<cpputil::Variant> avars;
//      if (n_true_mut(real_muts) > 0) { // partially wrong
//        auto tmpvar = cpputil::var_atomize(var);
//        for (unsigned ai = 0; ai < real_muts.size(); ++ai) {
//          if (not real_muts[ai]) {
//            avars.push_back(tmpvar[ai]);
//          }
//        }
//      } else { // fully wrong
//        avars.push_back(var);
//      }
      //R1R2 filter
      if (readpair_var.size() == 1) { // var in only one read
        if (opt.count_read == 0) {
          if (readpair_var[0].isIndel()) ++errorstat.indel_R1R2_disagree;
          else ++errorstat.snv_R1R2_disagree;
          continue;
        }
        if (opt.count_read == 1 && seg.size() == 2) {
          int ss =0, ee=0;
          std::tie(ss,ee) = cpputil::GetPairOverlapRStartAndRStop(seg.front(), seg.back());
          if (var.contig_start >= ss && var.contig_start < ee) {
            if (readpair_var[0].isIndel()) ++errorstat.indel_R1R2_disagree;
            else ++errorstat.snv_R1R2_disagree;
            continue;
          }
        }
      }

      int germ_support = 0, germ_depth = std::numeric_limits<int>::max();
      int site_depth = 0;
      std::pair<string, string> qtxt = {"NA", "NA"};
      if (var.isIndel() ) {// INDEL
        if (var.IndelLen() < opt.min_indel_len) {
          continue;
        }
        //INS seq should not contain 'N'
        if (var.contig_seq.find('N') != std::string::npos) {
          continue;
        }
        //INS at first or last base of a read may not be complete
        if (var.r1_start == -1 || var.r2_start == -1) {
          continue;
        }
        if (var.first_of_pair && var.r1_start + var.alt_seq.size() ==  orig_seqs[0].size())
          continue;
        if (var.second_of_pair && var.r2_start + var.alt_seq.size() == orig_seqs[1].size())
          continue;

        // INDEL near the end of a read
        int indel_near_r1_end = 0, indel_near_r2_end = 0;
        //std::cerr << var << ", " << var.r1_start << ", " << var.r2_start << ", " << orig_seqs[1].size() << std::endl;
        if (var.first_of_pair) {
          if (var.r1_start + opt.min_indel_dist_from_readend > (int) orig_seqs[0].size() || var.r1_start < opt.min_indel_dist_from_readend )
            indel_near_r1_end = 1;
          else
            indel_near_r1_end = -1;
        }
        if (var.second_of_pair) {
          if (var.r2_start + opt.min_indel_dist_from_readend > (int) orig_seqs[1].size() || var.r2_start < opt.min_indel_dist_from_readend )
            indel_near_r2_end = 1;
          else
            indel_near_r2_end = -1;
        }
        if ((indel_near_r1_end == 1 and indel_near_r2_end >= 0) or (indel_near_r2_end == 1 and indel_near_r1_end >= 0)) {
          ++errorstat.nindel_near_readend;
          continue;
        }


        if (var.first_of_pair) {
          int r1varend =  var.Type() == "INS" ? var.r1_start + var.alt_seq.length() : var.r1_start + 1;
          auto r1flank_f = orig_seqs[0].substr(std::max(var.r1_start - opt.indel_anchor_size + 1, 0), opt.indel_anchor_size);
          auto r1flank_b = orig_seqs[0].substr(r1varend, opt.indel_anchor_size);
          auto r1_pre_baseq = orig_qualities[0].substr(var.r1_start, 1);
          auto r1_suc_baseq = orig_qualities[0].substr(r1varend, 1);
          if (r1_pre_baseq[0] < opt.IndelAnchorBaseQ or r1_suc_baseq[0] < opt.IndelAnchorBaseQ )  {
            ++errorstat.nindel_filtered_adjbaseq;
            continue;
          }
          if (r1flank_f.find('N') != std::string::npos or r1flank_b.find('N') != std::string::npos) {
            ++errorstat.nindel_filtered_adjN;
            continue;
          }
        }
        if (var.second_of_pair) {
          int r2varend =  var.Type() == "INS" ? var.r2_start + var.alt_seq.length() : var.r2_start + 1;
          auto r2flank_f = orig_seqs[1].substr(std::max(var.r2_start - opt.indel_anchor_size + 1, 0), opt.indel_anchor_size);
          auto r2flank_b = orig_seqs[1].substr(r2varend, opt.indel_anchor_size);
          auto r2_pre_baseq = orig_qualities[1].substr(var.r2_start, 1);
          auto r2_suc_baseq = orig_qualities[1].substr(r2varend, 1);
//              std::cerr << "r2_start: " << var.r2_start << ", " << r2_pre_baseq <<", "<< r2_suc_baseq <<", " << seg[0].Qname() << std::endl;
          if (r2_pre_baseq[0] < opt.IndelAnchorBaseQ or r2_suc_baseq[0] < opt.IndelAnchorBaseQ )  {
            ++errorstat.nindel_filtered_adjbaseq;
            continue;
          }
          if (r2flank_f.find('N') != std::string::npos or r2flank_b.find('N') != std::string::npos) {
            ++errorstat.nindel_filtered_adjN;
            continue;
          }
        }

        if (bcf_reader.IsOpen()) {
          vector<bool> vcf_found = cpputil::search_var_in_database(bcf_reader,
                                                  var,
                                                  known,
                                                  "PRI_VCF",
                                                  true,
                                                  opt.bqual_min,
                                                  opt.all_mutant_frags);
          if (vcf_found[0]) {
            errorstat.nindel_masked_by_vcf1 += 1;
            continue;
          }
        }
        if (bcf_reader2.IsOpen()) {
          vector<bool> vcf_found = cpputil::search_var_in_database(bcf_reader2,
                                                  var,
                                                  known,
                                                  "SEC_VCF",
                                                  true,
                                                  opt.bqual_min,
                                                  opt.all_mutant_frags);
          if (vcf_found[0]) continue;
        }
        if (mafr.IsOpen()) {
          vector<bool> vcf_found = cpputil::search_var_in_database(mafr,
                                                  var,
                                                  known,
                                                  "MAF",
                                                  true,
                                                  opt.bqual_min,
                                                  opt.all_mutant_frags);
          if (vcf_found[0]) {
            errorstat.nindel_masked_by_maf += 1;
            continue;
          }
        }

        if (opt.germline_bam.size() > 1) {
          int exact_match = 0, fuzzy_match = 0;
          germ_depth = cpputil::ScanIndel(&pileup, true, ref, var.contig, var.contig_start,
                                          (int) var.alt_seq.size() - (int) var.contig_seq.size(),
                                          var.alt_seq.substr(1), true, opt.germline_var_maxdist,
                                          opt.germline_cutoff_vaf, 10, opt.germline_minbq, exact_match, fuzzy_match);
          if (exact_match > 0 or fuzzy_match > 1) {
            ++errorstat.seen_in_germ;
            continue;
          }
          if (germ_depth < opt.germline_mindepth) {
            ++errorstat.low_germ_depth;
            continue;
          }
        }
        //self pileup
        site_depth = cpputil::GenotypeVariant(&tumorpileup,
                                              ref,
                                              var.contig,
                                              var.contig_start,
                                              (int) var.alt_seq.size() - (int) var.contig_seq.size(),
                                              var.alt_seq.substr(1),
                                              true,1,
                                              opt.germline_cutoff_vaf, 10, opt.germline_minbq);
        if (site_depth == 0) {
          ++errorstat.seen_in_germ;
          continue;
        }
        if (site_depth == -1) {
          ++errorstat.nindel_filtered_overlap_snp;
          continue;
        }
        // search near by variant
        bool has_nearby = false;
        if (var.r1_start >0) {
          for (const auto &vvv : r1_vars) {
            if (vvv == var) continue;
            int cutoff = vvv.isIndel() ? opt.MIN_NEAREST_INDEL : opt.MIN_NEAREST_SNV;
            if (abs(vvv.r1_start - var.r1_start) < cutoff) {
              has_nearby = true;
              break;
            }
          }
        }
        if (has_nearby) {
          ++errorstat.nindel_filtered_adjvar;
          continue;
        }
        if (var.r2_start >0) {
          for (const auto &vvv : r2_vars) {
            if (vvv == var) continue;
            int cutoff = vvv.isIndel() ? opt.MIN_NEAREST_INDEL : opt.MIN_NEAREST_SNV;
            if (abs(vvv.r2_start - var.r2_start) < cutoff) {
              has_nearby = true;
              break;
            }
          }
        }
        if (has_nearby) {
          ++errorstat.nindel_filtered_adjvar;
          continue;
        }
        ++errorstat.nindel_error;
        errorstat.indel_nbase_error += abs((int) var.contig_seq.length() - (int) var.alt_seq.length());
        if (opt.count_read == 0 ) {
          qtxt = cpputil::QualContext(var, seg, 3);
        }
        ferr << var << '\t' << aux_prefix << '\t' << qtxt.first << '\t' <<qtxt.second << "\t" << site_depth << "\t" << germ_depth << '\n';
      } else { // SNV
        if (var.var_qual < opt.bqual_min * var.read_count) {
          errorstat.snv_filtered_baseq += var.alt_seq.size();
          continue;
        }

        if (bcf_reader.IsOpen()) {
          vector<bool> found = cpputil::search_var_in_database(bcf_reader,
                                                               var,
                                                               known,
                                                               "PRI_VCF",
                                                               true,
                                                               opt.bqual_min,
                                                               opt.all_mutant_frags);
          if (n_true_mut(found)> 0) {
            errorstat.nsnv_masked_by_vcf1 += var.alt_seq.size();
            continue;
          }
        }
        if (bcf_reader2.IsOpen()) {
          vector<bool> found = cpputil::search_var_in_database(bcf_reader2,
                                                               var,
                                                               known,
                                                               "SEC_VCF",
                                                               true,
                                                               opt.bqual_min,
                                                               opt.all_mutant_frags);
          if (n_true_mut(found)> 0) continue;
        }
        if (mafr.IsOpen()) {
          vector<bool> found = cpputil::search_var_in_database(mafr,
                                                               var,
                                                               known,
                                                               "MAF",
                                                               true,
                                                               opt.bqual_min,
                                                               opt.all_mutant_frags);
          if (n_true_mut(found) > 0) {
            errorstat.nsnv_masked_by_maf += 1;
            continue;
          }
        }
        //germline filter
        if (opt.germline_bam.size() > 1) {
          germ_depth = cpputil::ScanAllele(&pileup, var.contig, var.contig_start, var.alt_seq[0], true, germ_support, opt.germline_minbq);
        }
        if (germ_support >= opt.germline_minalt) {
          ++errorstat.seen_in_germ;
          //std::cerr << var << "\tseen in germ " << germ_support << std::endl;
          continue;
        }
        if (germ_depth < opt.germline_mindepth) {
          ++errorstat.low_germ_depth;
          //std::cerr << var << "\tunder covered in germ " << germ_support << std::endl;
          continue;
        }

        //selfpileup
        site_depth = cpputil::GenotypeVariant(&tumorpileup,
                                              ref,
                                              var.contig,
                                              var.contig_start,
                                              0,
                                              var.alt_seq,
                                              true,0,
                                              opt.germline_cutoff_vaf,
                                              10,
                                              opt.germline_minbq);
        if (site_depth == 0) {
          ++errorstat.seen_in_germ;
          continue;
        }

        //filter if family size > 1
        bool snv_family_disagree = false;
        for (unsigned ss = 0; ss < seg.size(); ss++) {
          if (seg[ss].FirstFlag()) {
            if (not var.first_of_pair) continue;
          } else {
            if (not var.second_of_pair) continue;
          }
          int32_t fsize = 1;
          bool cDstat = seg[ss].GetIntTag("cD", fsize);
          //std::cerr << "cD " << fsize << std::endl;
          if (cDstat and fsize > 1) {
            vector<int64_t> fam_depth;
            vector<int64_t> fam_error;
            bool has_depth = cpputil::GetBTag(seg[ss], "cd", fam_depth);
            if (has_depth) {
              int readpos = seg[ss].FirstFlag() ? var.r1_start : var.r2_start;
              //std::cerr << seg[ss] << " read pos " << readpos << std::endl;
              if (fam_depth[readpos] < round(fsize * opt.family_agree_rate)) {
                std::cerr << "family not agree on depth " << seg[ss].Qname() << "\t" << seg[ss].Position() << "\t" << readpos << "\t" << std::endl;
                snv_family_disagree = true;
                break;
              }
            }
            bool has_error = cpputil::GetBTag(seg[ss], "ce", fam_error);
            if (has_error) {
              int readpos = seg[ss].FirstFlag() ? var.r1_start : var.r2_start;
              if (fam_error[readpos] > round((1-opt.family_agree_rate) * fsize)) {
                std::cerr << "family not agree on bases " << seg[ss].Qname() << "\t" << seg[ss].Position() << "\t" << readpos << "\t" << std::endl;
                snv_family_disagree = true;
                break;
              }
            }
          }
        }
        if (snv_family_disagree) {
          ++errorstat.snv_family_disagree;
          continue;
        }


        if (not opt.allow_indel_near_snv && (cpputil::IndelLen(seg.front())> 0 ||
                    (seg.size() == 2 && cpputil::IndelLen(seg.back()) > 0))) {
          ++errorstat.mismatch_filtered_by_indel;
          continue;
        }
        if (var.MutType() == "T>G" and opt.min_passQ_frac_TT > 0) {
              if (var.DoubletContext(ref) == "TT" and nqpass < opt.min_passQ_frac_TT) {
                ++errorstat.lowconf_t2g;
                continue;
              }
        }
        nerr += var.alt_seq.size();
        if (opt.count_read == 0) {
          qtxt = cpputil::QualContext(var, seg, 3);
        }
        ferr << var << '\t' << aux_prefix << '\t' <<  qtxt.first << '\t' <<qtxt.second << "\t" << site_depth << "\t" << germ_depth << '\n';
      }
      if (opt.count_read == 2) {
        errorstat.nsnv_error += nerr * var.read_count;
      } else {
        errorstat.nsnv_error += nerr;
      }

      // other error profiles
      if (var.Type() == "SNV") {
        if (!opt.read_level_stat.empty()) {
          if (var.var_qual >= opt.bqual_min * var.read_count && (seg.size() == 1 || readpair_var.size() == 2)) {
            if (var.read_count == 2) {
              r1_nerror += nerr;
              r2_nerror += nerr;
            }
            else {
              if (var.first_of_pair ) r1_nerror += nerr;
              else r2_nerror += nerr;
            }
          }
        }
        //per cycle profile
        if (!opt.cycle_level_stat.empty()) {
          int tpos = var.alt_start < 0 ?  abs(var.alt_start + 1): var.alt_start;
          var.first_of_pair ? errorstat.R1_q0_error[tpos] += nerr : errorstat.R2_q0_error[tpos] += nerr;
          var.second_of_pair ? errorstat.R2_q0_error[tpos] += nerr : errorstat.R1_q0_error[tpos] += nerr;
          if (var.var_qual >= opt.bqual_min * var.read_count) {
            var.first_of_pair ? errorstat.R1_q30_error[tpos] += nerr: errorstat.R2_q30_error[tpos] += nerr;
            var.second_of_pair ? errorstat.R2_q30_error[tpos] += nerr: errorstat.R1_q30_error[tpos] += nerr;
          }
        }
      } //end other profiles
    }
    if (fail_alignment_filter) {
      ++errorstat.AS_filter;
    }

    if (readlevel.is_open()) {
      readlevel << seg.front().Qname() << '\t'
                << r1_nerror << '\t'
                << r2_nerror << '\t'
                << r1_den  << '\t'
                << r2_den << '\t'
                << pair_nmismatch << '\t'
                << abs(seg.front().InsertSize()) << '\t'
                << olen << '\n';
    }
  } // end for
  return;
}

int codec_accuracy(int argc, char ** argv) {
  string cmdline;
  for(int i=0; i<argc; i++){
    cmdline += argv[i];
    cmdline+=" ";
  }
  AccuOptions opt;
  int parse_ret =  accuracy_parse_options(argc, argv, opt);
  if (parse_ret) return 1;
  if (argc == 1) {
    accuracy_print_help();
    return 1;
  }
  if (not opt.preset.empty()) {
    if (opt.preset == "lenient") {
      opt.mapq = 60;
      opt.bqual_min = 30;
      opt.fragend_dist_filter = 12;
      opt.min_passQ_frac = 0.5;
      opt.filter_5endclip = true;
      opt.max_snv_filter = 10;
      opt.max_frac_prim_AS = 0.6;
      opt.germline_minmapq = 60;
      opt.max_pair_mismatch_frac = 0.1;
    } else if(opt.preset == "stringent") {
      opt.mapq = 60;
      opt.bqual_min = 30;
      opt.fragend_dist_filter = 12;
      opt.min_passQ_frac = 0.7;
      opt.filter_5endclip = true;
      opt.max_snv_filter = 4;
      opt.clustered_mut_cutoff = 3;
      opt.max_frac_prim_AS = 0.5;
      opt.germline_minmapq = 60;
      opt.max_pair_mismatch_frac = 0.03;
    } else if(opt.preset == "null") {
      opt.mapq = 0;
      opt.bqual_min = 0;
      opt.min_fraglen = 0;
      opt.germline_mindepth = 0;
      opt.germline_var_maxdist = 0;
      opt.allow_indel_near_snv = true;
      opt.min_indel_dist_from_readend = 0;
      opt.family_agree_rate = 0.0;
      opt.MIN_NEAREST_SNV = 0;
      opt.MIN_NEAREST_INDEL = 0;
    } else {
      std::cerr << "preset must by one of the [lenient, stringent and null]\n";
      exit(1);
    }
  }
  SeqLib::RefGenome ref;
  ref.LoadIndex(opt.reference);
  if (opt.mutation_metrics.empty()) {
    if (opt.sample.empty()) {
      std::cerr << "must specify output. -a or -o is not used\n";
      exit(0);
    }
    opt.mutation_metrics = opt.sample + ".mutation_metrics.txt";
  }
  if (opt.variants_out.empty()) {
    if (opt.sample.empty()) {
      std::cerr << "must specify output. -e or -o is not used\n";
      exit(0);
    }
    opt.variants_out = opt.sample + ".variants_called.txt";
  }
  if (opt.context_count.empty()) {
    if (opt.sample.empty()) {
      std::cerr << "must specify output. -C or -o is not used\n";
      exit(0);
    }
    opt.context_count = opt.sample + ".context_count.txt";
  }
  std::ofstream stat(opt.mutation_metrics);
  std::ofstream ferr(opt.variants_out);
  std::ofstream context(opt.context_count);
  if (!stat.is_open()) {
    std::cerr << opt.mutation_metrics << " cannot be opened\n";
    exit(1);
  }
  if (!ferr.is_open()) {
    std::cerr << opt.variants_out << " cannot be opened\n";
    exit(1);
  }
  if (!context.is_open()) {
    std::cerr << opt.context_count << " cannot be opened\n";
    exit(1);
  }

  cpputil::BCFReader bcf_reader;
  cpputil::BCFReader bcf_reader2;
  if (!opt.vcfs.empty()) {
    bcf_reader.Open(opt.vcfs[0].c_str());
    if (opt.vcfs.size() == 2) {
      bcf_reader2.Open(opt.vcfs[1].c_str());
    }
    else if (opt.vcfs.size() > 2){
      std::cerr << "Cannot read more than 2 vcfs currently\n";
      return 1;
    }
  }
  SeqLib::BWAWrapper  bwa;
  bwa.LoadIndex(opt.reference);
  cpputil::PileHandler pileup;
  if (opt.germline_bam.size() > 1) {
    pileup = cpputil::PileHandler(opt.germline_bam, opt.germline_minmapq);
  }
  std::ofstream known;
  std::ofstream readlevel;
  std::ofstream cyclelevel;
  string error_profile_header =
      "chrom\tref_pos\tref\talt\ttype\tdist_to_fragend\tsnv_base_qual\tread_count\tread_name\tfamily_size\tnumN\tnumQpass\tclen\tflen\tqual3p\tqual5p\tsite_depth\tgerm_depth";
  ferr << "#" << cmdline << std::endl;
  ferr << error_profile_header << std::endl;
  if (not opt.known_var_out.empty()) {
    known.open(opt.known_var_out);
    string known_var_header =
        "chrom\tref_pos\tref\talt\ttype\tread_pos\tsnv_base_qual\tfirst_of_pair\tread_name\tevidence";
    known << known_var_header << std::endl;
  }
  if (not opt.read_level_stat.empty()) {
    readlevel.open(opt.read_level_stat);
    string read_level_header =
        "read_name\tR1_nerror\tR2_nerror\tR1_efflen\tR2_efflen\tpair_nmismatch\tfraglen\toverlap_len";
    readlevel << read_level_header << std::endl;
  }
  if (not opt.cycle_level_stat.empty()) {
    cyclelevel.open(opt.cycle_level_stat);
    cyclelevel << "cycle\tR1_q0_error\tR1_q0_cov\tR1_q30_error\tR1_q30_cov\t"
            "R2_q0_error\tR2_q0_cov\tR2_q30_error\tR2_q30_cov\t"
            "R1_q0_erate\tR1_q30_erate\tR2_q0_erate\tR2_q30_erate\n";
  }
  cpputil::InsertSeqFactory isf(opt.bam,
                                10,
                                opt.load_supplementary,
                                opt.load_secondary,
                                opt.load_duplicate,
                                opt.load_proper_pair_only,
                                false);

  cpputil::PileHandler selfpile(opt.bam, 30);
  cpputil::MAFReader mafr;
  if (!opt.maf_file.empty()) {
    mafr.Open(opt.maf_file);
  }
  const int L = 250;
  auto errorstat = cpputil::ErrorStat(L, opt.bqual_min);
  if (opt.detail_qscore_prof) {
    errorstat = cpputil::ErrorStat({0,10,20,30}, L);
  }

  cpputil::TargetLayout tl(isf.bamheader(), opt.bed_file);
  // SeqLib::GenomicRegion gr;
  //while(tl.NextRegion(gr)) {
  vector<vector<cpputil::Segments>> chunk;
  std::unordered_set<std::string> pass_qnames;
  int last_chrom = -1;
  int last_end = 0;
  cpputil::UniqueQueue readpair_cnt(100000);
  cpputil::UniqueQueue failed_cnt(100000);
  for (unsigned i = 0; i < tl.NumRegion(); ++i) {
    const auto& gr = tl[i];
    if (gr.chr != last_chrom || gr.pos1 - last_end > 1e4) {
      last_chrom = gr.chr;
      readpair_cnt.clearQueue();
      failed_cnt.clearQueue();
    }
    last_end = gr.pos2;
    if (i % 100000 == 0) {
      auto timenow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
      std::cerr << i + 1 << " region processed. Last position: " << gr << std::ctime(&timenow) << std::endl;
    }
    if (isf.ReadByRegion(cpputil::ArePEAlignmentOverlapAtLeastK, gr, chunk, opt.pair_min_overlap, opt.MID_tag, opt.load_unpair)) {
      std::set<int32_t> blacklist;
      if (bcf_reader.IsOpen()) {
        bcf_reader.vcf_to_blacklist_snv(gr, isf.bamheader(), blacklist);
      }
      for (vector<cpputil::Segments>& frag : chunk) { //frag will be a family if -D is true
        //std::cerr << frag[0][0] << " load" << std::endl;
        if (opt.min_indel_len) {
          bool has_indel = false;
          for (const auto &seg: frag) {
            for (const auto &r : seg)
              if (cpputil::IndelLen(r) >= opt.min_indel_len) {
                has_indel = true;
                break;
              }
          }
          if (not has_indel) continue;
        }
        for (const auto& ff : frag) {
          readpair_cnt.add(ff[0].Qname());
        }
        if (failed_cnt.exist(frag[0][0].Qname())) {
          //std::cerr << frag[0][0] << " already faile" << std::endl;
          continue;
        }

        int pair_nmismatch = 0, olen = 0;
        float nqpass;
        int fail = cpputil::FailFilter(frag, isf.bamheader(), ref, opt, errorstat, isf.IsPairEndLib() & !opt.load_unpair, pair_nmismatch, nqpass, olen);
        if (fail) {
          failed_cnt.add(frag[0][0].Qname());
          errorstat.discard_frag_counter += frag.size();
          //std::cerr << fail << ": " << frag[0].size() << ", " << frag[0][0].Qname() << std::endl;
        }
        else {
          ErrorRateDriver(frag,
                          pair_nmismatch, nqpass, olen,
                          isf.bamheader(),
                          ref,
                          bwa,
                          &gr,
                          blacklist,
                          opt,
                          bcf_reader,
                          bcf_reader2,
                          mafr,
                          pileup,
                          selfpile,
                          ferr,
                          known,
                          readlevel,
                          errorstat);
        }
      } // end for
    } //end if
  } //end for

  std::cerr << "All region processed \n";

  vector<string> header = {"qcutoff",
                           "n_trim",
                           "n_bases_eval",
                           "n_A_eval",
                           "n_C_eval",
                           "n_G_eval",
                           "n_T_eval",
                           "n_snv",
                           "snv_rate",
                           "n_indel",
                           "indel_rate",
                           "n_indel_bases",
                           "indel_bases_rate",
                           "n_totalfrag",
                           "n_pass_filter_pair_seg",
                           "n_pass_filter_single_seg",
                           "n_filtered_total",
                           "n_snv_masked_by_vcf1",
                           "n_indel_masked_by_vcf1",
                           "n_snv_masked_by_maf",
                           "n_indel_masked_by_maf",
                           "n_filtered_smallfrag",
                           "n_filtered_mapq",
                           "n_filtered_scalip",
                           "n_filtered_pairmismatch_rate",
                           "n_filtered_passBQ_rate",
                           "n_filtered_numN",
                           "n_filtered_largefrag",
                           "n_filtered_edit",
                           "n_filtered_clustered",
                           "n_filtered_AS",
                           "n_filtered_badcigar",
                           "mut_germ_lowdepth",
                           "mut_germ_seen",
                           "nsnv_R1R2_disagree",
                           "nsnv_low_baseq",
                           "nsnv_family_disagree",
                           "nsnv_filtered_by_indel",
                           "nsnv_t2g_low_conf",
                           "nindel_R1R2_disagree",
                           "nindel_filtered_by_near_readend",
                           "nindel_filtered_by_adjbaeq",
                           "nindel_filtered_by_adjN",
                           "nindel_filtered_by_adjvar",
                           "nindel_filtered_by_overlap_SNP",
                           "qpass_rate"};

  //for (const auto& duo : errorstat.qcut_nerrors) {
//    header.push_back("all_n_bases_eval" );
//    header.push_back("all_n_errors" );
//    header.push_back("all_erate" );
    header.push_back("q0_R1_n_bases_eval" );
    header.push_back("q0_R2_n_bases_eval" );
  //}
  stat << cpputil::join(header, "\t") << std::endl;
  stat << std::to_string(opt.bqual_min) + "/" + std::to_string(opt.count_read)  << '\t'
       << opt.fragend_dist_filter << '\t'
       << errorstat.neval << '\t'
      << errorstat.base_counter['A'] << '\t'
      << errorstat.base_counter['C'] << '\t'
      << errorstat.base_counter['G'] << '\t'
      << errorstat.base_counter['T'] << '\t'
      << errorstat.nsnv_error << '\t'
       << (float) errorstat.nsnv_error / errorstat.neval << '\t'
      << errorstat.nindel_error << '\t'
      << (float) errorstat.nindel_error / errorstat.neval << '\t'
      << errorstat.indel_nbase_error << '\t'
      << (float) errorstat.indel_nbase_error / errorstat.neval << '\t'
       << readpair_cnt.NumAdded() << '\t'
       << errorstat.n_pass_filter_pairs << '\t'
       << errorstat.n_pass_filter_singles << '\t'
       << errorstat.discard_frag_counter << '\t'
      << errorstat.nsnv_masked_by_vcf1 << '\t'
      << errorstat.nindel_masked_by_vcf1 << '\t'
      << errorstat.nsnv_masked_by_maf << '\t'
      << errorstat.nindel_masked_by_maf << '\t'
       << errorstat.n_filtered_smallfrag << '\t'
      << errorstat.n_filtered_by_mapq << '\t'
      << errorstat.n_filtered_sclip << '\t'
      << errorstat.n_filtered_pairmismatch_rate << '\t'
      << errorstat.n_filtered_q30rate << '\t'
      << errorstat.n_filtered_numN << '\t'
       << errorstat.n_filtered_largefrag << '\t'
       << errorstat.n_filtered_edit << '\t'
      << errorstat.n_filtered_clustered << '\t'
      << errorstat.AS_filter << '\t'
      << errorstat.n_filtered_badcigar << '\t'
      << errorstat.low_germ_depth << '\t'
       << errorstat.seen_in_germ << '\t'
      << errorstat.snv_R1R2_disagree << '\t'
      << errorstat.snv_filtered_baseq << '\t'
      << errorstat.snv_family_disagree << '\t'
      << errorstat.mismatch_filtered_by_indel << '\t'
       << errorstat.lowconf_t2g << '\t'
      << errorstat.indel_R1R2_disagree << '\t'
      << errorstat.nindel_near_readend << '\t'
      << errorstat.nindel_filtered_adjbaseq << '\t'
      << errorstat.nindel_filtered_adjN << '\t'
      << errorstat.nindel_filtered_adjvar << '\t'
      << errorstat.nindel_filtered_overlap_snp << '\t'
       << (float) (errorstat.qcut_neval[opt.bqual_min].first +  errorstat.qcut_neval[opt.bqual_min].second) /
       (errorstat.qcut_neval[0].first + errorstat.qcut_neval[0].second)  << '\t';

  unsigned ii = 0;
  for (const auto& qcut : errorstat.cutoffs) {
    ++ii;
    int64_t r1_den = errorstat.qcut_neval[qcut].first;
    int64_t r2_den = errorstat.qcut_neval[qcut].second;
    stat << r1_den << '\t';
    stat << r2_den;
    if (ii == errorstat.cutoffs.size())  stat << '\n';
    else stat << '\t';
  }

  for (const auto& it : errorstat.triplet_counter) {
    context << it.first << "\t" << it.second << std::endl;
  }

  for (const auto& it : errorstat.doublet_counter) {
    context << it.first << "\t" << it.second << std::endl;
  }

  if (!opt.cycle_level_stat.empty()) {
    for (int i = 0; i < L; ++i) {
     cyclelevel  << i << '\t';
     cyclelevel  << errorstat.R1_q0_error[i] << '\t' << errorstat.R1_q0_cov[i] << '\t';
     cyclelevel  << errorstat.R1_q30_error[i] << '\t' << errorstat.R1_q30_cov[i] << '\t';
     cyclelevel  << errorstat.R2_q0_error[i] << '\t' << errorstat.R2_q0_cov[i] << '\t';
     cyclelevel  << errorstat.R2_q30_error[i] << '\t' << errorstat.R2_q30_cov[i] << '\t';
     cyclelevel  << (float) errorstat.R1_q0_error[i] / errorstat.R1_q0_cov[i] << '\t';
     cyclelevel  << (float) errorstat.R1_q30_error[i] / errorstat.R1_q30_cov[i] << '\t';
     cyclelevel  << (float) errorstat.R2_q0_error[i] / errorstat.R2_q0_cov[i] << '\t';
     cyclelevel  << (float) errorstat.R2_q30_error[i] / errorstat.R2_q30_cov[i] << '\n';
    }
  }
  return 0;
}
