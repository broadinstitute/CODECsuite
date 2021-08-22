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

//#include "SeqLib/BWAWrapper.h"
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
#include "BamRecordExt.h"
#include "Alignment.h"
#include "StringUtils.h"
#include "MAF.h"
#include "TargetLayout.h"
#include "MutCounter.h"

#define OPT_QSCORE_PROF   261
#define OPT_READ_LEVEL_STAT   262
#define OPT_CYCLE_LEVEL_STAT   263
using std::string;
using std::vector;
int MYINT_MAX = std::numeric_limits<int>::max();

struct AccuOptions {
  string bam;
  string accuracy_stat = "accuracy_stat.txt";
  string error_prof_out = "error_prof_out.txt";
  int mapq = 10;
  int bqual_min = 0;
  bool load_supplementary = false;
  bool load_secondary = false;
  bool load_unpair = false;
  bool load_duplicate = false;
  bool all_mutant_frags = false;
  int max_nm_filter = MYINT_MAX;
  int max_mismatch_filter = MYINT_MAX;
  int verbose = 0;
  vector<string> vcfs;
  string sample = "";
  string reference;
  int fragend_dist_filter = 0;
  string known_var_out = "known_var_out.txt";
  string context_count = "context_count.txt";
  bool detail_qscore_prof = false;
  string read_level_stat;
  string cycle_level_stat;
  //double min_frac_bqual_pass = 0.1;
  int pair_min_overlap = 0;
  bool count_nonoverlap = false;
  string maf_file;
  string bed_file;
  string trim_bam = "trimmed.bam";
};


static struct option  accuracy_long_options[] = {
    {"bam",                      required_argument,      0,        'b'},
    {"bed",                      required_argument,      0,        'L'},
    {"accuracy_stat",            required_argument,      0,        'a'},
    {"load_unpair",              no_argument,            0,        'u'},
    {"load_supplementary",       no_argument,            0,        'S'},
    {"load_secondary",           no_argument,            0,        '2'},
    {"load_duplicate",           no_argument,            0,        'D'},
    {"mapq",                     required_argument ,     0,        'm'},
    {"vcfs",                     required_argument ,     0,        'V'},
    {"maf",                      required_argument ,     0,        'M'},
    {"sample",                   required_argument,      0,        's'},
    {"reference",                required_argument,      0,        'r'},
    {"error_prof_out",           required_argument,      0,        'e'},
    {"bqual_min",                required_argument,      0,        'q'},
    //{"min_frac_bqual_pass",      required_argument,      0,        'F'},
    {"known_var_out",            required_argument,      0,        'k'},
    {"context_count",            required_argument,      0,        'C'},
    {"pair_min_overlap",         required_argument,      0,        'p'},
    {"count_nonoverlap",         no_argument,            0,        'c'},
    {"fragend_dist_filter",      required_argument,      0,        'd'},
    {"max_nm_filter",            required_argument,      0,        'n'},
    {"max_mismatch_filter",      required_argument,      0,        'x'},
    {"verbose",                  required_argument,      0,        'v'},
    {"all_mutant_frags",         no_argument,             0,        'A'},
    {"detail_qscore_prof",       no_argument,            0,        OPT_QSCORE_PROF},
    {"read_level_stat",          required_argument,      0,        OPT_READ_LEVEL_STAT},
    {"cycle_level_stat",         required_argument,      0,        OPT_CYCLE_LEVEL_STAT},
    {0,0,0,0}
};

const char* accuracy_short_options = "b:a:m:v:S2us:r:e:q:k:cM:p:d:n:x:V:L:DAC:";

void accuracy_print_help()
{
  std::cerr<< "---------------------------------------------------\n";
  std::cerr<< "Usage: codec accuracy [options]\n";
  std::cerr<< "General Options:\n";
  std::cerr<< "-v/--verbose,                          [default 0]\n";
  std::cerr<< "-b/--bam,                              input bam\n";
  std::cerr<< "-L/--bed,                              targeted region\n";
  std::cerr<< "-m/--mapq,                             min mapping quality [10].\n";
  std::cerr<< "-S/--load_supplementary,               include supplementary alignment [false].\n";
  std::cerr<< "-2/--load_secondary,                   include secondary alignment [false].\n";
  std::cerr<< "-u/--load_unpair,                      include unpaired alignment [false].\n";
  std::cerr<< "-D/--load_duplicate,                   include alignment marked as duplicates [false].\n";
  std::cerr<< "-V/--vcfs,                             comma separated VCF file(s) for germline variants or whitelist variants[null].\n";
  std::cerr<< "-M/--maf,                              MAF file for somatic variants [null].\n";
  std::cerr<< "-s/--sample,                           sample from the VCF file [null].\n";
  std::cerr<< "-r/--reference,                        reference sequence in fasta format [null].\n";
  std::cerr<< "-a/--accuracy_stat,                    output reporting accuracy for each alignment [accuracy_stat.txt].\n";
  std::cerr<< "-e/--error_prof_out,                   Error profile output in plain txt format [error_prof_out.txt].\n";
  std::cerr<< "-k/--known_var_out,                    Output for known var. [known_var_out.txt].\n";
  std::cerr<< "-C/--context_count,                    Output for trinucleotide and dinucleotide context context in the reference. [context_count.txt].\n";
  //std::cerr<< "--qscore_prof,                         Output qscore prof. First column is qscore cutoff; second column is number of bases in denominator\n";
  std::cerr<< "--detail_qscore_prof,                  Output finer scale qscore cutoffs, error rates profile. The default is only q0, q30 [false]. \n";
  std::cerr<< "--read_level_stat,                     Output read level error metrics.\n";
  std::cerr<< "--cycle_level_stat,                    Output cycle level error metrics.\n";
  std::cerr<< "\nFiltering Options:\n";
  std::cerr<< "-q/--bqual_min,                        Skip bases if min(q1, q2) < this when calculating error rate. q1, q2 are baseQ from R1 and R2 respectively [0].\n";
  std::cerr<< "-n/--max_edit_filter,                  Skip a read if its NM tag is larger than this value [INT_MAX].\n";
  std::cerr<< "-x/--max_nonNedit_filter,              Skip a read if the number of non-N bases edits is larger than this value [INT_MAX].\n";
  std::cerr<< "-d/--fragend_dist_filter,              Consider a variant if its distance to the fragment end is at least this value [0].\n";
  //std::cerr<< "-F/--min_frac_bqual_pass,              If bqual_min is none 0. The number of bases passing bqual filter has to be larger than or equal to this frac [0.1].\n";
  std::cerr<< "-p/--pair_min_overlap,                 When using selector, the minimum overlap between the two ends of the pair. -1 for complete overlap, 0 no overlap required [0].\n";
  std::cerr<< "-c/--count_nonoverlap,                 Count overhang (non overlapping) region of a read pair. Default [false].\n";
  std::cerr<< "-A/--all_mutant_frags,                  Output all mutant fragments even if failed by concordant or qscore filters Default [false].\n";
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
      case 'L':
        opt.bed_file = optarg;
        break;
      case 'a':
        opt.accuracy_stat = optarg;
        break;
      case 'm':
        opt.mapq = atoi(optarg);
        break;
      case 'q':
        opt.bqual_min = atoi(optarg);
        break;
      case 'n':
        opt.max_nm_filter = atoi(optarg);
        break;
      case 'x':
        opt.max_mismatch_filter = atoi(optarg);
        break;
      case 'd':
        opt.fragend_dist_filter = atoi(optarg);
        break;
      case 'S':
        opt.load_supplementary = true;
        break;
      case '2':
        opt.load_secondary = true;
        break;
      case 'u':
        opt.load_unpair = true;
        break;
      case 'D':
        opt.load_duplicate = true;
      case 'c':
        opt.count_nonoverlap = true;
        break;
      case 's':
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
        opt.error_prof_out = optarg;
        break;
      case 'k':
        opt.known_var_out = optarg;
        break;
      case 'C':
        opt.context_count = optarg;
        break;
      case OPT_QSCORE_PROF:
        opt.detail_qscore_prof = true;
        break;
      case OPT_READ_LEVEL_STAT:
        opt.read_level_stat = optarg;
        break;
      case OPT_CYCLE_LEVEL_STAT:
        opt.cycle_level_stat = optarg;
        break;
      case 'p':
        opt.pair_min_overlap = atoi(optarg);
        break;
      case 'A':
        opt.all_mutant_frags = true;
        break;
      default:accuracy_print_help();
        return 1;
    }
  } while (next_option != -1);

  return 0;
}


std::pair<int,int> NumEffBases(const cpputil::Segments& seg,
    const SeqLib::GenomicRegion* const gr,
    const SeqLib::RefGenome& ref,
    const string& chrname,
    const std::set<int>& blacklist,
    cpputil::ErrorStat& es,
    int minbq,
    bool count_context,
    bool count_overhang){
  // blacklist represents a set of SNV positions that will not be counted in the error rate calculation
  int r1;
  int r2;
  std::set<int> baseqblack;
  if (seg.size() == 1) {
    std::pair<int,int> range;
    if (gr) {
      range = cpputil::GetBamOverlapQStartAndQStop(seg.front(), *gr);
    } else {
      range.first = seg.front().AlignmentPosition();
      range.second = seg.front().AlignmentEndPosition();
    }
    if (seg.front().FirstFlag()) {
      if (count_context)
        r1 = cpputil::CountValidBaseAndContextInMatchedBases(seg.front(), blacklist, chrname, ref, minbq, baseqblack, es, range.first, range.second);
      else
        r1 = cpputil::CountValidBaseInMatchedBases(seg.front(), blacklist, minbq, baseqblack,range.first, range.second);
    }
    else {
      if (count_context)
        r2 = cpputil::CountValidBaseAndContextInMatchedBases(seg.front(), blacklist, chrname, ref, minbq, baseqblack, es, range.first, range.second);
      else
        r2 = cpputil::CountValidBaseInMatchedBases(seg.front(), blacklist, minbq, baseqblack, range.first, range.second);
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
    std::pair<int, int> range_front, range_back;
    if (count_overhang) {
      range_front = overlap_front;
      range_back = overlap_back;
    } else {
      std::tie(range_front, range_back) = cpputil::GetPairOverlapQStartAndQStop(seg.front(), seg.back());
      range_front.first = std::max(range_front.first, overlap_front.first);
      range_front.second = std::min(range_front.second, overlap_front.second);
      range_back.first = std::max(range_back.first, overlap_back.first);
      range_back.second = std::min(range_back.second, overlap_back.second);
    }
    if (seg.front().FirstFlag()) {
      if (count_context) {
        r1 = cpputil::CountValidBaseAndContextInMatchedBases(seg.front(), blacklist, chrname, ref, minbq, baseqblack, es, range_front.first, range_front.second);
        r2 = cpputil::CountValidBaseAndContextInMatchedBases(seg.back(), blacklist, chrname, ref, minbq, baseqblack, es, range_back.first, range_back.second);
      }
      else {
        r1 = cpputil::CountValidBaseInMatchedBases(seg.front(), blacklist, minbq, baseqblack, range_front.first, range_front.second);
        r2 = cpputil::CountValidBaseInMatchedBases(seg.back(), blacklist, minbq, baseqblack, range_back.first, range_back.second);

      }
    } else {
      if (count_context) {
        r2 = cpputil::CountValidBaseAndContextInMatchedBases(seg.front(), blacklist, chrname, ref, minbq, baseqblack, es, range_front.first, range_front.second);
        r1 = cpputil::CountValidBaseAndContextInMatchedBases(seg.back(), blacklist, chrname, ref, minbq, baseqblack, es, range_back.first, range_back.second);
      }
      else {
        r2 = cpputil::CountValidBaseInMatchedBases(seg.front(), blacklist, minbq, baseqblack, range_front.first, range_front.second);
        r1 = cpputil::CountValidBaseInMatchedBases(seg.back(), blacklist, minbq,  baseqblack, range_back.first, range_back.second);
      }
    }
  }
  return std::make_pair(r1, r2);
}

void CycleBaseCount(const cpputil::Segments& seg,
                    const SeqLib::GenomicRegion* const gr,
                    cpputil::ErrorStat& es) {
  auto res = std::make_pair(0,0);
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

void ErrorRateDriver(const vector<cpputil::Segments>& frag,
                    const SeqLib::BamHeader& bamheader,
                    const SeqLib::RefGenome& ref,
                    const SeqLib::GenomicRegion* const gr,
                    const std::set<int> blacklist,
                    const AccuOptions& opt,
                    const cpputil::BCFReader& bcf_reader,
                    const cpputil::BCFReader& bcf_reader2,
                    const cpputil::MAFReader& mafr,
                    std::ofstream& ferr,
                    std::ofstream& known,
                    std::ofstream& readlevel,
                    cpputil::ErrorStat& errorstat) {
  for (auto seg : frag) { // a fragment may have multiple segs due to supplementary alignments
    if (seg.size() > 2) {
      throw std::runtime_error("seg size cannot be more than 2");
    }
  //          if (seg.size() == 2 && seg.front().ReverseFlag()) { // put forward read first
  //            std::iter_swap(seg.begin(), seg.begin() + 1);
  //          }
    // Read level filtering step
    bool keep = true;
    for (const auto &s: seg) {
      //filtering reads before masking because the eof filtering will not update NM tag
      if ( cpputil::GetNM(s) - cpputil::CountNBasesInAlignment(s) > opt.max_mismatch_filter ||
          cpputil::GetNM(s) > opt.max_nm_filter) {
        if (opt.verbose) {
          std::cerr << "Discard read by edit distance filter " << s.Qname() <<"\n";
        }
        keep = false;
      }
    }
    if (not keep) {
      ++errorstat.discard_frag_counter;
      continue;
    }
    //EOF filter
    if (seg.size() == 2) {
      ++errorstat.pair_counter;
      cpputil::TrimPairFromFragEnd(seg.front(), seg.back(), opt.fragend_dist_filter);
    } else if(seg.size() == 1) {
      ++errorstat.single_counter;
      cpputil::TrimSingleFromFragEnd(seg.front(), opt.fragend_dist_filter);
    }

    //for readlevel output
    int r1_q0_den = 0, r2_q0_den = 0;
    int r1_q30_den =0, r2_q30_den = 0;
    int r1_q0_nerror = 0, r2_q0_nerror = 0;
    int r1_q30_nerror = 0, r2_q30_nerror = 0;
    string chrname = seg.front().ChrName(bamheader);
    for (int qcut : errorstat.cutoffs) {
      auto res = NumEffBases(seg, gr, ref, chrname, blacklist, errorstat, qcut, false, opt.count_nonoverlap);
      errorstat.qcut_neval[qcut].first += res.first;
      errorstat.qcut_neval[qcut].second += res.second;
      if (!opt.read_level_stat.empty()) {
        if (qcut == 0) {
          r1_q0_den = res.first;
          r2_q0_den = res.second;
        } else if (qcut == 30) {
          r1_q30_den = res.first;
          r2_q30_den = res.second;
        }
      }
    }

    int num = 0;
    // get denominator
    auto den = NumEffBases(seg, gr, ref, chrname, blacklist, errorstat, opt.bqual_min, true, opt.count_nonoverlap);
    errorstat.neval += den.first + den.second;
    if (!opt.cycle_level_stat.empty()) {
      CycleBaseCount(seg, gr, errorstat);
    }

    // first pass gets all variants in ROI
    std::map<cpputil::Variant, vector<cpputil::Variant>> var_vars;
    int rstart = 0;
    int rend = std::numeric_limits<int>::max();
    if (opt.count_nonoverlap) {
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

    for (auto &s : seg) {
      auto vars = cpputil::GetVar(s, bamheader, ref);
      for (auto &var : vars) {
        if (var.contig_start >= rstart && var.contig_start < rend)
          var_vars[var].push_back(var);
      }
    }

    // Currenlty, always assume the variants are concordance between R1 and R2
    //TODO: if count_nonoverlap == true, allow single variants in the overhang regions.
    std::vector<std::vector<cpputil::Variant>> quality_masked;
    quality_masked.reserve(var_vars.size());
    //Consolidate SNPs which has too purposes:
    // 1. change the variant quality of both reads to the lowest one so that they can be filtered of keep as a whole
    // 2. break doublets and etc. if not all bases are satisfying baseq cutoff.
    for (const auto& it: var_vars) {
      if (it.first.isIndel()) {
        quality_masked.push_back(it.second);
        continue;
      }
      if(it.second.size() == 1) {
        quality_masked.push_back(it.second);
      } else {
        auto vars = cpputil::snppair_consolidate(it.second[0], it.second[1], opt.bqual_min);
        quality_masked.insert(quality_masked.end(), vars.begin(), vars.end());
      }
    }
    for (const auto& readpair: quality_masked) {
      for (const auto& var : readpair) {
        vector<bool> real_muts;
        bool found;
        if (var.Type() == "SNV") {
          real_muts.resize(var.alt_seq.size(), false);
        }
        if (bcf_reader.IsOpen()) {
           found = cpputil::search_var_in_database(bcf_reader, var, known, "PRI_VCF", real_muts, true, opt.bqual_min, opt.all_mutant_frags);
        }
        if (bcf_reader2.IsOpen() && not found) {
          int maskbyvcf1 = n_true_mut(real_muts);
          found = cpputil::search_var_in_database(bcf_reader2, var, known, "SEC_VCF", real_muts, true, opt.bqual_min, opt.all_mutant_frags);
          if (var.Type() == "SNV" and var.var_qual >= opt.bqual_min)
            errorstat.nerror_masked_by_vcf2 += n_true_mut(real_muts) - maskbyvcf1;
        }
        if (mafr.IsOpen() && not found) {
          found = cpputil::search_var_in_database(mafr, var, known, "MAF", real_muts, true, opt.bqual_min, opt.all_mutant_frags);
        }
        if (not found) { // error site
          int nerr = n_false_mut(real_muts);
          if (var.Type() == "SNV" && var.var_qual >= opt.bqual_min && readpair.size() == 2) {
            num += nerr;
          }
          // detail error profiles
          if (n_true_mut(real_muts) > 0) { // partially match
            auto avars = cpputil::var_atomize(var);
            for (unsigned ai = 0; ai < real_muts.size(); ++ai) {
              if (not real_muts[ai]) {
                if (opt.all_mutant_frags) ferr << avars[ai] << '\n';
                else {
                  if (avars[ai].var_qual >= opt.bqual_min && readpair.size() == 2)
                    ferr << avars[ai] << '\n';
                }
              }
            }
          } else {
            if (opt.all_mutant_frags) ferr << var << '\n';
            else {
              if ((var.isIndel() || var.var_qual >= opt.bqual_min) && readpair.size() == 2)
                ferr << var << '\n';
            }
          }

          // other error profiles
          if (var.Type() == "SNV") {
            for (int qcut : errorstat.cutoffs) {
              if (var.var_qual >= qcut) {
                if (var.first_of_pair) {
                  errorstat.qcut_nerrors[qcut].first += nerr;
                } else {
                  errorstat.qcut_nerrors[qcut].second += nerr;
                }
                if (!opt.read_level_stat.empty()) {
                  if (qcut == 0) {
                    if (var.first_of_pair ) r1_q0_nerror += nerr;
                    else r2_q0_nerror += nerr;
                  }
                  if (qcut == 30) {
                    if (var.first_of_pair ) r1_q30_nerror += nerr;
                    else r2_q30_nerror += nerr;
                  }
                }
              }
            }
            //per cycle profile
            if (!opt.cycle_level_stat.empty()) {
              int tpos = var.alt_start < 0 ?  abs(var.alt_start + 1): var.alt_start;
              var.first_of_pair ? errorstat.R1_q0_error[tpos] += nerr : errorstat.R2_q0_error[tpos] += nerr;
              if (var.var_qual >= 30) {
                var.first_of_pair ? errorstat.R1_q30_error[tpos] += nerr: errorstat.R2_q30_error[tpos] += nerr;
              }
            }
          } //end other profiles
        }
      }
    }
    errorstat.nerror += num;

    if (readlevel.is_open()) {
      readlevel << seg.front().Qname() << '\t' << r1_q0_nerror << '\t' << r1_q0_den  << '\t'
                << r1_q30_nerror << '\t' << r1_q30_den << '\t'
                << r2_q0_nerror << '\t' << r2_q0_den << '\t'
                << r2_q30_nerror << '\t' << r2_q30_den << '\t'
                << abs(seg.front().InsertSize()) << '\n';
    }
  } // end for
  return;
}

int codec_accuracy(int argc, char ** argv) {
  AccuOptions opt;
  int parse_ret =  accuracy_parse_options(argc, argv, opt);
  if (parse_ret) return 1;
  if (argc == 1) {
    accuracy_print_help();
    return 1;
  }
  //if (!opt.vcf.empty()) {
  cpputil::BCFReader bcf_reader;
  cpputil::BCFReader bcf_reader2;
  if (!opt.vcfs.empty()) {
//    for (const auto& a : opt.vcfs) {
//      std::cerr << a << std::endl;
//    }
    bcf_reader.Open(opt.vcfs[0].c_str(), opt.sample);
    if (opt.vcfs.size() == 2) {
      bcf_reader2.Open(opt.vcfs[1].c_str());
    }
    else if (opt.vcfs.size() > 2){
      std::cerr << "Cannot read more than 2 vcfs currently\n";
      return 1;
    }
  }
  SeqLib::RefGenome ref;
  ref.LoadIndex(opt.reference);
  std::ofstream stat(opt.accuracy_stat);
  std::ofstream ferr(opt.error_prof_out);
  std::ofstream known(opt.known_var_out);
  std::ofstream context(opt.context_count);
  std::ofstream readlevel;
  std::ofstream cyclelevel;
  string error_profile_header =
      "chrom\tref_pos\tref\talt\ttype\tread_pos\tsnv_base_qual\tfirst_of_pair\tread_name";
  ferr << error_profile_header << std::endl;

  string known_var_header =
      "chrom\tref_pos\tref\talt\ttype\tread_pos\tsnv_base_qual\tfirst_of_pair\tread_name\tevidence";
  known << known_var_header << std::endl;
  if (not opt.read_level_stat.empty()) {
    readlevel.open(opt.read_level_stat);
    string read_level_header =
        "read_name\tR1_q0_nerror\tR1_q0_efflen\tR1_q30_nerror\tR1_q30_efflen\tR2_q0_nerror\tR2_q0_efflen\tR2_q30_nerror\tR2_q30_efflen\tfraglen";
    readlevel << read_level_header << std::endl;
  }
  if (not opt.cycle_level_stat.empty()) {
    cyclelevel.open(opt.cycle_level_stat);
    cyclelevel << "cycle\tR1_q0_error\tR1_q0_cov\tR1_q30_error\tR1_q30_cov\t"
            "R2_q0_error\tR2_q0_cov\tR2_q30_error\tR2_q30_cov\t"
            "R1_q0_erate\tR1_q30_erate\tR2_q0_erate\tR2_q30_erate\n";
  }
  cpputil::InsertSeqFactory isf(opt.bam, opt.mapq, opt.load_supplementary, opt.load_secondary, opt.load_duplicate, false);
//  SeqLib::BamWriter trim_bam_writer;
//  if (!opt.trim_bam.empty()) {
//    trim_bam_writer.SetHeader(isf.bamheader());
//    trim_bam_writer.Open(opt.trim_bam);
//    trim_bam_writer.WriteHeader();
//  }
  cpputil::MAFReader mafr;
  if (!opt.maf_file.empty()) {
    mafr.Open(opt.maf_file);
  }
  const int L = 250;
  auto errorstat = cpputil::ErrorStat(L);
  if (opt.detail_qscore_prof) {
    errorstat = cpputil::ErrorStat({0,10,20,30}, 250);
  }

  cpputil::TargetLayout tl(isf.bamheader(), opt.bed_file);
  // SeqLib::GenomicRegion gr;
  //while(tl.NextRegion(gr)) {
  vector<vector<cpputil::Segments>> chunk;
  for (int i = 0; i < tl.NumRegion(); ++i) {
    const auto& gr = tl[i];
    if (i % 100000 == 0) {
      auto timenow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
      std::cerr << i + 1 << " region processed. Last position: " << gr << std::ctime(&timenow) << std::endl;
    }
    if (isf.ReadByRegion(cpputil::ArePEAlignmentOverlapAtLeastK, gr, chunk, opt.pair_min_overlap, "", opt.load_unpair)) {
      std::set<int32_t> blacklist;
      if (bcf_reader.IsOpen()) {
        bcf_reader.vcf_to_blacklist_snv(gr, isf.bamheader(), blacklist);
      }
      if (bcf_reader2.IsOpen()) {
        bcf_reader.vcf_to_blacklist_snv(gr, isf.bamheader(), blacklist);
      }
      for (auto frag : chunk) {
        ErrorRateDriver(frag, isf.bamheader(), ref, &gr, blacklist, opt,
                        bcf_reader, bcf_reader2,
                        mafr, ferr, known, readlevel, errorstat);
      } // end for
    } //end if
  } //end for

  std::cerr << "All region processed \n";
  vector<string> header = {"qcutoff",
                                     "n_bases_eval",
                                     "n_A_eval",
                                     "n_C_eval",
                                     "n_G_eval",
                                     "n_T_eval",
                                     "n_errors",
                                     "erate",
                                     "n_errors_masked_by_vcf2",
                                     "num_pairs",
                                     "num_unpaireds",
                                     "n_discarded"};

  for (const auto& duo : errorstat.qcut_nerrors) {
    header.push_back("q" + std::to_string(duo.first) + "_n_bases_eval" );
    header.push_back("q" + std::to_string(duo.first) + "_n_errors" );
    header.push_back("q" + std::to_string(duo.first) + "_erate" );
    header.push_back("q" + std::to_string(duo.first) + "_R1_n_bases_eval" );
    header.push_back("q" + std::to_string(duo.first) + "_R1_n_errors" );
    header.push_back("q" + std::to_string(duo.first) + "_R1_erate" );
    header.push_back("q" + std::to_string(duo.first) + "_R2_n_bases_eval" );
    header.push_back("q" + std::to_string(duo.first) + "_R2_n_errors" );
    header.push_back("q" + std::to_string(duo.first) + "_R2_erate" );
  }

  stat << cpputil::join(header, "\t") << std::endl;
  stat << opt.bqual_min << '\t'
       << errorstat.neval << '\t'
      << errorstat.base_counter['A'] << '\t'
      << errorstat.base_counter['C'] << '\t'
      << errorstat.base_counter['G'] << '\t'
      << errorstat.base_counter['T'] << '\t'
      << errorstat.nerror << '\t'
       << (float) errorstat.nerror / errorstat.neval << '\t'
       << errorstat.nerror_masked_by_vcf2 << '\t'
       << errorstat.pair_counter << '\t'
       << errorstat.single_counter << '\t'
       << errorstat.discard_frag_counter << '\t';

  int ii = 0;
  for (const auto& qcut : errorstat.cutoffs) {
    ++ii;
    int64_t r1_den = errorstat.qcut_neval[qcut].first;
    int64_t r2_den = errorstat.qcut_neval[qcut].second;
    int64_t r1_num = errorstat.qcut_nerrors[qcut].first;
    int64_t r2_num = errorstat.qcut_nerrors[qcut].second;
    stat << r1_den + r2_den  << '\t';
    stat << r1_num + r2_num << '\t';
    stat << (float) (r1_num + r2_num) / (r1_den + r2_den) << '\t';
    stat << r1_den << '\t';
    stat << r1_num << '\t';
    stat << (float) r1_num / r1_den<< '\t';
    stat << r2_den << '\t';
    stat << r2_num << '\t';
    stat << (float) r2_num / r2_den;
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
