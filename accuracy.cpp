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

#define OPT_QSCORE_PROF   261
#define OPT_READ_LEVEL_STAT   262
using std::string;
using std::vector;
int MYINT_MAX = std::numeric_limits<int>::max();

struct Options {
  string bam;
  string accuracy_stat = "accuracy_stat.txt";
  string error_prof_out = "error_prof_out.txt";
  int mapq = 10;
  int bqual_min = 0;
  bool load_supplementary = false;
  bool load_secondary = false;
  bool load_unpair = false;
  int max_nm_filter = MYINT_MAX;
  int max_mismatch_filter = MYINT_MAX;
  int verbose = 0;
  vector<string> vcfs;
  string sample = "";
  string reference;
  int fragend_dist_filter = 0;
  string known_var_out = "known_var_out.txt";
  //string qscore_prof;
  bool detail_qscore_prof = false;
  string read_level_stat;
  //double min_frac_bqual_pass = 0.1;
  int pair_min_overlap = 0;
  bool overlap_only = false;
  string maf_file;
  string bed_file;
  string trim_bam = "trimmed.bam";
};


static struct option  long_options[] = {
    {"bam",                      required_argument,      0,        'b'},
    {"bed",                      required_argument,      0,        'L'},
    {"accuracy_stat",            required_argument,      0,        'a'},
    {"load_unpair",              no_argument,            0,        'u'},
    {"load_supplementary",       no_argument,            0,        'S'},
    {"load_secondary",           no_argument,            0,        '2'},
    {"mapq",                     required_argument ,     0,        'm'},
    {"vcfs",                     required_argument ,     0,        'V'},
    {"maf",                      required_argument ,     0,        'M'},
    {"sample",                   required_argument,      0,        's'},
    {"reference",                required_argument,      0,        'r'},
    {"error_prof_out",           required_argument,      0,        'e'},
    {"bqual_min",                required_argument,      0,        'q'},
    //{"min_frac_bqual_pass",      required_argument,      0,        'F'},
    {"known_var_out",            required_argument,      0,        'k'},
    {"pair_min_overlap",         required_argument,      0,        'p'},
    {"overlap_only",             no_argument,            0,        'O'},
    {"fragend_dist_filter",      required_argument,      0,        'd'},
    {"max_nm_filter",            required_argument,      0,        'n'},
    {"max_mismatch_filter",      required_argument,      0,        'x'},
    {"verbose",                  required_argument,      0,        'v'},
    //{"qscore_prof",              required_argument,      0,        OPT_QSCORE_PROF},
    {"detail_qscore_prof",       no_argument,            0,        OPT_QSCORE_PROF},
    {"read_level_stat",          required_argument,      0,        OPT_READ_LEVEL_STAT},
    {0,0,0,0}
};

const char* accuracy_short_options = "b:a:m:v:S2us:r:e:q:k:OM:p:d:n:x:V:L:";

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
  std::cerr<< "-V/--vcfs,                             comma separated VCF file(s) for germline variants or whitelist variants[null].\n";
  std::cerr<< "-M/--maf,                              MAF file for somatic variants [null].\n";
  std::cerr<< "-s/--sample,                           sample from the VCF file [null].\n";
  std::cerr<< "-r/--reference,                        reference sequence in fasta format [null].\n";
  std::cerr<< "-a/--accuracy_stat,                    output reporting accuracy for each alignment [accuracy_stat.txt].\n";
  std::cerr<< "-e/--error_prof_out,                   Error profile output in plain txt format [error_prof_out.txt].\n";
  std::cerr<< "-k/--known_var_out,                    Output for known var. [known_var_out.txt].\n";
  //std::cerr<< "--qscore_prof,                         Output qscore prof. First column is qscore cutoff; second column is number of bases in denominator\n";
  std::cerr<< "--detail_qscore_prof,                  Output finer scale qscore cutoffs, error rates profile. The default is only q0, q30 [false]. \n";
  std::cerr<< "--read_level_stat,                     Output read level stat.\n";
  std::cerr<< "\nFiltering Options:\n";
  std::cerr<< "-q/--bqual_min,                        Skip bases with baseQ smaller than this when calculating error rate [0].\n";
  std::cerr<< "-n/--max_edit_filter,                  Skip a read if its NM tag is larger than this value [INT_MAX].\n";
  std::cerr<< "-x/--max_nonNedit_filter,              Skip a read if the number of non-N bases edits is larger than this value [INT_MAX].\n";
  std::cerr<< "-d/--fragend_dist_filter,              Consider a variant if its distance to the fragment end is at least this value [0].\n";
  //std::cerr<< "-F/--min_frac_bqual_pass,              If bqual_min is none 0. The number of bases passing bqual filter has to be larger than or equal to this frac [0.1].\n";
  std::cerr<< "-p/--pair_min_overlap,                 When using selector, the minimum overlap between the two ends of the pair. -1 for complete overlap, 0 no overlap required [0].\n";
  std::cerr<< "-O/--overlap_only,                     Count only overlapped region of a read pair. Default [false].\n";
}

int accuracy_parse_options(int argc, char* argv[], Options& opt) {
  int option_index;
  int next_option = 0;
  do {
    next_option = getopt_long(argc, argv, accuracy_short_options, long_options, &option_index);
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
      case 'O':
        opt.overlap_only = true;
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
      case OPT_QSCORE_PROF:
        opt.detail_qscore_prof = true;
        break;
      case OPT_READ_LEVEL_STAT:
        opt.read_level_stat = optarg;
        break;
      case 'p':
        opt.pair_min_overlap = atoi(optarg);
        break;
      default:accuracy_print_help();
        return 1;
    }
  } while (next_option != -1);

  return 0;
}

struct ErrorStat {
  int64_t neval = 0;
  int64_t nerror = 0;
  int64_t nerror_masked_by_vcf2 = 0;
  int64_t discard_frag_counter = 0;
  int64_t pair_counter = 0;
  int64_t single_counter = 0;
  vector<int> cutoffs;
  vector<int64_t> R1_q0_cov;
  vector<int64_t> R1_q30_cov;
  vector<int64_t> R1_q0_error;
  vector<int64_t> R1_q30_error;
  vector<int64_t> R2_q0_cov;
  vector<int64_t> R2_q30_cov;
  vector<int64_t> R2_q0_error;
  vector<int64_t> R2_q30_error;
  vector<int64_t> base_counter;
  //qscore cutoff -> read1, read2
  std::map<int, std::pair<int64_t, int64_t>> qcut_neval;
  std::map<int, std::pair<int64_t, int64_t>> qcut_nerrors;

  //ErrorStat() = delete;
  ErrorStat(int L) {
    cutoffs = {0, 30};
    qcut_neval[0] = std::make_pair(0,0);
    qcut_neval[0] = std::make_pair(0,0);
    qcut_nerrors[30] = std::make_pair(0,0);
    qcut_nerrors[30] = std::make_pair(0,0);
    Init(L);
  }
  ErrorStat(const vector<int>& qcuts, int L): cutoffs(qcuts){
    assert(std::find(qcuts.begin(), qcuts.end(), 0) != qcuts.end());
    assert(std::find(qcuts.begin(), qcuts.end(), 30) != qcuts.end());
    for (auto q : qcuts) {
      qcut_neval[q] = std::make_pair(0,0);
      qcut_nerrors[q] = std::make_pair(0,0);
    }
    Init(L);
  }

 private:
  void Init(int L) {
    R1_q0_cov.resize(L);
    R1_q30_cov.resize(L);
    R1_q0_error.resize(L);
    R1_q30_error.resize(L);
    R2_q0_cov.resize(L);
    R2_q30_cov.resize(L);
    R2_q0_error.resize(L);
    R2_q30_error.resize(L);
    base_counter.resize(4);
  }
};

std::vector<int> NumEffBases(const cpputil::Segments& seg,
    const SeqLib::GenomicRegion* const gr,
    int minbq,
    bool count_overhang){
  //return a vector of size 8. A,C,G,T count, first read1, then read2
  std::vector<int> r1;
  std::vector<int> r2;
  if (seg.size() == 1) {
    std::pair<int,int> range;
    if (gr) {
      range = cpputil::GetBamOverlapQStartAndQStop(seg.front(), *gr);
    } else {
      range.first = seg.front().AlignmentPosition();
      range.second = seg.front().AlignmentEndPosition();
    }
    if (seg.front().FirstFlag()) {
      r1 = cpputil::CountValidBaseInMatchedBases(seg.front(), minbq, range.first, range.second);
    }
    else {
      r2 = cpputil::CountValidBaseInMatchedBases(seg.front(), minbq, range.first, range.second);
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
      r1 = cpputil::CountValidBaseInMatchedBases(seg.front(), minbq, range_front.first, range_front.second);
      r2 = cpputil::CountValidBaseInMatchedBases(seg.back(), minbq, range_back.first, range_back.second);
    } else {
      r2 = cpputil::CountValidBaseInMatchedBases(seg.front(), minbq, range_front.first, range_front.second);
      r1 = cpputil::CountValidBaseInMatchedBases(seg.back(),  minbq, range_back.first, range_back.second);
    }
  }
  r1.insert(r1.end(), r2.begin(), r2.end());
  return r1;
}

void CycleBaseCount(const cpputil::Segments& seg,
                    const SeqLib::GenomicRegion* const gr,
                    ErrorStat& es) {
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

void ErrorRateDriver(const vector<cpputil::Segments>& frag,
                    const SeqLib::BamHeader& bamheader,
                    const SeqLib::RefGenome& ref,
                    const SeqLib::GenomicRegion* const gr,
                    const Options& opt,
                    const cpputil::BCFReader& bcf_reader,
                    const cpputil::BCFReader& bcf_reader2,
                    const cpputil::MAFReader& mafr,
                    std::ofstream& ferr,
                    std::ofstream& known,
                    std::ofstream& readlevel,
                    ErrorStat& errorstat) {
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


    int r1_q0_den = 0, r2_q0_den = 0;
    int r1_q30_den =0, r2_q30_den = 0;
    for (int qcut : errorstat.cutoffs) {
      auto res = NumEffBases(seg, gr, qcut, not opt.overlap_only);
      errorstat.qcut_neval[qcut].first += std::accumulate(res.begin(), res.begin() + 4, 0);
      errorstat.qcut_neval[qcut].second += std::accumulate(res.begin() + 4, res.end(), 0);
      if (qcut == 0) {
        r1_q0_den = std::accumulate(res.begin(), res.begin() + 4, 0);
        r2_q0_den = std::accumulate(res.begin() + 4, res.end(), 0);
      } else if (qcut == 30) {
        r1_q30_den = std::accumulate(res.begin(), res.begin() + 4, 0);
        r2_q30_den = std::accumulate(res.begin() + 4, res.end(), 0);
      }
    }

    int num = 0;
    int r1_q0_nerror = 0, r2_q0_nerror = 0;
    int r1_q30_nerror = 0, r2_q30_nerror = 0;
    // get denominator
    auto den = NumEffBases(seg, gr, opt.bqual_min, not opt.overlap_only);
    for (unsigned bi = 0; bi < 4; ++bi) {
      errorstat.base_counter[bi] += (den[bi] + den[bi + 4]);
    }
    errorstat.neval += std::accumulate(den.begin(), den.end(), 0);
    CycleBaseCount(seg, gr, errorstat);

    // first pass gets all variants in ROI
    std::map<cpputil::Variant, vector<cpputil::Variant>> var_vars;
    int rstart = 0;
    int rend = std::numeric_limits<int>::max();
    if (not opt.overlap_only) {
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
    // second pass classifies the variants
    for (const auto& duo: var_vars) {
      for (const auto& var : duo.second) {
        if (bcf_reader.IsOpen()
            && bcf_reader.var_exist(var.contig, var.contig_start, var.alt_seq.c_str())) { // known var
          known << var << '\t' << "PRI_VCF" << '\n';
        }
        else if (bcf_reader2.IsOpen()
            && bcf_reader2.var_exist(var.contig, var.contig_start, var.alt_seq.c_str())) { // second known var
          known << var << '\t' << "SEC_VCF" << '\n';
          if (var.Type() == "SNV") {
            errorstat.nerror_masked_by_vcf2 += 1;
          }
        }
        else if (mafr.IsOpen()
            && mafr.VarExist(var.contig, var.contig_start + 1, var.alt_seq)) { //known maf site
          known << var << '\t' << "MAF" << '\n';
        }
        else { // error site
          if (var.Type() == "SNV" && var.var_qual >= opt.bqual_min) {
            num += 1;
            ferr << var << '\t' << (duo.second.size() == 2 ? 1 : 0) << '\n'; // concordant or not
          }
          if (var.Type() == "SNV") {
            for (int qcut : errorstat.cutoffs) {
              if (var.var_qual >= qcut) {
                var.first_of_pair ? ++errorstat.qcut_nerrors[qcut].first: ++errorstat.qcut_nerrors[qcut].second;
                if (qcut == 0) {var.first_of_pair ? ++r1_q0_nerror : ++r2_q0_nerror;}
                if (qcut == 30) {var.first_of_pair ? ++r1_q30_nerror : ++r2_q30_nerror;}
              }
            }
            int tpos = var.alt_start < 0 ?  abs(var.alt_start + 1): var.alt_start;
            var.first_of_pair ? errorstat.R1_q0_error[tpos] += 1 : errorstat.R2_q0_error[tpos] += 1;
            if (var.var_qual >= 30) {
              var.first_of_pair ? errorstat.R1_q30_error[tpos] += 1: errorstat.R2_q30_error[tpos] += 1;
            }
          } else { // indel
            ferr << var << '\t' << (duo.second.size() == 2 ? 1 : 0) << '\n'; // concordant or not
          }
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
  Options opt;
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
  std::ofstream readlevel;
  string error_profile_header =
      "chrom\tref_pos\tref\talt\ttype\tread_pos\tsnv_base_qual\tfirst_of_pair\tread_name\tconcordant";
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
  cpputil::InsertSeqFactory isf(opt.bam, opt.mapq, opt.load_supplementary, opt.load_secondary, false);
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
  auto errorstat = ErrorStat(L);
  if (opt.detail_qscore_prof) {
    errorstat = ErrorStat({0,10,20,30}, 250);
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
    if (isf.ReadByRegion(cpputil::ArePEAlignmentOverlapAtLeastK, gr, chunk, opt.pair_min_overlap, opt.load_unpair)) {
      for (auto frag : chunk) {
        ErrorRateDriver(frag, isf.bamheader(), ref, &gr, opt,
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
      << errorstat.base_counter[0] << '\t'
      << errorstat.base_counter[1] << '\t'
      << errorstat.base_counter[2] << '\t'
      << errorstat.base_counter[3] << '\t'
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

  stat << "cycle\tR1_q0_error\tR1_q0_cov\tR1_q30_error\tR1_q30_cov\t"
          "R2_q0_error\tR2_q0_cov\tR2_q30_error\tR2_q30_cov\t"
          "R1_q0_erate\tR1_q30_erate\tR2_q0_erate\tR2_q30_erate\n";
  for (int i = 0; i < L; ++i) {
    stat << i << '\t';
    stat << errorstat.R1_q0_error[i] << '\t' << errorstat.R1_q0_cov[i] << '\t';
    stat << errorstat.R1_q30_error[i] << '\t' << errorstat.R1_q30_cov[i] << '\t';
    stat << errorstat.R2_q0_error[i] << '\t' << errorstat.R2_q0_cov[i] << '\t';
    stat << errorstat.R2_q30_error[i] << '\t' << errorstat.R2_q30_cov[i] << '\t';
    stat << (float) errorstat.R1_q0_error[i] / errorstat.R1_q0_cov[i] << '\t';
    stat << (float) errorstat.R1_q30_error[i] / errorstat.R1_q30_cov[i] << '\t';
    stat << (float) errorstat.R2_q0_error[i] / errorstat.R2_q0_cov[i] << '\t';
    stat << (float) errorstat.R2_q30_error[i] / errorstat.R2_q30_cov[i] << '\n';
  }
  return 0;
}
