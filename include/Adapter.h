//
// Created by Ruolin Liu on 8/17/20.
//

#ifndef ADAPTERTRIM_INCLUDE_ADAPTER_H_
#define ADAPTERTRIM_INCLUDE_ADAPTER_H_
#include <string>
#include <map>
#include "FastxIO.h"
#include "DNAUtils.h"
#include <deque>
namespace CDS {
using FxRecord = cpputil::ExtFastxRecord;


const std::string O2N_LINKER = "AGATCAGTCGTTCACCGACTGCCCACAAGCTGCTATTGAGAGTAAGCACGTACGACTGATCT"; // first 20bp
/*
 * Legacy design
 */
std::map<const std::string, const std::string> LINKERS = {
    {"medium_v1", "AGATCGGAAGAGCTTCATCATTAGATCCATTACACATCAGATTAGTACCAGCTTCGAGGATCAACACGTCAGAGTCTAGCTGGTGATAGGAAGTGTAGAATGTGTAACTGACTTAACGCTCTTCCGATCT"},
    {"medium_v2", "AGATCGGAAGAGCTTCATCATTAGATCCATTACACATCAGATTAGTACCAGCTTCGAGGATCAACACGTCAGAGTCTAGCTGGTGATAGGAAGTGTAGAATGTGTAACTGACTTATGGCTCTTCCGATCT"},
    {"long_v1", "AGATCGGAAGAGCTTCATCATTAGATCCATTAATGTTACACTTCAACTCTTCACCCACATCAGATTAGTACCAGCTTCGAGGATCAACACGTCAGAGTCTAGCTGGTGATAGGAAGTGTAGGTAACATAGACGAAGTTATCAACAATGTGTAACTGACTTAACGCTCTTCCGATCT"},
    {"o2n", O2N_LINKER}
};

const std::string ILLUMINA_PAIR1_ADAPTER_SUFFIX = "ACACGTCTGAACTCCAGTCA"; // unique region R1 after 13bp
const std::string ILLUMINA_PAIR2_ADAPTER_SUFFIX = "GTCGTGTAGGGAAAGAGTGT"; // unique region R2 after 13bp
const std::string COMMON_13BP="AGATCGGAAGAGC";
const std::string ILLUMIA_PAIR1_BARCODE_ADAPTER_FULL = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC88888888ATCTCGTATGCCGTCTTCTGCTTG";
const std::string ILLUMIA_PAIR2_BARCODE_ADAPTER_FULL = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT66666666GTGTAGATCTCGGTGGTCGCCGTATCATT";
constexpr int TRUESEQ_PREFIX_LEN = 13;
const std::string R1_PREFIX_CUSTOM = "GCTCTTCCGATCTGACCT";
const std::string R1_SUFFIX_CUSTOM = "ACGCTAGATCGGAAGAGC";
const std::string R2_PREFIX_CUSTOM = "GCTCTTCCGATCTAGCGT";
const std::string R2_SUFFIX_CUSTOM = "AGGTCAGATCGGAAGAGC";


int MY_INT_LOWEST = std::numeric_limits<int>::lowest();
typedef enum {TRIM_INSUF = 0, TRIM, TRIM13, SINGLE, DOUBLE_LIGATION, TRIMEND, NOTRIM} TrimMode;
typedef enum {AdapFwd = 0, AdapRev} AdapDirct;
/*
 * In the output bam the trim_insuf contains both untrimmed and trimend
 */
inline std::string ToString(TrimMode tm) {
  std::vector<std::string> to_string = {"trim_insuf", "trim", "trim13", "single", "double_ligation", "trimend", "notrim"};
  return to_string[tm];
}

struct TrimStatus {
  TrimMode tm= NOTRIM;
  int adap_end = -1;
  AdapDirct adap_dirct = AdapFwd;
};
typedef std::pair<TrimStatus, TrimStatus> FragTrimS;


struct DP_Matrix {
  struct MatchBoard {
    // i for row index, j for col index
    // max_i_or_j is the end position of the best match between two strings
    // start_i_or_j is the start position of the best match between two strings
    // max_score is the score of max_i, max_j in the DP matrix
    unsigned max_i;
    unsigned max_j;
    unsigned start_i;
    unsigned start_j;
    int max_score;
    int last_row_max_score;
    unsigned last_row_max_j;
    int nm;
    int nmatch;
    std::string refgap;
    std::string querygap;
    std::string symbols;
    int global_max_score = std::numeric_limits<int>::lowest();
    MatchBoard() : max_i(0), max_j(0), start_i(0), start_j(0), max_score(MY_INT_LOWEST),
      last_row_max_score(MY_INT_LOWEST), last_row_max_j(0), nm(0), nmatch(0) {}
    void Print(bool q_no3_clip = false) {
      if (q_no3_clip) {
        std::cerr << "score, nmatch, nm: " << last_row_max_score << ", " << nmatch << ", " << nm << std::endl;
        std::cerr << "start_i: " << start_i << std::endl;
        std::cerr << "start_j, max_j: " << start_j << ", " << last_row_max_j << std::endl;
      } else {
        std::cerr << "score, nmatch, nm: " << max_score << ", " << nmatch << ", " << nm << std::endl;
        std::cerr << "start_i, max_i: " << start_i << ", " << max_i << std::endl;
        std::cerr << "start_j, max_j: " << start_j << ", " << max_j << std::endl;
      }
      std::cerr << refgap << "\n";
      std::cerr << symbols << "\n";
      std::cerr << querygap << "\n";
      std::cerr << '\n';
    }

    float PI() const {
      return (float) nmatch / (nm + nmatch);
    }
  };

  struct Element {
    int16_t dir; // 0 diag, -1 horizontal, 1 vertical
    int16_t score;
    Element() : dir(-1), score(0) {}
    Element(int d, int s) : dir(d), score(s) {}
  };

  DP_Matrix(std::string query, std::string ref, int stop_at_ref_basen = -1) : query_(query), ref_(ref) {
    nrow_ = query.size() + 1;
    ncol_ = (stop_at_ref_basen == -1 ? ref.size() : stop_at_ref_basen) + 1;
  }

  void Reset() {
    matrix_ = std::vector<std::vector<Element>>(nrow_, std::vector<Element>(ncol_));
  }

  void Reset(std::string query, std::string ref) {
    query_ = query;
    ref_ = ref;
    nrow_ = query.size() + 1;
    ncol_ = ref.size() + 1;
    matrix_ = std::vector<std::vector<Element>>(nrow_, std::vector<Element>(ncol_));
  }

  MatchBoard SimpleSearch(const int matchs, const int mismatchs, const int gaps, const int score_cutoff) {
    // return immediately if score_cutoff is reached.
    Reset();
    MatchBoard res;
    for (int i = 1; i < nrow_; ++i) {
      matrix_[i][0].dir = 1;
      matrix_[i][0].score = i * gaps;
    }
    for (int j = 1; j < ncol_; ++j) {
      for (int i = 1; i < nrow_; ++i) {
        int best_score = matrix_[i - 1][j - 1].score + (query_[i - 1] == ref_[j - 1] ? matchs : mismatchs);
        int best_dir = 0;
        if (matrix_[i - 1][j].score + gaps > best_score) {
          best_score = matrix_[i - 1][j].score + gaps;
          best_dir = 1;
        }
        if (matrix_[i][j - 1].score + gaps > best_score) {
          best_score = matrix_[i][j - 1].score + gaps;
          best_dir = -1;
        }
        matrix_[i][j] = Element(best_dir, best_score);
        if (best_score >= score_cutoff) {
          res.max_i = i;
          res.max_j = j;
          res.max_score = best_score;
          return res;
        }
      }
    }
    return res;
  }

  MatchBoard Compute(const int matchs, const int mismatchs, const int gaps,
                      bool q_no_5clip, bool r_no_5clip) {
    Reset();
    MatchBoard res;
    // init first row and col
    for (int i = 1; i < nrow_; ++i) {
      if (q_no_5clip) {
        matrix_[i][0].score = i * gaps;
      }
      matrix_[i][0].dir = 1;
    }
    if (r_no_5clip) {
      for (int j = 1; j < ncol_; ++j) {
        matrix_[0][j].score = j *gaps;
      }
    }
    for (int j = 1; j < ncol_; ++j) {
      for (int i = 1; i < nrow_; ++i) {
        int best_score = matrix_[i - 1][j - 1].score + (query_[i - 1] == ref_[j - 1] ? matchs : mismatchs);
        int best_dir = 0;
        if (matrix_[i - 1][j].score + gaps > best_score) {
          best_score = matrix_[i - 1][j].score + gaps;
          best_dir = 1;
        }
        if (matrix_[i][j - 1].score + gaps > best_score) {
          best_score = matrix_[i][j - 1].score + gaps;
          best_dir = -1;
        }
        matrix_[i][j] = Element(best_dir, best_score);
        if (best_score > res.max_score) {
          res.max_i = i;
          res.max_j = j;
          res.max_score = best_score;
        }
      }
      if (matrix_[nrow_ - 1][j].score > res.last_row_max_score) {
        res.last_row_max_score = matrix_[nrow_ - 1][j].score;
        res.last_row_max_j = j;
      }
    }
    res.global_max_score = matrix_[nrow_-1][ncol_-1].score;
    return res;
  }

  void BackTrace(MatchBoard& res, bool q_no_3clip, bool q_no_5clip, bool r_no_5clip) {
    res.nm = 0;
    res.nmatch = 0;
    res.refgap.clear();
    res.querygap.clear();
    res.symbols.clear();

    unsigned ii = res.max_i;
    unsigned start_j = res.max_j;
    unsigned start_i = 0;

    if (q_no_3clip) {
      start_j = res.last_row_max_j;
      ii = nrow_ - 1;
    }
    //BackTrace
    // prepare for search, max_i, max_j are the end positions of the match
    while (ii > 0) {
      switch (matrix_[ii][start_j].dir) {
        case -1:start_j--;
          res.refgap += ref_[start_j];
          res.querygap += '-';
          res.symbols += ' ';
          ++res.nm;
          break;
        case 1:ii--;
          res.querygap += query_[ii];
          res.refgap += '-';
          res.symbols += ' ';
          ++res.nm;
          break;
        case 0:start_j--; ii--;
          res.querygap += query_[ii];
          res.refgap += ref_[start_j];
          if (query_[ii] != ref_[start_j]) {
            ++res.nm;
            res.symbols += ' ';
          }
          else {
            ++res.nmatch;
            res.symbols += '|';
          }
          break;
      }
      if (start_j == 0) {
        start_i = ii;
        break;
      }
    }
    if (q_no_5clip && start_i > 0) {
      res.nm += start_i;
      res.refgap += std::string(start_i, '-');
      res.symbols += std::string(start_i, ' ');
      auto x = query_.substr(0, start_i);
      std::reverse(x.begin(), x.end());
      res.querygap += x;
    }
    if (r_no_5clip && start_j > 0) {
      res.nm += start_j;
      res.querygap += std::string(start_j, '-');
      res.symbols += std::string(start_j, ' ');
      auto x = ref_.substr(0, start_j);
      std::reverse(x.begin(), x.end());
      res.refgap += x;
    }
    std::reverse(res.refgap.begin(), res.refgap.end());
    std::reverse(res.querygap.begin(), res.querygap.end());
    std::reverse(res.symbols.begin(), res.symbols.end());
    res.start_i = start_i;
    res.start_j = start_j;
  }

  static MatchBoard SemiGlobalAlign(const std::string& pattern, const std::string& ref, int match, int mismatch, int gap) {
    DP_Matrix full = DP_Matrix(pattern, ref); // aligned full linker
    auto full_mb = full.Compute(match, mismatch, gap, true, false);
    full.BackTrace(full_mb, true, true, false);
    return full_mb;
  }
  static MatchBoard LocalAlign(const std::string& pattern, const std::string& ref, int match, int mismatch, int gap) {
    DP_Matrix semi(pattern, ref);
    DP_Matrix::MatchBoard semi_mb = semi.Compute(match, mismatch, gap, false, false);
    semi.BackTrace(semi_mb, false, false, false);
    return semi_mb;
  }

  std::string query_;
  std::string ref_;
  int nrow_;
  int ncol_;
  std::vector<std::vector<Element>> matrix_;
};

class AdapterMatch {

  // Query (linker) as row in the DP, reference (read) as col in the DP
  //std::string adap_;
  std::string read_;
  unsigned rl_;
  int gap_s_;
  int match_s_;
  int mis_match_s_;
  unsigned three_end_adapter_min_len_;
  unsigned junction_adapter_min_len_;
  unsigned ligation_min_len_;
  int fast_search_score_cutoff_;

 public:
  AdapterMatch (std::string r, int matchs, int mismatchs, int gaps, unsigned min3, unsigned minj, unsigned lml, int fastscore):
      read_(r),
      rl_(r.size()),
      gap_s_(gaps), match_s_(matchs), mis_match_s_(mismatchs),
      three_end_adapter_min_len_(min3), junction_adapter_min_len_(minj), ligation_min_len_(lml),
      fast_search_score_cutoff_(fastscore)
  {}
  // V1 adapter
  AdapterMatch (std::string r, int matchs, int mismatchs, int gaps, unsigned min3, unsigned minj):
      AdapterMatch(r, matchs, mismatchs, gaps, min3, minj, 15, 20)
  {}
  // V2 adapter
  AdapterMatch (std::string r, int matchs, int mismatchs, int gaps) :
      AdapterMatch(r, matchs, mismatchs, gaps, 0, 0, 15, 16)
  {}

  unsigned TrimBeforeInsert(const std::string& pattern, TrimStatus& ts, int match_till, int debug = 0) const {
    const int MAX_ALLOWED_MISMATCH = 4;
    auto mb = DP_Matrix::SemiGlobalAlign(pattern, read_.substr(0, match_till), match_s_, mis_match_s_, gap_s_);
    if (debug == 1) {
      std::cerr << "Trim Before Insert\n";
      mb.Print(true);
    }
    if (mb.nm <= MAX_ALLOWED_MISMATCH) {
      ts.tm = TRIM;
      return mb.last_row_max_j;
    }
    else {
      ts.tm = NOTRIM;
      return 0;
    }
  }

  std::pair<unsigned, unsigned> TrimAfterInsert( const std::string& pattern, const int pl,
      TrimStatus& ts, const unsigned margin_after_first_trim, unsigned nbase_skip, bool allow_rc_search,
      int allow_diff = 2, int debug = 0) const { // heuristics
    //return adapter start, end
    // fast mapping full adpt
    DP_Matrix dp(pattern, read_);
    bool is_rc = false;
    std::string rc_pattern = pattern;
    DP_Matrix::MatchBoard mb = dp.SimpleSearch(match_s_, mis_match_s_, gap_s_, fast_search_score_cutoff_);
    if (allow_rc_search && mb.max_score == MY_INT_LOWEST) { // search reverse comp.
      cpputil::reverse_complement(rc_pattern);
      dp.Reset(rc_pattern, read_);
      mb = dp.SimpleSearch(match_s_, mis_match_s_, gap_s_, fast_search_score_cutoff_);
      if (mb.max_score != MY_INT_LOWEST) is_rc = true;
    }

    if (mb.max_score != MY_INT_LOWEST) {
      dp.BackTrace(mb, false, true, false); // for simple search
      if (debug == 1) {
        std::cerr << "Trim After Insert. First search\n";
        std::cerr << read_ << std::endl;
        std::cerr << pattern << std::endl;
        mb.Print();
      }
      if (mb.start_j < ligation_min_len_) { // Adapter dimer or double ligation
        DP_Matrix::MatchBoard full_mb;
        if (is_rc) {
          full_mb = DP_Matrix::SemiGlobalAlign(rc_pattern, read_, match_s_, mis_match_s_, gap_s_);
        } else {
          full_mb = DP_Matrix::SemiGlobalAlign(pattern, read_, match_s_, mis_match_s_, gap_s_);
        }
        if (full_mb.PI() > 0.80) {
          if (debug == 1) {
            full_mb.Print(true);
          }
          ts.tm = DOUBLE_LIGATION;
          if (is_rc) ts.adap_dirct = AdapRev;
          return std::make_pair(full_mb.start_j, full_mb.last_row_max_j);
        }
      }
      if (not is_rc) { //regular cds
        ts.tm = TRIM;
        return std::make_pair(mb.start_j, read_.size());
      }
    }
    // the fast search has not found the adapter
    // Now we try match the prefix of pattern of length pl.
    // Search at the end of the read first
    dp.Reset(pattern.substr(0, pl), read_.substr(nbase_skip));
    mb  = dp.Compute(match_s_, mis_match_s_, gap_s_,true, false);
    dp.BackTrace(mb, false, true, false);
    if (debug == 1) {
      std::cerr << "Trim After Insert\n";
      std::cerr << read_ << std::endl;
      std::cerr << pattern << std::endl;
      mb.Print();
    }
    //If adapter at the end of the reads
    if (mb.max_j >= rl_ - nbase_skip - margin_after_first_trim && mb.start_i == 0) {
      if (mb.max_j - mb.start_j >= three_end_adapter_min_len_) {
        ts.tm = TRIMEND;
        return std::make_pair(nbase_skip + mb.start_j, nbase_skip + mb.max_j);
      }
    }
    // Search else where in the read then
    dp.BackTrace(mb, true, true, false);
    if (debug == 1) {
      std::cerr << "Full adapter trim\n";
      mb.Print(true);
    }
    if (mb.nm <= allow_diff) {
      ts.tm = mb.last_row_max_j >= rl_ - nbase_skip - margin_after_first_trim ? TRIMEND : NOTRIM;
      return std::make_pair(nbase_skip + mb.start_j, nbase_skip + mb.last_row_max_j); // trim mode will be decided by next function
    }
    //Canot find adapter anywhere
    return std::make_pair(read_.size(), read_.size()); // notrim
  }

  unsigned TwoStepTrim(const std::string& linker, const std::string& trueseq_suffix, TrimStatus& ts, const int nbase_skip, int debug = 0) const {

    unsigned match_start, match_end;
    std::tie(match_start, match_end) = TrimAfterInsert(linker,
        TRUESEQ_PREFIX_LEN, ts, 4, nbase_skip, true, 1, debug);
    if (ts.tm == DOUBLE_LIGATION) { //
      return match_end;
    }
    ts.adap_end = match_end;
    if (ts.tm == TRIM) {
      return match_start;
    }
    if (match_start == read_.size()) { // notrim
      return match_start;
    }
    if (ts.tm == TRIMEND) {
      return match_start;
    }
    DecideSingleInsert(match_end, linker, trueseq_suffix, ts, debug);
    return match_start;
  }

  void DecideSingleInsert(const unsigned start_pos, const std::string& linker,
                          const std::string& illu_suffix, TrimStatus& ts, int debug = 0) const {

    unsigned illu_min_len = std::min(illu_suffix.size(), read_.size() - start_pos);

    //global alignment for matching illumina adapter
    DP_Matrix dpill(illu_suffix.substr(0, illu_min_len), read_.substr(start_pos, illu_min_len));
    DP_Matrix::MatchBoard illu_mb = dpill.Compute(match_s_, mis_match_s_, gap_s_, true, true);
    dpill.BackTrace(illu_mb, true, true, true);
    //local alignment for matching cds adapter, due to dualF
    std::string linker_suffix = linker.substr(TRUESEQ_PREFIX_LEN);
    unsigned cds_min_len = std::min(linker_suffix.size(), read_.size() - start_pos);
    auto cds_mb = DP_Matrix::LocalAlign(linker_suffix.substr(0, cds_min_len), read_.substr(start_pos, cds_min_len),
        match_s_, mis_match_s_, gap_s_);
    auto cds_mb_sgb = DP_Matrix::SemiGlobalAlign(linker_suffix.substr(0, cds_min_len), read_.substr(start_pos, cds_min_len),
                                            match_s_, mis_match_s_, gap_s_);

    if (debug == 1) {
      std::cerr << "trueseq match\n";
      std::cerr << read_.substr(start_pos, illu_min_len) << std::endl;
      std::cerr << illu_suffix.substr(0, illu_min_len) << std::endl;
      illu_mb.Print(true);
      std::cerr << "linker semiglobal match\n";
      std::cerr << read_.substr(start_pos, cds_min_len) << std::endl;
      std::cerr << linker_suffix.substr(0, cds_min_len) << std::endl;
      cds_mb_sgb.Print(true);
    }

    if ((cds_mb.max_j - cds_mb.start_j >= junction_adapter_min_len_ &&
        abs(cds_mb.start_j - cds_mb.start_i) <= junction_adapter_min_len_ * 0.1) // matched in similar distance to beginning
      //  || (float) cds_mb.nmatch / cds_min_len >= 0.8) {
      || cds_mb_sgb.PI() >= 0.8) {
      ts.tm = TRIM;
      return;
    }

    if(illu_mb.PI() >= 0.75) {
    //if(illu_mb.PI() >= 0.95) { // o2n
      ts.tm = SINGLE;
      return;
    }
    // try matching read with RC Adpt
    //TODO: Compute rc linker just once.
    std::string rc_linker = linker;
    cpputil::reverse_complement(rc_linker);
    auto rc_linker_suffix = rc_linker.substr(TRUESEQ_PREFIX_LEN);
    auto rc_cds_mb = DP_Matrix::LocalAlign(rc_linker_suffix.substr(0, cds_min_len), read_.substr(start_pos, cds_min_len),
        match_s_, mis_match_s_, gap_s_);
    auto rc_cds_mb_sgb = DP_Matrix::SemiGlobalAlign(rc_linker_suffix.substr(0, cds_min_len), read_.substr(start_pos, cds_min_len),
                                           match_s_, mis_match_s_, gap_s_);
    if (debug == 1) {
      std::cerr << "rc linker semiglobalmatch\n";
      rc_cds_mb_sgb.Print(true);
    }
    if ((rc_cds_mb.max_j - rc_cds_mb.start_j >= junction_adapter_min_len_ &&
        abs(rc_cds_mb.start_j - rc_cds_mb.start_i) <= junction_adapter_min_len_ * 0.1) // matched in similar distance to beginning
        //|| (float) rc_cds_mb.nmatch / cds_min_len >= 0.8) {
        || rc_cds_mb_sgb.PI() >= 0.8) {
      ts.tm = TRIM;
      ts.adap_dirct = AdapRev;
      return;
    }

    ts.tm = TRIM13;
    return;
  }

};

std::pair<std::string, std::string> QualTrim(std::string seq, std::string qual, int minbq= 3) {
  std::deque<unsigned> keeps;
  std::string cleaned_seq;
  std::string cleaned_qual;
  for(unsigned i = 0; i < qual.size(); ++i) {
    if (qual[i] >= minbq + 33) {
      keeps.push_back(i);
    }
  }
  for (unsigned i = 0; i < qual.size(); ++i) {
    if (i == keeps.front()) {
      keeps.pop_front();
      cleaned_seq += seq[i];
      cleaned_qual += qual[i];
    }
  }
  return std::make_pair(cleaned_seq, cleaned_qual);
}

void TrimUmi(const int r1_umi_len, const int r2_umi_len, FxRecord& read1, FxRecord& read2) {
  std::string umi;
  std::string umiq;
  if (r1_umi_len != 0) {
    std::string r1_umi = read1.seq.substr(0, r1_umi_len);
    std::string r1_umiq = read1.qual.substr(0, r1_umi_len);
    read1.seq = read1.seq.substr(r1_umi_len);
    read1.qual = read1.qual.substr(r1_umi_len);
    umi = r1_umi;
    umiq = r1_umiq;
  }
  if (r2_umi_len != 0 ) {
    std::string r2_umi = read2.seq.substr(0, r2_umi_len);
    std::string r2_umiq = read2.qual.substr(0, r2_umi_len);
    read2.seq = read2.seq.substr(r2_umi_len);
    read2.qual = read2.qual.substr(r2_umi_len);
    umi +=  "-" + r2_umi;
    umiq += "-" + r2_umiq;
  }
  read1.umi = cpputil::AnnotatedSeq(umi, umiq);
  read2.umi = cpputil::AnnotatedSeq(umi, umiq);
}

template<typename Options>
TrimStatus TrimFrontAndBack(const std::string& front_adap, const std::string& back_adap, const Options& opt, FxRecord& read,
      const int match_front_till) {
  //Return only two status, TRIM or NOTRIM based on the trimming of front_adap
  //The trimming status of front_adap is attached to the read name
  AdapterMatch am(read.seq, opt.MATCH_SCORE,
      opt.MISMATCH_SCORE, opt.GAP_SCORE);
  TrimStatus ts;
  unsigned newstart = am.TrimBeforeInsert(front_adap, ts, match_front_till, opt.debug);
  read.adap5 = cpputil::AnnotatedSeq(read.seq.substr(0, newstart), read.qual.substr(0, newstart));
  read.seq = read.seq.substr(newstart);
  read.qual = read.qual.substr(newstart);
  read.barcode = front_adap;

  //trim after insert
  AdapterMatch am2(read.seq, opt.MATCH_SCORE, opt.MISMATCH_SCORE, opt.GAP_SCORE);
  TrimStatus ts2;
  unsigned mstart, mend;
  std::tie(mstart, mend)  = am2.TrimAfterInsert(back_adap, back_adap.size(),
      ts2, 0, 0, false, 4, opt.debug);
  //if (ts2.tm != TRIMEND && mstart != read.seq.size()) ts2.tm = TRIM;
  if (mstart < read.seq.size()) {
    read.adap3 = cpputil::AnnotatedSeq(read.seq.substr(mstart, mend - mstart), read.qual.substr(mstart, mend - mstart));
    if (mend < read.seq.size()) {
      read.trim3 = cpputil::AnnotatedSeq(read.seq.substr(mend), read.qual.substr(mend));
    }
  }
  read.seq = read.seq.substr(0, mstart);
  read.qual = read.qual.substr(0, mstart);
  return ts;
}

template<typename Options>
TrimStatus TrimJunctions(const std::string& junc, const Options& opt, const std::string& illu_adap,
                         FxRecord& read) {
  AdapterMatch am(read.seq, opt.MATCH_SCORE, opt.MISMATCH_SCORE, opt.GAP_SCORE, opt.THREE_END_ADAPTER_MIN_LEN, opt.JUNCTION_ADAPTER_MIN_LEN);
  TrimStatus ts;
  // The trimmode are one of the following:
  // [notrim]: No trim or Less than 13bp are trimmed off the 3'end..
  // [trim]: The most confident trim. Linker is found at the middle of the read.
  // [trim13]:

  auto read_len =  am.TwoStepTrim(junc, illu_adap, ts, opt.nbase_skip, opt.debug);
  if (ts.tm == DOUBLE_LIGATION) {
    read.adap5 = cpputil::AnnotatedSeq(read.seq.substr(0, read_len), read.qual.substr(0, read_len));
    if (read_len == read.seq.size()) {
      read.seq = "";
      read.qual = "";
      return ts;
    }
    // if the leftover is not empty, try to trim off any adapter..
    read.seq = read.seq.substr(read_len);
    read.qual = read.qual.substr(read_len);
    DP_Matrix dp(junc.substr(0, TRUESEQ_PREFIX_LEN), read.seq);
    auto mb = dp.Compute(opt.MATCH_SCORE, opt.MISMATCH_SCORE, opt.GAP_SCORE,true, false);
    dp.BackTrace(mb, true, true, false);
    if (opt.debug == 1) {
      std::cerr << "double ligation\n";
      mb.Print(true);
    }

    bool found_13bp_adpt = mb.start_i == 0 && mb.nm <= 2 ? true : false;
    if (!found_13bp_adpt) {
      dp.BackTrace(mb, false, true, false);
      if (mb.start_i == 0 && (int) mb.max_j == dp.ncol_ -1 && mb.max_score >= (int) opt.THREE_END_ADAPTER_MIN_LEN * opt.MATCH_SCORE) {
        found_13bp_adpt = true;
        if (opt.debug == 1) {
          std::cerr << "found 13bp adpater at the end the read\n";
          mb.Print();
        }
      }
    }
    if (found_13bp_adpt) {
      if (opt.debug == 1) {
        std::cerr << "trim 13bp adapter..\n";
      }
      read.adap3 = cpputil::AnnotatedSeq(read.seq.substr(mb.start_j), read.qual.substr(mb.start_j));
      read.seq = read.seq.substr(0, mb.start_j);
      read.qual = read.qual.substr(0, mb.start_j);
    }
    return ts;
  }

  if (read_len < read.seq.size()) {
    auto shared_region_seq = read.seq.substr(read_len, ts.adap_end - read_len);
    auto shared_region_qual = read.qual.substr(read_len, ts.adap_end - read_len);
    read.adap3 = cpputil::AnnotatedSeq(shared_region_seq, shared_region_qual);
    if (ts.adap_end < (int) read.seq.size()) {
      auto unique_region_seq = read.seq.substr(ts.adap_end);
      auto unique_region_qual = read.qual.substr(ts.adap_end);
      read.trim3 = cpputil::AnnotatedSeq(unique_region_seq, unique_region_qual);
    }
  }
  if (read_len < read.seq.size()) {
    read.seq = read.seq.substr(0, read_len);
    read.qual = read.qual.substr(0, read_len);
  }
  if (ts.tm == TRIMEND || ts.tm == NOTRIM) {
    ts.tm = TRIM_INSUF;
  }
  return ts;
}

template<typename Options>
FragTrimS TrimPairedRead(const Options& opt, FxRecord& read1, FxRecord& read2) {
  //Pre Adapter Global Trim, 5end only
  if (opt.R1_UMI_LEN > 0 || opt.R2_UMI_LEN > 0) {
    TrimUmi(opt.R1_UMI_LEN, opt.R2_UMI_LEN, read1, read2);
  }
  FragTrimS fts;
  if (opt.linker_type != "adapter_v2") { //for obsolete adapter type
    if (LINKERS.find(opt.linker_type) == LINKERS.end()) throw std::runtime_error("Unknown Linker type");
    if (opt.debug == 1) std::cerr << read1.id <<"\n";
    std::string junc = LINKERS[opt.linker_type];
    fts.first = TrimJunctions(junc, opt, ILLUMINA_PAIR1_ADAPTER_SUFFIX, read1);
    if (fts.first.adap_dirct == AdapFwd) {
      cpputil::reverse_complement(junc);
    }
    if (opt.debug == 1) std::cerr << read2.id << "\n";
    fts.second = TrimJunctions(junc, opt, ILLUMINA_PAIR2_ADAPTER_SUFFIX, read2);
    if (fts.first.adap_dirct == AdapRev || fts.second.adap_dirct == AdapRev) {
      read1.rc_adpt = 1;
      read2.rc_adpt = 1;
    }
  }
  else {
    const int match_front_till = 30;
    std::string r1_pre = read1.index_barcode();
    std::string r1_post = read2.index_barcode();
    cpputil::reverse_complement(r1_post);
    std::string r2_pre = read2.index_barcode();
    std::string r2_post = read1.index_barcode();
    cpputil::reverse_complement(r2_post);
    if (opt.debug == 1) std::cerr << read1.id <<"\n";
    fts.first = TrimFrontAndBack(r1_pre, r1_post, opt, read1, match_front_till);
    if (opt.debug == 1) std::cerr << read2.id << "\n";
    fts.second = TrimFrontAndBack(r2_pre, r2_post, opt, read2, match_front_till);
  }

  //Post Adapter Global Trim
  if ((int) read1.seq.size() > opt.POST_TRIM3 + opt.POST_TRIM5){
    read1.seq = read1.seq.substr(opt.POST_TRIM5, read1.seq.size() - opt.POST_TRIM5 - opt.POST_TRIM3);
    read1.qual = read1.qual.substr(opt.POST_TRIM5, read1.qual.size() - opt.POST_TRIM5 - opt.POST_TRIM3);
  }
  if ((int) read2.seq.size() > opt.POST_TRIM3 + opt.POST_TRIM5){
    read2.seq = read2.seq.substr(opt.POST_TRIM5, read2.seq.size() - opt.POST_TRIM5 - opt.POST_TRIM3);
    read2.qual = read2.qual.substr(opt.POST_TRIM5, read2.qual.size() - opt.POST_TRIM5 - opt.POST_TRIM3);
  }

  return fts;
}

};

#endif //ADAPTERTRIM_INCLUDE_ADAPTER_H_
