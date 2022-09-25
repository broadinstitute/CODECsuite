//
// Created by Ruolin Liu on 8/16/21.
//

#ifndef CODECSUITE_BBCPPUTIL_INCLUDE_MUTCOUNTER_H_
#define CODECSUITE_BBCPPUTIL_INCLUDE_MUTCOUNTER_H_

#include <vector>
#include <map>
#include "SeqLib/BamRecord.h"
#include "SeqLib/RefGenome.h"
#include "DNAUtils.h"
#include "BamRecordExt.h"

namespace cpputil{
static int PAD_5 = 5;
struct ErrorStat {
  int64_t neval = 0;
  int64_t nsnv_error = 0;
  int64_t nindel_error = 0;
  int64_t indel_nbase_error = 0;
  int64_t nsnv_masked_by_vcf1 = 0;
  int64_t nindel_masked_by_vcf1 = 0;
  int64_t nsnv_masked_by_maf = 0;
  int64_t nindel_masked_by_maf = 0;
  int64_t discard_frag_counter = 0;
  int64_t n_pass_filter_pairs = 0;
  int64_t n_pass_filter_singles = 0;
  std::vector<int> cutoffs;
  std::vector<int64_t> R1_q0_cov;
  std::vector<int64_t> R1_q30_cov;
  std::vector<int64_t> R1_q0_error;
  std::vector<int64_t> R1_q30_error;
  std::vector<int64_t> R2_q0_cov;
  std::vector<int64_t> R2_q30_cov;
  std::vector<int64_t> R2_q0_error;
  std::vector<int64_t> R2_q30_error;
  std::map<char, int64_t> base_counter;
  std::map<std::string, int64_t> triplet_counter;
  std::map<std::string, int64_t> doublet_counter;
  //qscore cutoff -> read1, read2
  std::map<int, std::pair<int64_t, int64_t>> qcut_neval;
  std::map<int, std::pair<int64_t, int64_t>> qcut_nerrors;
  int n_filtered_sclip = 0;
  int n_filtered_Nrate = 0;
  int n_filtered_q30rate = 0;
  int n_filtered_smallfrag = 0;
  int n_filtered_largefrag = 0;
  int n_filtered_edit = 0;
  int n_filtered_clustered = 0;
  int nindel_filtered_adjbaseq = 0;
  int nindel_filtered_adjN = 0;
  int nindel_filtered_adjvar = 0;
  int nindel_filtered_overlap_snp = 0;
  int mismatch_filtered_by_indel = 0;
  int snv_family_disagree = 0;
  int snv_R1R2_disagree = 0;
  int snv_filtered_baseq = 0;
  int indel_R1R2_disagree = 0;
  int lowconf_t2g = 0;
  int low_germ_depth = 0;
  int seen_in_germ = 0;
  int AS_filter = 0;

  //ErrorStat() = delete;
  ErrorStat(int L, int q) {
    cutoffs = {0};
    //cutoffs = {0, 30};
    qcut_neval[0] = std::make_pair(0,0);
    qcut_neval[q] = std::make_pair(0,0);
    qcut_nerrors[0] = std::make_pair(0,0);
    qcut_nerrors[q] = std::make_pair(0,0);
    Init(L);
  }
  ErrorStat(const vector<int>& qcuts, int L): cutoffs(qcuts){
    assert(std::find(qcuts.begin(), qcuts.end(), 0) != qcuts.end());
    //assert(std::find(qcuts.begin(), qcuts.end(), 30) != qcuts.end());
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
  }
};

std::pair<int,int> CountValidBaseInMatchedBases(const SeqLib::BamRecord &b,
                                                const std::set<int>& site_blacklist,
                                                const int minbq,
                                                std::set<int>& bq_blacklist,
                                                int qstart,
                                                int qend,
                                                bool N_is_valid) {
  //Valid bases are not N and baseq >= qcutoff
  if (qend <= qstart) return std::make_pair(0, 0);
  if (qend > b.AlignmentEndPosition()) {
    qend = b.AlignmentEndPosition();
  }
  if (qstart < b.AlignmentPosition()) {
    qstart = b.AlignmentPosition();
  }
  uint32_t *c = bam_get_cigar(b.raw());
  int32_t readpos = 0;
  int32_t refpos = b.Position();
  size_t i = 0;

  for (; i < b.raw()->core.n_cigar; ++i) {
    if (bam_cigar_opchr(c[i]) == 'S')
      readpos += bam_cigar_oplen(c[i]);
    else if (bam_cigar_opchr(c[i]) != 'H')
      break;
  }
  int res = 0, q0res=0;
  const uint8_t *p = bam_get_seq(b.raw());
  const uint8_t *bq = bam_get_qual(b.raw());
  int cur = 0;
  for (; i < b.raw()->core.n_cigar; ++i) {
    char cigar = bam_cigar_opchr(c[i]);
    if (cigar == 'M' or cigar == 'X' or cigar == '=') {
      for(unsigned ww = 0; ww < bam_cigar_oplen(c[i]); ++ww) {
        cur = ww + readpos;
        if (cur <  qstart) continue;
        if (cur >= qend) break;
        if (not site_blacklist.empty() and site_blacklist.find(refpos + ww) != site_blacklist.end()) continue;
        //if ((N_is_valid && bam_seqi(p,cur) == 15) || bam_seqi(p,cur) != 15) ++q0res;
        ++q0res;
        if ((not N_is_valid && bam_seqi(p,cur) == 15) or bq[cur] < minbq) {
          bq_blacklist.insert(refpos + ww);
          continue;
        }
        if (not bq_blacklist.empty() and bq_blacklist.find(refpos + ww) != bq_blacklist.end()) continue;
        ++res;
      }
      readpos += bam_cigar_oplen(c[i]);
      refpos += bam_cigar_oplen(c[i]);
    } else if(cigar == 'I') {
      readpos += bam_cigar_oplen(c[i]);
    } else if(cigar == 'D') {
      if (readpos >= qend) break;
      int dl = bam_cigar_oplen(c[i]);
      if (readpos >= qstart) {
        for (int l = 0; l < dl; ++l ) {
          bq_blacklist.insert(refpos + l);
        }
      }
      refpos += dl;
    }
  }
  return std::make_pair(res, q0res);
}

std::pair<int, int> CountValidBaseAndContextInMatchedBases(const SeqLib::BamRecord &b,
                                                           const std::set<int>& site_blacklist,
                                                           const std::string chrname,
                                                           const SeqLib::RefGenome& refgenome,
                                                           const int minbq,
                                                           std::set<int>& bq_blacklist,
                                                           ErrorStat& es,
                                                           int qstart,
                                                           int qend,
                                                           bool N_is_valid) {
  //Valid bases are not N and baseq >= qcutoff
  if (qend <= qstart) return std::make_pair(0, 0);
  if (qend > b.AlignmentEndPosition()) {
    qend = b.AlignmentEndPosition();
  }
  if (qstart < b.AlignmentPosition()) {
    qstart = b.AlignmentPosition();
  }
  uint32_t *c = bam_get_cigar(b.raw());
  int32_t readpos = 0;
  int32_t refpos = b.Position();
  int32_t rstart = std::max(refpos - 1, 0);
  std::string ref = refgenome.QueryRegion(chrname, rstart, b.PositionEnd()); // 1bp padding in both ends
  size_t i = 0;
  for (; i < b.raw()->core.n_cigar; ++i) {
    if (bam_cigar_opchr(c[i]) == 'S')
      readpos += bam_cigar_oplen(c[i]);
    else if (bam_cigar_opchr(c[i]) != 'H')
      break;
  }

  int res = 0, q0res=0;
  const uint8_t *p = bam_get_seq(b.raw());
  const uint8_t *bq = bam_get_qual(b.raw());
  int cur = 0;
  for (; i < b.raw()->core.n_cigar; ++i) {
    char cigar = bam_cigar_opchr(c[i]);
    if (cigar == 'M' or cigar == '=' or cigar == 'X') {
      for(unsigned ww = 0; ww < bam_cigar_oplen(c[i]); ++ww) {
        cur = ww + readpos;
        if (cur <  qstart) continue;
        if (cur >= qend) break;
        if (not site_blacklist.empty() and site_blacklist.find(refpos + ww) != site_blacklist.end()) continue;
        if ((N_is_valid && bam_seqi(p,cur) == 15) || bam_seqi(p,cur) != 15) ++q0res;
        if ((not N_is_valid && bam_seqi(p,cur) == 15) or bq[cur] < minbq) {
          bq_blacklist.insert(refpos + ww);
          continue;
        }
        if (not bq_blacklist.empty() and bq_blacklist.find(refpos + ww) != bq_blacklist.end()) continue;
        ++res;
        es.base_counter[toupper(ref[refpos + ww - rstart])]++;
        if (refpos + ww - rstart -1 >= 0) {
          auto tri = ref.substr(refpos + ww - rstart - 1, 3);
          std::transform(tri.begin(), tri.end(), tri.begin(), ::toupper);
          ++es.triplet_counter[tri];
        }

        if (bq[cur + 1] >= minbq &&
            site_blacklist.find(refpos + ww + 1) == site_blacklist.end() &&
            bq_blacklist.find(refpos + ww + 1) == bq_blacklist.end()) {
          auto di = ref.substr(refpos + ww - rstart, 2);
          std::transform(di.begin(), di.end(), di.begin(), ::toupper);
          ++es.doublet_counter[di];
        }
      }
      readpos += bam_cigar_oplen(c[i]);
      refpos += bam_cigar_oplen(c[i]);
    } else if(cigar == 'I') {
      readpos += bam_cigar_oplen(c[i]);
    } else if(cigar == 'D') {
      if (readpos >= qend) break;
      int dl = bam_cigar_oplen(c[i]);
      if (readpos >= qstart) {
        for (int l = 0; l < dl; ++l ) {
          bq_blacklist.insert(refpos + l);
        }
      }
      refpos += dl;
    }
  }
  return std::make_pair(res, q0res);
}

static std::pair<int,int> CountDenom(const cpputil::Segments& seg,
                              const SeqLib::GenomicRegion* const gr,
                              const SeqLib::RefGenome& ref,
                              const string& chrname,
                              const std::set<int>& blacklist,
                              cpputil::ErrorStat& es,
                              int minbq,
                              std::pair<int, int>& nq0,
                              bool count_context,
                              int count_overhang,
                              bool N_is_valid){
  // blacklist represents a set of SNV positions that will not be counted in the error rate calculation
  int r1 = 0, r2 = 0, r1q0 = 0, r2q0 = 0;
  std::set<int> baseqblack;
  if (seg.size() == 1) {
    std::pair<int,int> range;
    const auto& br = seg.front();
    if (gr) {
      range = cpputil::GetBamOverlapQStartAndQStop(br, *gr);
    } else {
      range.first = br.AlignmentPosition();
      range.second = br.AlignmentEndPosition();
    }
    if (not br.PairedFlag() or br.FirstFlag()) {
      if (count_context)
        std::tie(r1,r1q0) = cpputil::CountValidBaseAndContextInMatchedBases(br, blacklist, chrname, ref, minbq, baseqblack, es, range.first, range.second, N_is_valid);
      else
        std::tie(r1, r1q0) = cpputil::CountValidBaseInMatchedBases(br, blacklist, minbq, baseqblack,range.first, range.second, N_is_valid);
    }
    else {
      if (count_context)
        std::tie(r2, r2q0) = cpputil::CountValidBaseAndContextInMatchedBases(br, blacklist, chrname, ref, minbq, baseqblack, es, range.first, range.second, N_is_valid);
      else
        std::tie(r2, r2q0) = cpputil::CountValidBaseInMatchedBases(br, blacklist, minbq, baseqblack, range.first, range.second, N_is_valid);
    }
  } else {
    std::pair<int, int> overlap_front, overlap_back;
    if (count_overhang == 1 && !IsPairOverlap(seg[0], seg[1])) {
      count_overhang = 2;
    }
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
    if (count_overhang == 2) {
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
      size_t r1_nmask;
      if (count_context) {
        if (count_overhang == 2) {
          std::tie(r1, r1q0) = cpputil::CountValidBaseAndContextInMatchedBases(seg.front(), blacklist, chrname, ref, minbq, baseqblack, es, range_front.first, range_front.second, N_is_valid);
        } else if (count_overhang == 1){
          if (seg.front().ReverseFlag()) {
            std::tie(r2, r2q0) = cpputil::CountValidBaseAndContextInMatchedBases(seg.back(), blacklist, chrname, ref, minbq, baseqblack, es, overlap_back.first, range_back.first, N_is_valid);
            baseqblack.clear();
            cpputil::CountValidBaseInMatchedBases(seg.back(), blacklist, minbq, baseqblack, range_back.first, range_back.second, N_is_valid); // for blacklist only
          } else {
            std::tie(r1, r1q0) = cpputil::CountValidBaseAndContextInMatchedBases(seg.front(), blacklist, chrname, ref, minbq, baseqblack, es, overlap_front.first, range_front.first, N_is_valid);
            baseqblack.clear();
            cpputil::CountValidBaseInMatchedBases(seg.front(), blacklist, minbq, baseqblack, range_front.first, range_front.second, N_is_valid); // for blacklist only
          }
        } else {
          std::tie(r1, r1q0) = cpputil::CountValidBaseInMatchedBases(seg.front(), blacklist, minbq, baseqblack, range_front.first, range_front.second, N_is_valid);
        }
        r1_nmask = baseqblack.size();
        if (count_overhang == 2) baseqblack.clear();// when count_overhang is true, treat R1 and R2 independently
        if (count_overhang == 1) {
          if (seg.front().ReverseFlag()) {
            std::tie(r1, r1q0) = cpputil::CountValidBaseAndContextInMatchedBases(seg.front(), blacklist, chrname, ref, minbq, baseqblack, es, range_front.first, range_front.second, N_is_valid);
            std::set<int> baseqblack2;
            auto r1_res = cpputil::CountValidBaseAndContextInMatchedBases(seg.front(), blacklist, chrname, ref, minbq, baseqblack2, es, range_front.second, overlap_front.second, N_is_valid);
            r1 += r1_res.first;
            r1q0 += r1_res.second;
          } else {
            std::tie(r2, r2q0) = cpputil::CountValidBaseAndContextInMatchedBases(seg.back(), blacklist, chrname, ref, minbq, baseqblack, es, range_back.first, range_back.second, N_is_valid);
            std::set<int> baseqblack2;
            auto r2_res = cpputil::CountValidBaseAndContextInMatchedBases(seg.back(), blacklist, chrname, ref, minbq, baseqblack2, es, range_back.second, overlap_back.second, N_is_valid);
            r2 += r2_res.first;
            r2q0 += r2_res.second;
          }
        } else {
          std::tie(r2, r2q0) = cpputil::CountValidBaseAndContextInMatchedBases(seg.back(), blacklist, chrname, ref, minbq, baseqblack, es, range_back.first, range_back.second, N_is_valid);
        }
      }
      else { // not count context
        if (count_overhang == 1) {
          if (seg.front().ReverseFlag()) {
            std::tie(r2, r2q0) = cpputil::CountValidBaseInMatchedBases(seg.back(), blacklist, minbq, baseqblack, overlap_back.first, range_back.first, N_is_valid);
            baseqblack.clear();
            cpputil::CountValidBaseInMatchedBases(seg.back(), blacklist, minbq, baseqblack, range_back.first, range_back.second, N_is_valid);
          }
          else {
            std::tie(r1, r1q0) = cpputil::CountValidBaseInMatchedBases(seg.front(), blacklist, minbq, baseqblack, overlap_front.first, range_front.first, N_is_valid);
            // for blacklist only
            baseqblack.clear();
            cpputil::CountValidBaseInMatchedBases(seg.front(), blacklist, minbq, baseqblack, range_front.first, range_front.second, N_is_valid);
          }
        } else {
          std::tie(r1, r1q0) = cpputil::CountValidBaseInMatchedBases(seg.front(), blacklist, minbq, baseqblack, range_front.first, range_front.second, N_is_valid);
        }
        r1_nmask = baseqblack.size();
        if (count_overhang == 2) baseqblack.clear();
        if (count_overhang == 1) {
          if (seg.front().ReverseFlag()) {
            std::tie(r1, r1q0) = cpputil::CountValidBaseInMatchedBases(seg.front(), blacklist, minbq, baseqblack, range_front.first, range_front.second, N_is_valid);
            std::set<int> baseqblack2;
            auto r1_res = cpputil::CountValidBaseInMatchedBases(seg.front(), blacklist, minbq, baseqblack2, range_front.second, overlap_front.second, N_is_valid);
            r1 += r1_res.first;
            r1q0 += r1_res.second;
          } else {
            std::tie(r2, r2q0) = cpputil::CountValidBaseInMatchedBases(seg.back(), blacklist, minbq, baseqblack, range_back.first, range_back.second, N_is_valid);
            std::set<int> baseqblack2;
            auto r2_res = cpputil::CountValidBaseInMatchedBases(seg.back(), blacklist, minbq, baseqblack2, range_back.second, overlap_back.second, N_is_valid);
            r2 += r2_res.first;
            r2q0 += r2_res.second;
          }
        } else {
          std::tie(r2, r2q0) = cpputil::CountValidBaseInMatchedBases(seg.back(), blacklist, minbq, baseqblack, range_back.first, range_back.second, N_is_valid);
        }
      }
      if (count_overhang == 0) r1 -= baseqblack.size() - r1_nmask;
    } else {
      throw std::runtime_error("Read order wrong\n");
    }
  }
  nq0 = std::make_pair(r1q0, r2q0);
  return std::make_pair(r1, r2);
}

std::pair<int,int> NumEffectBases(const cpputil::Segments& seg, int minbq, bool count_overhang, bool N_is_valid) {
  // blacklist represents a set of SNV positions that will not be counted in the error rate calculation
  int r1 = 0, r2 = 0, r1q0 = 0, r2q0 = 0;
  std::set<int> baseqblack;
  std::set<int> blacklist;
  if (seg.size() == 1) {
    std::pair<int,int> range;
    const auto& br = seg.front();
    range.first = br.AlignmentPosition();
    range.second = br.AlignmentEndPosition();
    if (not br.PairedFlag() or br.FirstFlag()) {
      std::tie(r1, r1q0) = CountValidBaseInMatchedBases(br, blacklist, minbq, baseqblack,range.first, range.second, N_is_valid);
    }
    else {
      std::tie(r2, r2q0) = CountValidBaseInMatchedBases(br, blacklist, minbq, baseqblack, range.first, range.second, N_is_valid);
    }
  } else {
    std::pair<int, int> overlap_front, overlap_back;
    overlap_front.first = seg.front().AlignmentPosition();
    overlap_front.second = seg.front().AlignmentEndPosition();
    overlap_back.first = seg.back().AlignmentPosition();
    overlap_back.second = seg.back().AlignmentEndPosition();
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
      size_t r1_nmask;
      std::tie(r1, r1q0) = cpputil::CountValidBaseInMatchedBases(seg.front(), blacklist, minbq, baseqblack, range_front.first, range_front.second, N_is_valid);
      r1_nmask = baseqblack.size();
      if (count_overhang) baseqblack.clear();
      std::tie(r2, r2q0) = cpputil::CountValidBaseInMatchedBases(seg.back(), blacklist, minbq, baseqblack, range_back.first, range_back.second, N_is_valid);
      if (not count_overhang) r1 -= baseqblack.size() - r1_nmask;
    } else {
      throw std::runtime_error("Read order wrong\n");
    }
  }
  return std::make_pair(r1, r2);
}

bool HasBadCigar(const SeqLib::BamRecord& rec) {
  //If Indel is next to soft clipping
  const auto cigar = rec.GetCigar();
  auto s = cigar.size();
  if (s == 0) return true;
  if (s > 1 && cigar[0].Type() == 'S') {
    if (cigar[1].Type() == 'I' or cigar[1].Type() == 'D') return true;
  }
  if (s > 1 && cigar[s-1].Type() == 'S') {
    if (cigar[s-2].Type() == 'I' or cigar[s-2].Type() == 'D') return true;
  }
  return false;
}

template<typename Options>
int FailFilter(const vector<cpputil::Segments>& frag,
               const SeqLib::BamHeader& bamheader,
               const SeqLib::RefGenome& ref,
               const Options& opt,
               ErrorStat& errorstat,
                    bool paired_only,
                    bool indel_calls_only,
                    int& frag_numN,
                    int& nqpass,
                    int& olen) {
  /*
   * Pass =0, >0 fail mode
   */
  if (indel_calls_only) {
    bool has_indel = false;
    for (const auto &seg: frag) {
      for (const auto &r : seg)
      if (IndelLen(r) > 0) {
        has_indel = true;
        break;
      }
    }
    if (not has_indel) {
      return 888;
    }
  }
  Segments const *seg = nullptr;
  if (frag.size() == 1) {
    seg = &frag[0];
  } else {
    for (const auto& ss : frag) {
      if (not ss[0].DuplicateFlag() and not ss[0].SecondaryFlag()) {
        seg = &ss;
      }
    }
  }
  if (not seg) {return 999;} // if all reads are duplicates
  if (seg->size() > 2 || (paired_only && seg->size() != 2)) {
    for (auto s : *seg) {
      std::cerr << s << std::endl;
    }
    throw std::runtime_error("Unexpected #reads for a fragment");
  }

  if (seg->size() == 1 || (*seg)[0].InsertSize() == 0) {
    if ((int) (*seg)[0].NumMatchBases() < std::max(1, opt.min_fraglen)) {
      ++errorstat.n_filtered_smallfrag;
      return 1;
    }
  } else {
    if (abs((*seg)[0].InsertSize()) < std::max(1, opt.min_fraglen)) {
      ++errorstat.n_filtered_smallfrag;
      return 1;
    }
  }

  olen = EffFragLen(*seg, opt.count_read);
  const std::set<int> blacklist;
  std::pair<int, int> q0den(0, 0);
  auto qpass = CountDenom(*seg, nullptr, ref, "", blacklist, errorstat, opt.bqual_min,
                          q0den, false, opt.count_read, false);

  if (seg->size() == 2) {
    if (opt.filter_5endclip && (cpputil::NumSoftClip5End((*seg)[0]) > 0 || cpputil::NumSoftClip5End((*seg)[1]) > 0)) {
      ++errorstat.n_filtered_sclip;
      return 2;
    }
    nqpass = opt.count_read? qpass.first + qpass.second : qpass.first;
    frag_numN = seg->front().CountNBases();
    frag_numN = std::max(frag_numN, seg->back().CountNBases());
  } else {
    if (opt.filter_5endclip && cpputil::NumSoftClip5End((*seg)[0]) > 1) {
      ++errorstat.n_filtered_sclip;
      return 2;
    }
    nqpass = std::max(qpass.first, qpass.second);
    frag_numN = seg->front().CountNBases();
  }

  if (nqpass < olen * opt.min_passQ_frac) {
    ++errorstat.n_filtered_q30rate;
    return 3;
  }

  if (frag_numN > abs(seg->front().InsertSize()) * opt.max_N_frac || frag_numN > opt.max_N_filter) {
    ++errorstat.n_filtered_Nrate;
    return 4;
  }

  for (const auto &s: *seg) {

    if ((int) s.NumMatchBases() > opt.max_fraglen || abs(s.InsertSize()) > opt.max_fraglen) {
      ++errorstat.n_filtered_largefrag;
      return 5;
    }

    if (cpputil::GetNMismatch(s) > opt.max_snv_filter) {
      ++errorstat.n_filtered_edit;
      return 6;
    }

    if (opt.clustered_mut_cutoff <  std::numeric_limits<int>::max() && cpputil::HasClusteredMuts(s, bamheader, ref, opt.clustered_mut_cutoff)) {
      ++errorstat.n_filtered_clustered;
      return 7;
    }

    if (HasBadCigar(s)) {
      ++errorstat.AS_filter;
      return 9;
    }

    int XS, AS;
    bool xsstat = s.GetIntTag("XS", XS);
    bool asstat = s.GetIntTag("AS", AS);
    if (xsstat and asstat) {
      if (XS >= AS * opt.max_frac_prim_AS) {
        ++errorstat.AS_filter;
        return 8;
      }
    }
  }

  return 0;
}


}//end namespace

#endif //CODECSUITE_BBCPPUTIL_INCLUDE_MUTCOUNTER_H_
