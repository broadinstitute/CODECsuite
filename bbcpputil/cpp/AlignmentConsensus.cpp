//
// Created by Ruolin Liu on 3/19/20.
//

#include "AlignmentConsensus.h"


namespace cpputil {

std::string GetConsensusTemplate(const Segments& segs, int32_t& ref_most_left) {
  /*
   *  '.' : uninitialized
   *  '+' : insertion
   */
  ref_most_left = std::numeric_limits<int32_t>::max();
  int32_t ref_most_right = 0;
  for (auto & seg : segs) {
    ref_most_left = std::min(seg.PositionWithSClips(), ref_most_left);
    ref_most_right = std::max(seg.PositionEndWithSClips(), ref_most_right);
  }
  int32_t ref_span = ref_most_right - ref_most_left;
  assert(ref_span > 0);
  std::map<int,int> ins_len;
  for (auto & seg: segs) {
    const SeqLib::Cigar cigar = seg.GetCigar();
    int32_t relative_start = seg.PositionWithSClips() - ref_most_left;
    find_insert_(cigar, relative_start, ins_len);
  }
//  const SeqLib::Cigar left_cigar = seg.front().GetCigar();
//  const SeqLib::Cigar right_cigar = seg.back().GetCigar();
//  find_insert_(left_cigar, 0, ins_len);
//  find_insert_(right_cigar, seg.back().PositionWithSClips() - seg.front().PositionWithSClips(), ins_len);
  std::string consens;
  int b = 0;
  for (const auto pos_len : ins_len) {
    consens += std::string(pos_len.first - b, '.');
    consens += std::string(pos_len.second, '+');
    b = pos_len.first;
  }
  consens += std::string(ref_span - b, '.');
  return consens;
}

std::pair<std::string, std::string>
    GetGappedSeqAndQual(const SeqLib::BamRecord &r, const std::string& seq, const int start, const std::string consensus_template) {

  // seq should be same length as r.Sequence() but can with difference sequence
  // `.` is overhang
  // `+` is insertion
  // `-` is deletion
  int template_start = 0;
  std::string consns_templ = consensus_template;
  std::string quality_templ = std::string(consensus_template.size(), 33);
  int ref_count = 0;
  for (char c : consns_templ) {
    if (ref_count == start) break;
    template_start++;
    if (c == '.') ref_count++;
  }
  auto const &cigar = r.GetCigar();
  auto const &qual = r.Qualities();

  int read_start = 0;
  for (auto c = cigar.begin(); c != cigar.end(); ++c) {
    if (c->Type() == 'S' || c->Type() == 'M' || c->Type() == 'D') {
      for (unsigned ii = 0; ii < c->Length(); ++ii) {
        while ('+' == consns_templ[template_start]) { template_start++; };
        char s, q;
        if (c->Type() == 'D') {
          s = '-';
          q = qual[read_start];
        } else {
          s = seq[read_start];
          q = qual[read_start++];
        }
        if (consns_templ[template_start] != '+') {
          consns_templ[template_start] = s;
          quality_templ[template_start++] = q;
        }
      }
    } else if (c->Type() == 'I') {
      for (unsigned ii = 0; ii < c->Length(); ++ii) {
        consns_templ[template_start] = seq[read_start];
        quality_templ[template_start++] = qual[read_start++];
      }
      while (consns_templ[template_start] == '+') { template_start++; };
    }
  }
  return std::make_pair(consns_templ, quality_templ);
}

std::string MergePairSeq(const Segments &segs, const std::vector<std::string>& seqs, bool trim_overhang) {
  //seqs should hold the fastq seq for segs. They could be same length but different seqs
  //not to change indel
  assert (segs.size() == 2);
  int ref_most_left;
  const char NUL = 6;
  const std::string consns_templ = GetConsensusTemplate(segs, ref_most_left);
  std::string consns_seq1(consns_templ.size(), NUL);
  int32_t start = 0;
  std::string dnaseq, qual;
  std::vector<std::string> dna_pileup(segs.size());
  std::vector<std::string> qual_pileup(segs.size());
  for (unsigned sid = 0; sid < segs.size(); ++sid) {
    auto const& seg = segs[sid];
    start = seg.PositionWithSClips() - ref_most_left;
    std::tie(dnaseq, qual) = GetGappedSeqAndQual(seg, seqs[sid], start, consns_templ);
    dna_pileup[sid] = dnaseq;
    qual_pileup[sid] = qual;
  }
  for (unsigned jj = 0; jj < consns_templ.size(); ++jj) {
    // paired baseq calibration. If only one of the baseq < cutoff, make the other one baseq = cutoff -1
    // so that when later we filter by baseq by this cutoff, they either both stay or both out
    if (dna_pileup[0][jj] == '.' or dna_pileup[1][jj] == '.') { // overhang
      if (not trim_overhang) {
        if (dna_pileup[0][jj] != '.' and dna_pileup[0][jj] != '-' and dna_pileup[1][jj] == '.') {
          consns_seq1[jj] = dna_pileup[0][jj];
        } else if (dna_pileup[0][jj] == '.' and dna_pileup[1][jj] != '.' and dna_pileup[1][jj] != '-') {
          consns_seq1[jj] = dna_pileup[1][jj];
        }
      }
    } else {
      if (dna_pileup[0][jj] != dna_pileup[1][jj]) {
        assert(dna_pileup[0][jj] != '-' or dna_pileup[1][jj] != '+');
        assert(dna_pileup[0][jj] != '+' or dna_pileup[1][jj] != '-');
        if (dna_pileup[0][jj] == '-' or dna_pileup[0][jj] == '+') {
          consns_seq1[jj] = dna_pileup[1][jj];
        }
        else if (dna_pileup[1][jj] == '-' or dna_pileup[1][jj] == '+') {
          consns_seq1[jj] = dna_pileup[0][jj];
        }
        else {
          if (qual_pileup[0][jj] < qual_pileup[1][jj]) {
            consns_seq1[jj] = dna_pileup[1][jj];
          } else {
            consns_seq1[jj] = dna_pileup[0][jj];
          }
        }

      } else if (dna_pileup[0][jj] >= 'A') {
        consns_seq1[jj] = dna_pileup[0][jj];
      } else if (dna_pileup[0][jj] == '+') {
        assert(false);
      }
    }
  }
  consns_seq1.erase(std::remove(consns_seq1.begin(), consns_seq1.end(), NUL), consns_seq1.end());
  return consns_seq1;
}

std::string MergePair(const Segments &seg, bool trim_overhang) {
  std::vector<std::string> seqs;
  std::vector<std::string> dummy_quals;
  for (auto&s : seg) {
    seqs.push_back(s.Sequence());
  }
  auto seq = MergePairSeq(seg, seqs, trim_overhang);
  return seq;
}

std::pair<std::string, std::string> PairSeqConsensus(const Segments &seg, bool trim_overhang, int qcutoff) {
  std::vector<std::string> seqs;
  std::vector<std::string> dummy_quals;
  for (auto&s : seg) {
    seqs.push_back(s.Sequence());
  }
  auto seq = cpputil::PairConsensus(seg, seqs, trim_overhang, qcutoff, dummy_quals);
  return seq;
}

std::pair<std::vector<std::string>, std::vector<std::string>> GetPairPileup(const Segments &segs) {
  assert(segs.size() == 2);
  const char NUL = 6;
  int ref_most_left;
  const std::string consns_templ = GetConsensusTemplate(segs, ref_most_left);
  std::vector<std::string> seqs;
  for (auto&s : segs) {
    seqs.push_back(s.Sequence());
  }
  std::string consns_seq1(consns_templ.size(), NUL);
  std::string consns_seq2(consns_templ.size(), NUL);
  int32_t start = 0;
  std::string dnaseq, qual;
  std::vector<std::string> dna_pileup(segs.size());
  std::vector<std::string> qual_pileup(segs.size());
  for (unsigned sid = 0; sid < segs.size(); ++sid) {
    auto const& seg = segs[sid];
    start = seg.PositionWithSClips() - ref_most_left;
    std::tie(dnaseq, qual) = GetGappedSeqAndQual(seg, seqs[sid], start, consns_templ);
    dna_pileup[sid] = dnaseq;
    qual_pileup[sid] = qual;
  }
  return std::make_pair(dna_pileup, qual_pileup);
}

std::pair<std::string, std::string> PairConsensus(const Segments &segs, const std::vector<std::string>& seqs,
    bool trim_overhang, int qcutoff, std::vector<std::string>& out_quals) {
  //seqs should hold the fastq seq for segs. They could be same length but different seqs
  //not to change indel
  assert (segs.size() == 2);
  out_quals.resize(segs.size());
  int ref_most_left;
  const char NUL = 6;
  const std::string consns_templ = GetConsensusTemplate(segs, ref_most_left);
  std::string consns_seq1(consns_templ.size(), NUL);
  std::string consns_seq2(consns_templ.size(), NUL);
  for (auto & qual : out_quals) {
    qual = std::string(consns_templ.size(), 33);
  }
  int32_t start = 0;
  std::string dnaseq, qual;
  std::vector<std::string> dna_pileup(segs.size());
  std::vector<std::string> qual_pileup(segs.size());
  for (unsigned sid = 0; sid < segs.size(); ++sid) {
    auto const& seg = segs[sid];
    start = seg.PositionWithSClips() - ref_most_left;
    std::tie(dnaseq, qual) = GetGappedSeqAndQual(seg, seqs[sid], start, consns_templ);
    dna_pileup[sid] = dnaseq;
    qual_pileup[sid] = qual;
  }
  for (unsigned jj = 0; jj < consns_templ.size(); ++jj) {
    // paired baseq calibration. If one of the baseq < cutoff, make all baseq low enough so that VC will ingnore them
    if (dna_pileup[0][jj] >= 'A' && dna_pileup[1][jj] >= 'A'
        && std::min(qual_pileup[0][jj], qual_pileup[1][jj]) < static_cast<char>(33 + qcutoff)) {
      qual_pileup[0][jj] = static_cast<char>(35);
      qual_pileup[1][jj] = static_cast<char>(35);
    }
    if (dna_pileup[0][jj] == '.' or dna_pileup[1][jj] == '.') { // overhang
      if (not trim_overhang) {
        if (dna_pileup[0][jj] != '.' and dna_pileup[0][jj] != '-' and dna_pileup[1][jj] == '.') {
          consns_seq1[jj] = dna_pileup[0][jj];
          out_quals[0][jj] = qual_pileup[0][jj];
        } else if (dna_pileup[0][jj] == '.' and dna_pileup[1][jj] != '.' and dna_pileup[1][jj] != '-') {
          consns_seq2[jj] = dna_pileup[1][jj];
          out_quals[1][jj] = qual_pileup[1][jj];
        }
      }
    } else {
      if (dna_pileup[0][jj] != dna_pileup[1][jj]) {
        assert(dna_pileup[0][jj] != '-' or dna_pileup[1][jj] != '+');
        assert(dna_pileup[0][jj] != '+' or dna_pileup[1][jj] != '-');
        if (dna_pileup[0][jj] == '-' or dna_pileup[0][jj] == '+') {
          consns_seq2[jj] = dna_pileup[1][jj];
        }
        else if (dna_pileup[1][jj] == '-' or dna_pileup[1][jj] == '+') {
          consns_seq1[jj] = dna_pileup[0][jj];
        }
        else {
          consns_seq1[jj] = 'N';
          consns_seq2[jj] = 'N';
        }
#if 0
    // DO CONSENSUS for INDEL
        if (std::min(dna_pileup[0][jj], dna_pileup[1][jj]) == '-') {
          // the consensus of del and non-del is N
          //consns_seq[jj] = std::max(dna_pileup[0][jj], dna_pileup[1][jj]);
          consns_seq1[jj] = 'N';
          consns_seq2[jj] = 'N';
        } else if (std::min(dna_pileup[0][jj], dna_pileup[1][jj]) != '+') {
          // the consensus of ins and non-ins is N
          consns_seq1[jj] = 'N';
          consns_seq2[jj] = 'N';
        }
#endif

      } else if (dna_pileup[0][jj] >= 'A') {
        consns_seq1[jj] = dna_pileup[0][jj];
        consns_seq2[jj] = dna_pileup[0][jj];
      } else if (dna_pileup[0][jj] == '+') {
        assert(false);
      }
    }
    out_quals[0][jj] = consns_seq1[jj] == NUL ? NUL : qual_pileup[0][jj];
    out_quals[1][jj] = consns_seq2[jj] == NUL ? NUL : qual_pileup[1][jj];
  }
  consns_seq1.erase(std::remove(consns_seq1.begin(), consns_seq1.end(), NUL), consns_seq1.end());
  consns_seq2.erase(std::remove(consns_seq2.begin(), consns_seq2.end(), NUL), consns_seq2.end());
  for (auto& qual : out_quals) {
    qual.erase(std::remove(qual.begin(),qual.end(), NUL ), qual.end());
  }
  return std::make_pair(consns_seq1, consns_seq2);
}


// Simple consensus by requiring all bases to be the same, otherwise 'N'.
// Same for INS and DEL. The Inserted seq has to be the same. The insertions with different lengths are truncated the the
// smallest length.
// The consensus of a deleted base and undeleted base is N
std::pair<std::string, std::string> MergeSegs(const Segments &segs, const std::vector<std::string> seqs,
                                              bool trim_overhang, int qcutoff,
                                              std::vector<std::string>& ori_quals) {
  //seqs should hold the fastq seq for segs. They could be same length but different seqs
  //Generate a consensus template by opening gaps(MSA format)
  //Each read in ungapped space will be converted to gap space and iteratively updates the consensus
  if (segs.empty()) {
    assert(false);
  }
  int ref_most_left;
  const int NUL = 6;
  const std::string consns_templ = GetConsensusTemplate(segs, ref_most_left);
  std::string consns_seq(consns_templ.size(), NUL);
  std::string consns_qual(consns_templ.size(), 33);
  ori_quals.resize(segs.size());
  for (auto & qual : ori_quals) {
    qual = std::string(consns_templ.size(), 33);
  }

  int32_t start = 0;
  std::string dnaseq, qual;
  for (unsigned sid = 0; sid < segs.size(); ++sid) {
    auto const& seg = segs[sid];
    start = seg.PositionWithSClips() - ref_most_left;
    std::tie(dnaseq, qual) = GetGappedSeqAndQual(seg, seqs[sid], start, consns_templ);
    for (unsigned ii = 0; ii < consns_templ.size(); ++ii) {
      if (trim_overhang && dnaseq[ii] == '.') {
        consns_seq[ii] = '+'; // to be removed
      }
      if (dnaseq[ii] == '.') continue;
      if (consns_seq[ii] == NUL) {
        consns_seq[ii] = dnaseq[ii];
      } else if (consns_seq[ii] != 'N' && consns_seq[ii] != '+') {
         if (dnaseq[ii] != consns_seq[ii]) {
           consns_seq[ii] = dnaseq[ii] == '+' ? '+' : 'N';
         }
      }
      consns_qual[ii] = std::max(qual[ii], consns_qual[ii]);
      ori_quals[sid][ii] = qual[ii];
    }
  }
  for (unsigned ii = 0; ii < consns_templ.size(); ++ii) {
    if (consns_seq[ii] == NUL || consns_seq[ii] == '+' || consns_seq[ii] == '-') {
      consns_seq[ii] = NUL;
      consns_qual[ii] = NUL;
      for (auto& qual : ori_quals) {
        qual[ii] = NUL;
      }
    }
  }
  consns_seq.erase(std::remove(consns_seq.begin(), consns_seq.end(), NUL), consns_seq.end());
  consns_qual.erase(std::remove(consns_qual.begin(), consns_qual.end(), NUL), consns_qual.end());
  for (auto& qual : ori_quals) {
    qual.erase(std::remove(qual.begin(),qual.end(), NUL ), qual.end());
  }
  return std::make_pair(consns_seq, consns_qual);
}

}