//
// Created by Ruolin Liu on 3/19/20.
//

#ifndef CPPUTIL_INCLUDE_ALIGNMENTCONSENSUS_H_
#define CPPUTIL_INCLUDE_ALIGNMENTCONSENSUS_H_

#include <string>
#include <vector>
#include <map>
#include <htslib/sam.h>
#include "BamRecordExt.h"
#include "Alignment.h"

namespace cpputil {

inline void find_insert_(const SeqLib::Cigar &cigar, int left_cursor, std::map<int, int> &ins_len) {
  for (auto it = cigar.begin(); it != cigar.end(); ++it) {
    if (it->Type() == 'H') {
      continue;
    } else if (it->Type() == 'I') {
      ins_len[left_cursor] = std::max(ins_len[left_cursor], (int) it->Length());
    } else {
      left_cursor += it->Length();
    }
  }
}


std::string GetConsensusTemplate(const Segments& segs, int32_t& ref_most_left);

std::pair<std::string, std::string>
    GetGappedSeqAndQual(const SeqLib::BamRecord &r, const int start, const std::string consensus_template);

std::pair<std::string, std::string> MergePair(const Segments &segs, const std::vector<std::string>& seqs,
                                              bool trim_overhang, int qcutoff, std::vector<std::string>& out_quals);

std::pair<std::string, std::string> MergePairSeq(const Segments &seg, bool trim_overhang, int qcutoff);
}

#endif //CPPUTIL_INCLUDE_ALIGNMENTCONSENSUS_H_
