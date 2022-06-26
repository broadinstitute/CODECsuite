//
// Created by Ruolin Liu on 2/26/20.
// This file holds old stuff for finding multiple segments (e.g. oligos) in
// reads.
//

#ifndef CPPUTIL_INCLUDE_ALIGNMENT_H_
#define CPPUTIL_INCLUDE_ALIGNMENT_H_

#include "BamRecordExt.h"

namespace cpputil {

typedef std::vector<SeqLib::BamRecord> Segments; // sequenced parts of a paired-end read for example.
                                                 // or group of duplicated reads

// check all BamReocrd has same start and stop
inline bool AreSegsCompleteOverlap(const Segments &segs) {
  if (segs.empty()) return false;
  int32_t s = segs.front().PositionWithSClips();
  int32_t e = segs.front().PositionEndWithSClips();
  for (auto const & seg : segs) {
    if (seg.PositionWithSClips() != s) return false;
    if (seg.PositionEndWithSClips() != e) return false;
  }
  return true;
}

inline std::string GetUid(const Segments &segs, const std::string& umi_tag) {
  std::string uid;
  bool status;
  if (umi_tag.empty()) status = false;
  else {
    status = segs.front().GetZTag(umi_tag, uid);
    if (!status and segs.size() == 2)  {
      status = segs.back().GetZTag(umi_tag, uid);
    }
  }
  if(!status) { // molecular identifier not exist
    std::string rxtag;
    int32_t s, e;
    if (segs.size() == 1) {
      s = segs.front().PositionWithSClips();
      e = segs.front().PositionEndWithSClips();
    } else if(segs.size() == 2){
      if (segs.front().ReverseFlag()) {
        s = segs.front().PositionEndWithSClips();
        e = segs.back().PositionWithSClips();
      } else {
        s = segs.front().PositionWithSClips();
        e = segs.back().PositionEndWithSClips();
      }
    }
    bool rxtag_status = segs.front().GetZTag("RX", rxtag);
    if (rxtag_status) {
      uid = std::to_string(s) + "," + std::to_string(e) + ":" + rxtag;
    } else {
    }
  }
  return uid;
}

inline bool AreSegsCompleteOverlapExcludingSclip(const Segments &segs) {
  if (segs.empty()) return false;
  int32_t s = segs.front().Position();
  int32_t e = segs.front().PositionEnd();
  for (auto const & seg : segs) {
    if (seg.Position() != s) return false;
    if (seg.PositionEnd() != e) return false;
  }
  return true;
}

inline int GetNumOverlapBasesPEAlignment(const Segments & segs, bool FR_only = true) {
  assert(segs.size() == 2);
  assert(segs.front().Qname() == segs.back().Qname());
  if (segs.front().Interchromosomal()) return 0;
  if (segs.front().PairOrientation() != 0 && FR_only) return 0;
  int left = std::max(segs.front().Position(),  segs.back().Position());
  int right = std::min(segs.front().PositionEnd(),  segs.back().PositionEnd());
  int ol = right - left;
  if (ol < 0) {
    ol = 0;
  }
  return ol;
}

inline bool ArePEAlignmentOverlapAtLeastK(const Segments & segs, int k) {
  if (k == -1) return AreSegsCompleteOverlapExcludingSclip(segs);
  int ol = GetNumOverlapBasesPEAlignment(segs);
  if (ol < k) return false;
//  if (segs.front().Position() < segs.back().PositionEnd() && segs.back().Position() < segs.front().PositionEnd()) {
//    return true;
//  }
  return true;
}

inline bool SegmentNotEmpty(const Segments &seg, int dummy) {
  return !seg.empty();
}

}
#endif //CPPUTIL_INCLUDE_ALIGNMENT_H_
