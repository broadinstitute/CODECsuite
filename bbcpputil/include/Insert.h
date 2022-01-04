#ifndef INSERT_SEQ_H
#define INSERT_SEQ_H
#include <map>
//#include <queue>
#include <string>
#include <algorithm>
#include <list>
#include <vector>
#include <htslib/sam.h>
#include <memory>
#include <iostream>
#include <SeqLib/BamRecord.h>
#include "AlignmentConsensus.h"


using std::list;
using std::vector;
using std::pair;
using std::string;

namespace cpputil {

class InsertSeq {

  using BamPointer = std::shared_ptr<bam1_t>;
  using BamRecord = SeqLib::BamRecord;
  //Switch to SeqLib::BamRecord
  typedef Segments::iterator iterator;
  typedef Segments::const_iterator const_iterator;

  Segments inprogress_, invalid_, ambiguous_;
  vector<Segments> paired_; // segment size is two
  string seqid_;

  struct SeqQual {
    string seq;
    string qual;
    SeqQual(string s, string q) : seq(s), qual(q) {}
  };

  bool IsMate_(const bam1_t *bam, const bam1_t *mate) const {
    /* Modified based on Rsamtool src/Template.h
    //https://github.com/Bioconductor/Rsamtools

    // is_mate checks the following bit flags:
    // 1. Bit 0x40 and 0x80: Segments are a pair of first/last OR
    //    neither segment is marked first/last
    // 2. Bit 0x100: Both segments are secondary OR both not secondary
    // 3. Bit 0x10 and 0x20: Strand flag 0x20 of one mate must match strand
    //                       flag 0x10 of the other mate and vice versa
    // 4. Bit 0x2: Both proper OR both not proper 
    // 5. mpos match:
    //      bit 0x10 of rec1 == bit 0x20 of rec2 AND
    //      bit 0x10 or rec2 == bit 0x20 of rec1
    //      segment2 mpos matches segment1 pos
    // 6. Both Mapped;
    */
    const bool bam_read1 = bam->core.flag & BAM_FREAD1;
    const bool mate_read1 = mate->core.flag & BAM_FREAD1;
    const bool bam_read2 = bam->core.flag & BAM_FREAD2;
    const bool mate_read2 = mate->core.flag & BAM_FREAD2;
    const bool bam_secondary = bam->core.flag & BAM_FSECONDARY;
    const bool mate_secondary = mate->core.flag & BAM_FSECONDARY;
    const bool bam_proper = bam->core.flag & BAM_FPROPER_PAIR;
    const bool mate_proper = mate->core.flag & BAM_FPROPER_PAIR;
    const bool bam_rev = bam->core.flag & BAM_FREVERSE;
    const bool mate_rev = mate->core.flag & BAM_FREVERSE;
    const bool bam_mrev = bam->core.flag & BAM_FMREVERSE;
    const bool mate_mrev = mate->core.flag & BAM_FMREVERSE;
    const bool bam_unmap = bam->core.flag & BAM_FUNMAP;
    const bool mate_unmap = bam->core.flag & BAM_FMUNMAP;
    const uint32_t
        pos = bam->core.pos,
        mpos = bam->core.mpos,
        mate_pos = mate->core.pos,
        mate_mpos = mate->core.mpos;
    return
        ((bam_read1 ^ bam_read2) && (mate_read1 ^ mate_read2)) &&
            (bam_read1 != mate_read1) &&
            (bam_secondary == mate_secondary) &&
            (((bam_rev != mate_mrev) && (bam_mrev != mate_rev)) ||
                ((bam_rev == mate_mrev) && (bam_mrev == mate_rev))) &&
            (bam_proper == mate_proper) &&
            !bam_unmap &&
            !mate_unmap &&
            (pos == mate_mpos) && (mpos == mate_pos) &&
            (bam->core.mtid == mate->core.tid);
  }

  void add_to_paired(BamRecord bam, BamRecord mate) {
    // keep the order of R1 then R2
    if (bam.FirstFlag()) {
      Segments tmp{bam, mate};
      paired_.emplace_back(tmp);
    } else {
      Segments tmp{mate, bam};
      paired_.emplace_back(tmp);
    }
  }

  static int32_t GetBreakPointCorrection_(const int32_t stop, const SeqLib::Cigar &cigar) {
    int32_t num_refbase_consumed = 0;
    int32_t correction = 0;
    assert(std::distance(cigar.begin(), cigar.end()) > 1);
    for (auto it = cigar.begin(); it != cigar.end(); ++it) {
      if (num_refbase_consumed > stop) break;
      if (it->Type() == 'I') {
        correction += it->Length();
      } else if (it->Type() == 'D') {
        num_refbase_consumed += it->Length();
        correction -= it->Length();
      } else if (it->Type() == 'H') {
        continue;
      } else {
        num_refbase_consumed += it->Length();
      }
    }
    return correction;
  }

  static pair<SeqQual, SeqQual> GetSegmentOverhang(const Segments &seg) {
    int32_t left_front = seg.front().PositionWithSClips();
    int32_t left_end = seg.front().PositionEndWithSClips();
    int32_t right_front = seg.back().PositionWithSClips();
    int32_t right_end = seg.back().PositionEndWithSClips();
//    if (left_front > right_front) {
//      DEBUG(seg.front())
//      DEBUG(seg.back());
//    }
    assert(left_front <= right_front);
    assert(left_end <= right_end);

    int32_t left_break = right_front - left_front;
    int32_t right_break = left_end - right_front + 1;
    left_break += GetBreakPointCorrection_(left_break, seg.front().GetCigar());
    right_break += GetBreakPointCorrection_(right_break, seg.back().GetCigar());
    SeqQual left_oh(seg.front().Sequence().substr(0, left_break), seg.front().Qualities().substr(0, left_break));
    SeqQual right_oh(seg.back().Sequence().substr(right_break - 1, seg.back().Sequence().size() - right_break + 1),
                     seg.back().Qualities().substr(right_break - 1, seg.back().Qualities().size() - right_break + 1));

    return std::make_pair(left_oh, right_oh);
  }

//  static bool IsSorted(const Segments &seg) {
//    return seg.front().PositionWithSClips() <= seg.back().PositionWithSClips();
//  }


//  static string GetConsensusTemplate(const Segments &seg) {
//    /*
//     *  '~' : uninitialized
//     *  '+' : insertion
//     *  '-' : deletion
//     */
//    assert(IsSorted(seg));
//    const SeqLib::Cigar left_cigar = seg.front().GetCigar();
//    const SeqLib::Cigar right_cigar = seg.back().GetCigar();
//    int ref_span = std::max(seg.back().PositionEndWithSClips(), seg.front().PositionEndWithSClips())
//        - seg.front().PositionWithSClips();
//    std::map<int, int> ins_len;
//    find_insert_(left_cigar, 0, ins_len);
//    find_insert_(right_cigar, seg.back().PositionWithSClips() - seg.front().PositionWithSClips(), ins_len);
//    string consens;
//    int b = 0;
//    for (const auto pos_len : ins_len) {
//      consens += string(pos_len.first - b, '~');
//      consens += string(pos_len.second, '+');
//      b = pos_len.first;
//    }
//    consens += string(ref_span - b, '~');
//    return consens;
//  }


 public:
  InsertSeq() = default;
  InsertSeq(BamRecord br) {
    add(br);
  }

  void add(BamRecord br) {
    inprogress_.push_back(br);
  }

  decltype(auto) paired() const {
    return (paired_);
  }

  decltype(auto) inprogress() const {
    return (inprogress_);
  }

  decltype(auto) forward_segments() const {
    Segments forward;
    for (const auto &read : inprogress_) {
      if (read.ReverseFlag()) {
        continue;
      } else {
        forward.push_back(read);
      }
    }
    return forward;
  }

  decltype(auto) reverse_segments() const {
    Segments reverse;
    for (const auto &read : inprogress_) {
      if (read.ReverseFlag()) {
        reverse.push_back(read);
      }
    }
    return reverse;
  }

  bool empty() const {
    return inprogress_.empty() && invalid_.empty() && ambiguous_.empty() && paired_.empty();
  }

  void Mate() {
    /* Adapted from Rsamtool src/Template.h
    //https://github.com/Bioconductor/Rsamtools
    */
    // This is O(n^2) where n is the number of reads. This does not work for large n.

    // Mate paired bam records to segments. Segments are non-overlap intervals on genome.
    const int unmated = -1, multiple = -2, processed = -3;
    vector<pair<int, BamRecord> >
        status(inprogress_.size(),
               pair<int, BamRecord>(unmated, BamRecord()));
    Segments::iterator it0;

    // identify unambiguous and ambiguous mates
    it0 = inprogress_.begin();
    for (unsigned int i = 0; i < inprogress_.size(); ++i) {
      status[i].second = *it0;
      Segments::iterator it1 = it0;
      for (unsigned int j = i + 1; j < inprogress_.size(); ++j) {
        ++it1;
        if (IsMate_(it0->raw(), it1->raw())) {
          status[i].first = status[i].first == unmated ? j : multiple;
          status[j].first = status[j].first == unmated ? i : multiple;
        }
      }
      ++it0;
    }

    // process unambiguous and ambiguous mates
    for (unsigned int i = 0; i < status.size(); ++i) {
      if (status[i].first == unmated)
        continue;
      if (status[i].first >= 0 && status[status[i].first].first >= 0) {
        // unambiguous mates
        add_to_paired(status[i].second, status[status[i].first].second);
        status[status[i].first].first = processed;
        status[i].first = processed;
      } else if (status[i].first != processed) {
        // ambiguous mates, added to 'ambiguous' queue
        ambiguous_.push_back(status[i].second);
        status[i].first = processed;
      }
      ++it0;
    }

    // remove segments that have been assigned to paired or
    // ambiguous queue
    it0 = inprogress_.begin();
    for (unsigned int i = 0; i != status.size(); ++i) {
      if (status[i].first == processed) {
        it0 = inprogress_.erase(it0);
      } else {
        ++it0;
      }
    }
  }
};
}
#endif
