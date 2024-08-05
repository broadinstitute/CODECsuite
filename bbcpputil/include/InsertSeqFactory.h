#ifndef INSERT_SEQ_FACTORY_H
#define INSERT_SEQ_FACTORY_H

#include <map>
#include <set>
#include <SeqLib/BamReader.h>
#include <SeqLib/GenomicRegion.h>
#include "DNAUtils.h"

#include "TargetLayout.h"
#include "Alignment.h"
#include "Insert.h"
#include "BamRecordExt.h"

#include "FastxIO.h"

#ifndef NDEBUG
#include <iostream>
#endif

#ifndef NDEBUG
#  define DEBUG(x) do {std::cerr << x << std::endl;} while(0);
#else
#  define DEBUG(x) do {} while (0)
#endif

using std::map;
using std::string;

namespace cpputil {
class InsertSeqFactory {

  SeqLib::BamReader bam_reader_;
  map<string, InsertSeq> id_inserts_;
  map<string, InsertSeq>::iterator it_;
  int min_mapq_;
  bool load_supp_;
  bool load_secondary_;
  bool load_duplicate_;
  bool load_proper_pair_only_;
  bool clip3_ = false;
  int last_chr_;
  bool finished_ = false;
  int paired_end_library_ = 1; // 0 single end, 1 paired end

  bool passfilter(const SeqLib::BamRecord& b) {
    if (!b.MappedFlag()) return false;
    if (!load_supp_ && b.SupplementaryFlag()) return false;
    if (!load_secondary_ && b.SecondaryFlag()) return false;
    if (!load_duplicate_ && b.DuplicateFlag()) return false;
    if (load_proper_pair_only_ && !b.ProperPair()) return false;
    if (b.MapQuality() < min_mapq_) return false;
    return true;
  }

  std::pair<string, int> loadrecord(int32_t stop_at = INT32_MAX) {
    //load a single record
    SeqLib::BamRecord b;
    string qname;
    int chrid;
    while (true) {
      auto ret = bam_reader_.GetNextRecord(b);
      if (!ret) {
        return std::make_pair(string(), -1);
      }
      if (!b.MappedFlag()) continue;
      if (b.Position() > stop_at) {
        //Stop and this read is not added
        return std::make_pair(b.Qname(), -1);
      }
      if (!passfilter(b)) continue;
      qname = b.Qname();
      chrid = b.ChrID();
      if (clip3_) {
        cpputil::SoftClip3end(b);
      }

      auto itfind = id_inserts_.find(b.Qname());
      if (itfind == id_inserts_.end()) {
        id_inserts_.emplace(b.Qname(), b);
      } else {
        itfind->second.add(b);
      }
      break;
    }
    return std::make_pair(qname, chrid);
  }

  bool yieldfrags(bool (*Selector)(const Segments &, int param),
                   const bool load_unpair,
                   const int pair_min_ol,
                   const std::string uid_tag_name,
                   std::vector<std::vector<Segments>>& ret) {
    if (id_inserts_.empty()) return false;
    it_ = id_inserts_.begin();
    std::map<std::string, std::vector<Segments>> uid_to_segs;
    int dummy = 1; // psudo uid
    while (true) {// Yield a group of read(s), either paired, single, or family until no more reads.
      std::vector<Segments> recs;
      if (paired_end_library_) {
        if (load_unpair) {
          recs = YieldPairAndUnpair();
        } else {
          recs = YieldPair(Selector, pair_min_ol);
        }
      } else {
        recs = YieldSingle();
      }
      if (recs.empty()) break;
      if (not load_duplicate_) {
        for (auto& rec : recs) {
          uid_to_segs[std::to_string(dummy++)].push_back(rec);
        }
      } else {
        for (auto& rec : recs) {
          std::string uid = GetUid(rec, uid_tag_name);
          uid_to_segs[uid].push_back(rec);
        }
      }
    }
    ret.reserve(uid_to_segs.size());
    for (auto& it: uid_to_segs) {
      ret.push_back(it.second);
    }
    id_inserts_.clear();
    return true;
  }

 public:
  //InsertSeqFactory() = delete;
  InsertSeqFactory() = default;

  InsertSeqFactory(const std::string &bam, int mapq, bool load_supp, bool load_sec, bool load_duplicate, bool load_proper_pair_only, bool clip3):
      min_mapq_(mapq),
      load_supp_(load_supp),
      load_secondary_(load_sec),
      load_duplicate_(load_duplicate),
      load_proper_pair_only_(load_proper_pair_only),
      clip3_(clip3),
      last_chr_(-1) {
    bam_reader_.Open(bam);
    SeqLib::BamRecord b;
  };

  bool IsPairEndLib() const {
    return (paired_end_library_);
  }

  bool ReadByRegion(bool (*Selector)(const Segments &, int param),  // selector
                    const SeqLib::GenomicRegion& gr,
                    std::vector<std::vector<Segments>>& ret,
                    int pair_min_ol, // mim overlap between read1 and read2
                    const std::string uid_tag_name,
                    bool load_unpair) {
    //clearance
    if (not ret.empty()) ret.clear();
    if (not id_inserts_.empty()) id_inserts_.clear();
    bool stat = bam_reader_.SetRegion(gr);
    if (not stat) {
      std::cerr << gr << " not found" << std::endl;
      return false;
    }
    std::string readid;
    int chrid;
    while(true) {
      std::tie(readid, chrid) = loadrecord(gr.pos2);
      if (chrid == -1) break;
    }
    bool status = yieldfrags(Selector, load_unpair, pair_min_ol, uid_tag_name, ret);
    return status;
  }

  std::vector<Segments> FetchReadNameSorted(bool load_unpair = false) {
    std::vector<Segments> ret;
    SeqLib::BamRecord b;
    while(true) {
      bool has_read = bam_reader_.GetNextRecord(b);
      if (!has_read) {
        break;
      }
      if (!passfilter(b)) continue;
      if (clip3_) {
        cpputil::SoftClip3end(b);
      }
      auto itfind = id_inserts_.find(b.Qname());
      if (itfind == id_inserts_.end()) {
        for(auto it  : id_inserts_) {
          it.second.Mate();
          for (auto seg: it.second.paired()) {
            ret.push_back(seg);
          }
          if (load_unpair) {
            for (auto bam: it.second.inprogress()) {
              Segments tmp(1, bam);
              ret.push_back(tmp);
            }
          }
        }
        id_inserts_.clear();
        id_inserts_.emplace(b.Qname(), b);
        return ret;
      } else {
        itfind->second.add(b);
      }
    }
    if (!id_inserts_.empty()) {
      for(auto it  : id_inserts_) {
        it.second.Mate();
        for (auto seg: it.second.paired()) {
          ret.push_back(seg);
        }
        if (load_unpair) {
          for (auto bam: it.second.inprogress()) {
            Segments tmp(1, bam);
            ret.push_back(tmp);
          }
        }
      }
      id_inserts_.clear();
    } else {
      finished_ = true;
    }
    return ret;
  }


  std::vector<std::vector<Segments>> ReadByChrom(bool (*Selector)(const Segments &, int param),  // selector
                                                 int pair_min_ol, // mim overlap between read1 and read2
                                                 const std::string uid_tag_name = "",
                                                 bool load_unpair = false) {
    std::vector<std::vector<Segments>> ret;
    if (finished()) return ret;
    while (true) {
      std::string readid;
      int chrid;
      // Load reads in to id_inserts_;
      std::tie(readid, chrid) = loadrecord();
      if (last_chr_ != -1 && chrid != last_chr_) { // first read in a different chromosome (2nd and above)
        InsertSeq save;
        if (!readid.empty()) { // save this read from processing and put it back to id_inserts_ after processing
          save = id_inserts_[readid];
          id_inserts_.erase(readid);
        }
        it_ = id_inserts_.begin();

        yieldfrags(Selector, load_unpair, pair_min_ol, uid_tag_name, ret);

        if (!readid.empty()) { // Put back the saved read
          id_inserts_.emplace(readid, save);
        }
        last_chr_ = chrid;
        break;
      }
      if (chrid == -1) {
        finished_ = true;
        break;
      }
      last_chr_ = chrid;
    }
    return ret;
  }

  bool finished() {
    return finished_;
  }

  decltype(auto) bamheader() const {
    return (bam_reader_.Header());
  }

  //Iterator
  std::vector<Segments> YieldSingle() {
    std::vector<Segments> res;
    if (it_ != id_inserts_.end()) {
      res.push_back(it_->second.inprogress());
      ++it_;
    }
    return res;
  }

  //Iterator
  std::vector<Segments> YieldPairAndUnpair() {
    std::vector<Segments> res;
    for(;it_ != id_inserts_.end();) {
      it_->second.Mate();
      for (auto seg: it_->second.paired()) {

          res.push_back(seg);
      }
      for (auto bam: it_->second.inprogress()) {
        Segments tmp(1, bam);
        res.push_back(tmp);
      }
      ++it_;
      if (!res.empty()) break;
    }
    return res;
  }

  //Iterator
  std::vector<Segments> YieldPair(bool (*Selector)(const Segments &, int k), int k) {
    /*
     * Always yield R1 and R2 in order
     */
    std::vector<Segments> res;
    for (; it_ != id_inserts_.end();) {
      it_->second.Mate();
      bool found = false;
      for (auto seg: it_->second.paired()) {
        if (Selector(seg, k)) {
          res.push_back(seg);
          found = true;
        }
      }
      ++it_;
      if (found) {
        break;
      }
    }
    return res;
  }

  // Only if both reads of FR pair pass the filter
//  std::vector<Segments> YieldFamily(bool (*Selector)(const Segments &)) {
//    std::vector<Segments> forward_reverse_bams;
//    for (; family_it_ != start_end_2_rnames_.end();) {
//      Segments forwards;
//      Segments reverses;
//      for (auto id : family_it_->second) {
//        auto f = id_inserts_[id].forward_segments();
//        auto r = id_inserts_[id].reverse_segments();
//        if (!f.empty() && !r.empty()) {
//          forwards.insert(forwards.end(), f.begin(), f.end());
//          reverses.insert(reverses.end(), r.begin(), r.end());
//        }
//      }
//      ++family_it_;
//      if (Selector(forwards) && Selector(reverses)) {
//        forward_reverse_bams.push_back(forwards);
//        forward_reverse_bams.push_back(reverses);
//        break;
//      }
//    }
//    return forward_reverse_bams;
//  }

};
}

#endif
