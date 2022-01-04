//
// Created by Ruolin Liu on 5/9/20.
//

#ifndef CPPUTIL_INCLUDE_TARGETLAYOUT_H_
#define CPPUTIL_INCLUDE_TARGETLAYOUT_H_

#include <vector>
#include <map>
#include <unordered_map>

#include <SeqLib/BamReader.h>
#include <SeqLib/GenomicRegion.h>
#include <SeqLib/GenomicRegionCollection.h>

namespace cpputil{

class TargetLayout {
  size_t idx_ = 0;
  SeqLib::GenomicRegionVector ginvs_;

  bool _Load(const SeqLib::BamHeader& header, const std::string& bed_path) {
    SeqLib::GRC grc;
    bool ret = grc.ReadBED(bed_path, header);
    if (!ret) {
      throw std::runtime_error(bed_path + " cannot be read!");
    }
    ginvs_ = grc.AsGenomicRegionVector();
    std::cerr << "read " << ginvs_.size() << " regions\n";
    return true;
  }

public:
  TargetLayout() = default;
  TargetLayout(const SeqLib::BamHeader& header, const std::string& bed_path) {
    _Load(header, bed_path);
  }

  size_t NumRegion() const {
    return ginvs_.size();
  }

  decltype(auto) operator[] (int i) const{
    return (ginvs_.at(i));
  }

  bool NextRegion(SeqLib::GenomicRegion & gr) {
    if (idx_ < ginvs_.size()) {
      gr = ginvs_[idx_];
      ++idx_;
      return true;
    } else {
      return false;
    }
  }


};

}
#endif //CPPUTIL_INCLUDE_TARGETLAYOUT_H_
