//
// Created by Ruolin Liu on 10/12/20.
//

#ifndef ADAPTERTRIM_CPPUTIL_INCLUDE_FASTXRECORD_H_
#define ADAPTERTRIM_CPPUTIL_INCLUDE_FASTXRECORD_H_

#include <SeqLib/BamRecord.h>
#include "StringUtils.h"
#include "DNAUtils.h"

namespace cpputil {

inline std::string broad_name(const std::string &name) {
  auto fields = split(name, ":");
  std::string res = fields[1];
  for (unsigned i = 3; i < fields.size(); ++i) {
    res += ":" + fields[i];
  }
  return res;
}

inline std::pair<std::string, uint64_t> split_instrument_id_from_broad_name(const std::string &broad_name) {
  auto fields = split(broad_name, ":");
  auto instrument_id = fields[0];
  std::string the_rest = "";
  for (unsigned i = 1; i < fields.size(); ++i) {
    the_rest += fields[i];
  }
  return std::make_pair(instrument_id, std::stoull(the_rest));
}

struct FastxRecord {
  std::string id; // full ID
  std::string seq;
  std::string qual;
  size_t name_idx = std::string::npos;
  FastxRecord() = default;
  FastxRecord(std::string i, std::string s, std::string q) : id(i), seq(s), qual(q) {
    name_idx = id.find(' ');
  };

  FastxRecord(const SeqLib::BamRecord &br, bool duplex_umi = false) {
    id = br.Qname();
    seq = br.Sequence();
    qual = br.Qualities();
    if (br.ReverseFlag()) {
      reverse_complement(seq);
      std::reverse(qual.begin(), qual.end());
    }
    std::string umiseq;
    std::string umiqual;
    if (br.GetZTag("RX", umiseq)) {
      bool status = br.GetZTag("QX", umiqual);
      if (!status) {
        umiqual = std::string('I', umiseq.size());
      }
      if (duplex_umi) {
        auto umiseqs = cpputil::split(umiseq, "-");
        auto umi1_qual = umiqual.substr(0, umiseqs[0].size());
        auto umi2_qual = umiqual.substr(umiseqs[0].size()+1, umiseqs[1].size());
        if (br.FirstFlag()) {
          seq = umiseqs[0] + seq;
          qual = umi1_qual + qual;
        } else {
          seq = umiseqs[1] + seq;
          qual = umi2_qual + qual;
        }
      } else {
        seq = umiseq + seq;
        qual = umiqual + qual;
      }
    }
  }

  void update_id_with_umi(const std::string &umi) {
    auto name = id.substr(0, name_idx);
    auto suffix = id.substr(name_idx);
    id = name + "_" + umi + suffix;
    name_idx += umi.size() + 1;
  }

  virtual void cleanup() {
    id.clear();
    seq.clear();
    qual.clear();
    name_idx = std::string::npos;
  }

  bool is_filtered() {
    std::size_t found = id.find_first_of(':', name_idx);
    if (found < id.size() - 1) {
      if (id[found+1] == 'Y') return true;
    }
    return false;
  }

  std::string index_barcode() const {
    std::size_t found = id.find_last_of(':');
    if (found != std::string::npos) {
      return id.substr(found + 1);
    } else {
      return "";
    }
  }

  std::string name() const {
    return (name_idx != std::string::npos ? id.substr(0, name_idx) : id);
  }

  //broad cannonicalized name,e.g.,
  //"D00203:HCY5YBCX3200606:HCY5YBCX3:1:1105:4656:14095"
  //"HCY5YBCX3200606:1:1105:4656:14095"
  // Warning, should only work for illumina fastq read convention

  std::string broad_id() const {
    if (name_idx != std::string::npos) {
      return broad_name(this->name()) + " " + id.substr(name_idx);
    } else {
      return broad_name(this->name());
    }
  }
};

class AnnotatedSeq {
  std::string seq_;
  std::string qual_;
 public:
  AnnotatedSeq() = default;
  AnnotatedSeq(std::string s, std::string q) : seq_(s), qual_(q) {
    assert(s.size() == q.size());
  };
  bool empty() const {
    return seq_.size() == 0;
  }
  decltype(auto) qual() const {
    return (qual_);
  }
  decltype(auto) seq() const {
    return (seq_);
  }
  void cleanup() {
    seq_.clear();
    qual_.clear();
  }
};

struct ExtFastxRecord : public FastxRecord {
  AnnotatedSeq adap5;
  AnnotatedSeq adap3;
  AnnotatedSeq umi;
  AnnotatedSeq trim3;
  std::string barcode;
  int tm = 255; //unsigned
  int rc_adpt = 0;

  void cleanup() override {
    tm = 255;
    rc_adpt = 0;
    FastxRecord::cleanup();
    adap3.cleanup();
    adap5.cleanup();
    umi.cleanup();
    trim3.cleanup();
    barcode.clear();
  }
};

template<typename Stream>
Stream &operator<<(Stream &os, const FastxRecord &fxr) {
  if (fxr.qual.empty()) {
    os << ">" + fxr.id << "\n" << fxr.seq << "\n";
  } else {
    os << "@" + fxr.id << "\n" << fxr.seq << "\n" << "+\n" << fxr.qual << "\n";
  }
  return os;
}

}
#endif //ADAPTERTRIM_CPPUTIL_INCLUDE_FASTXRECORD_H_
