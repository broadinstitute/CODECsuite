//
// Created by Ruolin Liu on 3/26/20.
//

#ifndef CPPUTIL_INCLUDE_MAF_H_
#define CPPUTIL_INCLUDE_MAF_H_
#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <algorithm>
#include "StringUtils.h"

namespace cpputil {
class MAFReader {
  std::ifstream in_;
  std::map<std::string, std::map<int32_t, std::vector<std::string> >> records_;
  bool isopen_ = false;

 public:
  MAFReader() = default;

  MAFReader(const std::string& maf) {
    Open(maf);
  }

  bool IsOpen() const {
    return isopen_;
  }

  bool Open(const std::string& maf) {
    isopen_ = true;
    in_ = std::ifstream(maf);
    string line;
    //get header
    getline(in_, line, '\n');
    if (in_.eof()) {
      return false;
    }
    std::vector<std::string> fields;
    split_by_char(line, '\t', fields);
    for(auto &s : fields) {
      std::transform(s.begin(), s.end(), s.begin(),
      [](unsigned char c) -> unsigned char { return std::toupper(c); });
    }
    int chr_idx =  std::distance(fields.begin(), std::find(fields.begin(), fields.end(), "CHROMOSOME"));
    int start_idx =  std::distance(fields.begin(), std::find(fields.begin(), fields.end(), "START_POSITION"));
    int alt_idx =  std::distance(fields.begin(), std::find(fields.begin(), fields.end(), "TUMOR_SEQ_ALLELE2"));

    while (true) {
      line.clear();
      fields.clear();
      getline(in_, line, '\n');
      if (in_.eof()) {
        break;
      }
      split_by_char(line, '\t', fields);
      auto it = records_.find(fields[chr_idx]);
      if ( it == records_.end()) {
        std::map<int32_t, std::vector<std::string>> key = {{std::stoi(fields[start_idx]), std::vector<std::string>(1, fields[alt_idx])}};
        records_[fields[chr_idx]] = key;
      } else {
        auto it2 = it->second.find(std::stoi(fields[start_idx]));
        if (it2 == it->second.end()) {
          it->second[std::stoi(fields[start_idx])] = std::vector<std::string>(1, fields[alt_idx]);
        } else {
          it2->second.push_back(fields[alt_idx]);
        }
      }
    }
    return true;
  }

  bool var_exist(const std::string& contig, const int32_t pos, std::string alt="") const {
    auto it = records_.find(contig);
    if (it != records_.end()) {
      auto it2 = it->second.find(pos + 1);
      if (it2 != it->second.end()) {
        if (alt.empty()) return true;
        if( std::find(it2->second.begin(), it2->second.end(), alt) != it2->second.end()) {
          return true;
        }
      }
    }
    return false;
  }

  void Print() const {
    for (auto const& rec : records_) {
      for (auto const& pos_alts : rec.second) {
        for (auto const& alt : pos_alts.second) {
          std::cerr << rec.first << "\t" << pos_alts.first << "\t" << alt << std::endl;
        }
      }
    }
  }
};
}

#endif //CPPUTIL_INCLUDE_MAF_H_
