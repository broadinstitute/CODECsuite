//
// Created by Ruolin Liu on 9/18/20.
//

#ifndef ADAPTERTRIM_INCLUDE_INDEX_H_
#define ADAPTERTRIM_INCLUDE_INDEX_H_
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <limits>
#include <SeqLib/ssw_cpp.h>
#include <SeqLib/BWAWrapper.h>
#include <SeqLib/BamRecord.h>
#include <memory>
#include <set>
#include "BamRecordExt.h"
#include "StringUtils.h"
#include "FastxIO.h"
#include "Gotoh.h"

namespace CDS {
int SSW(const std::string &ref, const std::string &query) {
  StripedSmithWaterman::Aligner aligner;
  StripedSmithWaterman::Filter filter;
  StripedSmithWaterman::Alignment alignment;
  aligner.Align(query.c_str(), ref.c_str(), ref.size(), filter, &alignment, 15);
  return alignment.mismatches + alignment.query_begin + query.length() - alignment.query_end;
}


bool IsIntermol(const SeqLib::BWAWrapper& bwa,
    const std::string& name1, const std::string& seq1, const std::string& name2, const std::string& seq2) {
  SeqLib::BamRecordVector read1_bam;
  SeqLib::BamRecordVector read2_bam;
  bwa.AlignSequence(seq1, name1, read1_bam, false, -1, 0);
  bwa.AlignSequence(seq2, name2, read2_bam, false, -1, 0);
  if (read1_bam.empty() and read2_bam.empty()) {
    return false;
  }
  if (read1_bam.empty() or read2_bam.empty()) {
    return true;
  }
//  if (verbose) {
//    std::cerr << read1_bam[0] << std::endl;
//    std::cerr << read2_bam[0] << std::endl;
//  }
  auto is = cpputil::InsertSize(read1_bam[0], read2_bam[0]);
  if (is > 1000 || is == 0)  return true;
  else return false;
}

class IndexBarcode {
  std::vector<std::string> index1s_;
  std::vector<std::string> index2s_;
  std::vector<std::string> snames_;
  std::vector<uint64_t> nmatched_;
  std::vector<cpputil::FastqWriter> fq1_writers_;
  std::vector<cpputil::FastqWriter> fq2_writers_;
  cpputil::FastqWriter unkfq1_writer_;
  cpputil::FastqWriter unkfq2_writer_;
  cpputil::FastqWriter hopfq1_writer_;
  cpputil::FastqWriter hopfq2_writer_;
  std::ifstream file_;
  SeqLib::BWAWrapper  bwa_;
  int max_ed_;
  bool out_unmatched_;
  bool out_hopped_;
  bool verbose_;

//  static const int INDEX_START = 3;
//  static const int INDEX_LEN = 18;

  int MatchIndex (const std::string& seq, const std::string& qual, const std::vector<std::string>& indexes, int& nm) {
    int lowest_nm = std::numeric_limits<int>::max();
    int second_lowest_nm = std::numeric_limits<int>::max();
    int best_idx = 0;
    for (unsigned i  = 0; i < indexes.size(); ++i) {
      //int s = SSW(seq, indexes[i]);
      AffineGap ag(seq, indexes[i]);
      Alignment align(seq, indexes[i], ag.Path());
      int s = align.NM();
      if (s < lowest_nm) {
        second_lowest_nm = lowest_nm;
        lowest_nm = s;
        best_idx = i;
      } else if (s < second_lowest_nm) {
        second_lowest_nm = s;
      }
    }
    nm = lowest_nm;
    if (lowest_nm <= max_ed_) {
      return best_idx;
    }
    if (lowest_nm == max_ed_ + 1 && second_lowest_nm > max_ed_ + 3) {
      return best_idx;
    }
    return -1;
  }

 public:
  IndexBarcode (const std::string& index_file, const std::string& outprefix, const int max_ed, bool out_unmatched, bool out_hopped, bool v):
        file_(index_file), max_ed_(max_ed), out_unmatched_(out_unmatched), out_hopped_(out_hopped), verbose_(v)
  {
    std::string header;
    std::string line;
    std::getline(file_, header);
    auto colnames = cpputil::split(header, ",");
    if (colnames.size() != 3 || colnames[0] != "SampleName" || colnames[1] != "IndexBarcode1" || colnames[2] != "IndexBarcode2") {
      throw std::runtime_error("Invalid index file\n Format required as three tab-delimited columns with header SampleName\tIndexBarcode1\tIndexBarcode2");
    }
    if (out_unmatched) {
      unkfq1_writer_.open(outprefix + ".unmatched.1.fastq.gz");
      unkfq2_writer_.open(outprefix + ".unmatched.2.fastq.gz");
    }
    if (out_unmatched) {
      hopfq1_writer_.open(outprefix + ".hopped.1.fastq.gz");
      hopfq2_writer_.open(outprefix + ".hopped.2.fastq.gz");
    }
    std::set<std::string> unique_sids;
    while(std::getline(file_, line)) {
      std::cerr << line << std::endl;
      auto fields = cpputil::split(line, ",");
      if (unique_sids.find(fields[0]) == unique_sids.end()) {
        unique_sids.insert(fields[0]);
        snames_.push_back(fields[0]);
        index1s_.push_back(fields[1]);
        index2s_.push_back(fields[2]);
        fq1_writers_.emplace_back(outprefix + "." + fields[0] + ".1.fastq.gz");
        fq2_writers_.emplace_back(outprefix + "." + fields[0] + ".2.fastq.gz");
      } else {
        std::cerr << "Warning: duplicated sample name in library_params. Ignore \"" << line << "\"\n";
      }
    }
    nmatched_.resize(snames_.size(), 0);
    // Print output header
    if (verbose_) {
      std::cout << "id\t"
                   "observed_1\t"
                   "barcode_1\t"
                   "nm1\t"
                   "observed_2\t"
                   "barcode_2\t"
                   "nm2\t"
                   "sample_1\t"
                   "sample_2\t"
                   "matched\t"
                   "conflicted\t"
                   "hopped" << std::endl;
    }
  }

  void LoadBwa(const std::string &refgenome) {
    std::cerr << "loading index " << refgenome << std::endl;
    bwa_.LoadIndex(refgenome);
    std::cerr << "finished load index " << refgenome << std::endl;
  }

  std::pair<std::string, std::string> ExtractIndex(const cpputil::FastxRecord& read, int index_begin, int index_len) {
    std::string s1 = read.seq.substr(index_begin, index_len);
    std::string q1 = read.qual.substr(index_begin, index_len);
    return std::make_pair(s1, q1);
  }

  void DecodePair(const cpputil::FastxRecord& r1, const cpputil::FastxRecord& r2, int index_begin, int index_len) {
    std::string ob1, qual1, ob2, qual2;
    std::tie(ob1, qual1) = ExtractIndex(r1, index_begin, index_len);
    std::tie(ob2, qual2) = ExtractIndex(r2, index_begin, index_len);
    int nm1, nm2;
    int idx1 = MatchIndex(ob1, qual1, index1s_, nm1);
    int idx2 = MatchIndex(ob2, qual2, index2s_, nm2);
    std::string r1b = idx1 == -1 ? "" : index1s_[idx1];
    std::string r2b = idx2 == -1 ? "" : index2s_[idx2];
    std::string r1s = idx1 == -1 ? "" : snames_[idx1];
    std::string r2s = idx2 == -1 ? "" : snames_[idx2];
    std::string match;
    std::string conflict = "0";
    if (idx1 == -1 || idx2 == -1 || idx1 != idx2 ) {
      if (out_unmatched_) {
        unkfq1_writer_.Write(r1.id, r1.seq, r1.qual);
        unkfq2_writer_.Write(r2.id, r2.seq, r2.qual);
      }
      match = "0";
      if (idx1 != -1 && idx2 != -1 && idx1 != idx2) conflict = "1";
    } else {
      size_t stop = r1.id[r1.id.size() - 1] == ':' ? r1.id.size() : r1.id.size() - 1;
      fq1_writers_[idx1].Write(r1.id.substr(0, stop) + r1b, r1.seq, r1.qual);
      fq2_writers_[idx2].Write(r2.id.substr(0, stop) + r2b, r2.seq, r2.qual);
      match = "1";
      ++nmatched_[idx1];
    }
    std::string hopped = "0";
    if (conflict == "1" and !bwa_.IsEmpty()) {
      hopped = IsIntermol(bwa_, r1.name(), r1.seq.substr(index_begin + index_len + 1),
                          r2.name(), r2.seq.substr(index_begin + index_len + 1)) ? "0" : "1";
      if (out_hopped_ && hopped == "1") {
        hopfq1_writer_.Write(r1.id, r1.seq, r1.qual);
        hopfq2_writer_.Write(r2.id, r2.seq, r2.qual);
      }
    }
    if (verbose_) {
      std::cout << r1.name() << "\t" << ob1 << "\t" << r1b << "\t" << nm1 << "\t" << ob2 << "\t" << \
                 r2b << "\t" << nm2 << "\t" << r1s << "\t" << r2s << "\t"
                << match << "\t" << conflict << "\t" << hopped << "\n";
    }
  }
  uint64_t total_matched() const {return std::accumulate(nmatched_.begin(), nmatched_.end(), (uint64_t) 0);}
  decltype(auto) samples() const {return snames_;}
  decltype(auto) nmatched() const {return nmatched_;}
};

}


#endif //ADAPTERTRIM_INCLUDE_INDEX_H_
