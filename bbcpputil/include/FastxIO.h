#ifndef FASTA_IO
#define FASTA_IO
#include <string>
#include <memory>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>
#include <seqan/stream/iostream_bgzf.h>
#include <fstream>

#include "FastxRecord.h"

namespace cpputil {

inline bool endswith(std::string const &value, std::string const &ending) {
  if (ending.size() > value.size()) return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

class FastqWriter {
  std::unique_ptr<std::ofstream> ofstream_;
  std::unique_ptr<seqan::bgzf_ostream> seqStream_;
  //std::mutex mtx_;
 public:
  FastqWriter() = default;
  FastqWriter(std::ofstream& f) : seqStream_(std::make_unique<seqan::bgzf_ostream>(f)) {}
  FastqWriter(const std::string& file) {
    open(file);
  }
  void open(std::ofstream& f) {
    seqStream_ = std::make_unique<seqan::bgzf_ostream>(f);
  }
  void open(const std::string& file) {
    ofstream_ = std::make_unique<std::ofstream>(file);
    seqStream_ = std::make_unique<seqan::bgzf_ostream>(*ofstream_);
  }
  void Write(const std::string &id, const std::string &seq, const std::string &qual) {
    seqan::CharString rid = id;
    seqan::CharString rseq = seq;
    seqan::CharString rqual = qual;
    //mtx_.lock();
    seqan::writeRecord(*seqStream_, rid, rseq, rqual, seqan::Fastq());
    //mtx_.unlock();
  }

  void Write(const std::string &id, const std::string &seq) {
    std::string qual = std::string(seq.size(), 'I');
    Write(id, seq, qual);
  }

  void Write(const FastxRecord &fxr) {
    if (fxr.qual.empty()) {
      Write(fxr.id, fxr.seq);
    } else {
      Write(fxr.id, fxr.seq, fxr.qual);
    }

  }
};

class FastxReader {
  seqan::SeqFileIn seqfilein_;
  int ftype_; //0 for fasta, 1 for fastq, 2 for unknown
 public:
  FastxReader(std::string fastx) : seqfilein_(fastx.c_str()) {
    if (endswith(fastx, ".fq") or endswith(fastx, ".fastq") or
        endswith(fastx, ".fq.gz") or endswith(fastx, ".fastq.gz")) {
      ftype_ = 1;
    } else if (endswith(fastx, ".fa") or endswith(fastx, ".fasta") or
        endswith(fastx, ".fa.gz") or endswith(fastx, ".fasta.gz")) {
      ftype_ = 0;
    } else {
      ftype_ = 2;
      throw std::runtime_error("unknown file format " + fastx);
    }
  }

  bool yield(FastxRecord &record) {
    record.cleanup();
    seqan::CharString id;
    seqan::CharString seq;
    seqan::CharString qual;

    if (seqan::atEnd(seqfilein_)) {
      return false;
    }
    if (ftype_ == 0) {
      seqan::readRecord(id, seq, seqfilein_);
    } else {
      seqan::readRecord(id, seq, qual, seqfilein_);
    }
    record = FastxRecord(seqan::toCString(id), seqan::toCString(seq), seqan::toCString(qual));
    return true;
  }
};

}//end namepsace
#endif
