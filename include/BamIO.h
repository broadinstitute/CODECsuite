//
// Created by Ruolin Liu on 7/13/20.
//

#ifndef ADAPTERTRIM_INCLUDE_BAMIO_H_
#define ADAPTERTRIM_INCLUDE_BAMIO_H_
#include <string>
#include <SeqLib/BamRecord.h>
#include <SeqLib/BamWriter.h>
#include <SeqLib/BamHeader.h>
#include "BamRecordExt.h"
#include "DNAUtils.h"
#include "FastxRecord.h"

namespace cpputil {

class UnMappedBamWriter {
  SeqLib::BamWriter bam_writer_;
  std::string sample_;
  std::string rgid_;

  const SeqLib::BamRecord CreateUBamRecord(const ExtFastxRecord& fxr, bool first_read) {
    SeqLib::BamRecord out;
    out.init();
    bam1_t* b = out.raw();
    b->core.tid = -1;
    b->core.pos = -1;
    b->core.qual = 0;
    b->core.flag = first_read ? 77: 141;

    // set dumy mate
    b->core.mtid = -1;
    b->core.mpos = -1;
    b->core.isize = 0;

    // allocate all the data
    b->core.l_qname = fxr.name().length() + 1;
    b->core.l_qseq = fxr.seq.length(); //(seq.length()>>1) + seq.length() % 2; // 4-bit encoding
    b->l_data = b->core.l_qname + ((b->core.l_qseq+1)>>1) + (b->core.l_qseq);
    b->data = (uint8_t*)malloc(b->l_data);

    // allocate the qname
    memcpy(b->data, fxr.name().c_str(), fxr.name().length() + 1);

    // allocate the sequence
    uint8_t* m_bases = b->data + b->core.l_qname;

    // TODO move this out of bigger loop
    int slen = fxr.seq.length();
    for (int i = 0; i < slen; ++i) {
      // bad idea but works for now
      uint8_t base = 15;
      if (fxr.seq.at(i) == 'A')
        base = 1;
      else if (fxr.seq.at(i) == 'C')
        base = 2;
      else if (fxr.seq.at(i) == 'G')
        base = 4;
      else if (fxr.seq.at(i) == 'T')
        base = 8;

      m_bases[i >> 1] &= ~(0xF << ((~i & 1) << 2));   ///< zero out previous 4-bit base encoding
      m_bases[i >> 1] |= base << ((~i & 1) << 2);  ///< insert new 4-bit base encoding

    }
    if (!fxr.qual.empty() && fxr.qual.length() != (unsigned) b->core.l_qseq)
      throw std::invalid_argument("New quality score should be same as seq length");

    // length of qual is always same as seq. If empty qual, just set first bit of qual to 0
    if (not fxr.qual.empty()) {
      char * q = strdup(fxr.qual.data());
      for (size_t i = 0; i < fxr.qual.length(); ++i)
        q[i] -= 33;
      memcpy(bam_get_qual(b), q, fxr.qual.length()); // dont copy /0 terminator
      free(q);
    }

    out.AddIntTag("td", fxr.rc_adpt);
    if (not fxr.umi.empty()) {
      out.AddZTag("RX", fxr.umi.seq());
      out.AddZTag("QX", fxr.umi.qual());
    }
    if (not fxr.adap5.empty()) {
      out.AddZTag("s5", fxr.adap5.seq());
      out.AddZTag("q5", fxr.adap5.qual());
    }
    if (not fxr.adap3.empty()) {
      out.AddZTag("s3", fxr.adap3.seq());
      out.AddZTag("q3", fxr.adap3.qual());
    }
    if (not fxr.trim3.empty()) {
      out.AddZTag("sl", fxr.trim3.seq());
      out.AddZTag("ql", fxr.trim3.qual());
    }
    if (not rgid_.empty()) {
      out.AddZTag("RG", rgid_);
    }
    if (not fxr.barcode.empty()) {
      out.AddZTag("bc", fxr.barcode);
    }
    if (fxr.tm != 255) {
      out.AddIntTag("tm", fxr.tm);
    }
    return out;
  }

  const SeqLib::BamRecord CreateUBamRecord(const SeqLib::BamRecord &bam_template, std::string seq, std::string qual, bool single_end) {
    // seq and qual are assumed to be PLUS strand
//    if (not ProperPair(bam_template)) {
//      throw std::runtime_error("not a proper pair");
//    }
    SeqLib::BamRecord out;
    out.init();
    SeqLib::Cigar c;
    out.SetCigar(c);
    if (bam_template.ReverseFlag()) {
      reverse_complement(seq);
      reverse(qual);
    }
    out.SetQname(bam_template.Qname());
    if (single_end) {
      out.raw()->core.flag = 4;
    }
    else {
      out.raw()->core.flag = bam_template.FirstFlag() ? 77: 141;
    }
    out.SetSequence(seq);
    out.SetQualities(qual, 33);
    out.SetChrID(-1);
    out.SetChrIDMate(-1);
    out.SetPosition(-1);
    out.SetPositionMate(-1);
    int32_t cD;
    std::string rg;
    int32_t cM;
    std::string mi;
    std::string rx;
    if (bam_template.GetIntTag("cD", cD)) {
      out.AddIntTag("cD", cD);
    }
    if (bam_template.GetZTag("RG", rg)) {
      out.AddZTag("RG", rg);
    }
    if (bam_template.GetIntTag("cM", cM)) {
      out.AddIntTag("cM", cM);
    }
    if (bam_template.GetZTag("MI", mi)) {
      out.AddZTag("MI", mi);
    }
    if (bam_template.GetZTag("RX", rx)) {
      out.AddZTag("RX", rx);
    }
    return out;
  }

 public:
  UnMappedBamWriter() = default;
  UnMappedBamWriter(std::string path, std::string rgid, std::string sample) : sample_(sample), rgid_(rgid) {
    std::string header_str = "@HD\tVN:1.5\tGO:none\n";
    header_str += "@RG\tID:" + rgid_ + "\tSM:" + sample_ + "\n";
    SeqLib::BamHeader bh(header_str);
    bam_writer_.Open(path);
    bam_writer_.SetHeader(bh);
    bam_writer_.WriteHeader();
  }

  UnMappedBamWriter(std::string path, const SeqLib::BamHeader& tpl) {
    std::istringstream iss(tpl.AsString());
    std::string line;
    std::string newhdr;
    while (std::getline(iss, line, '\n')) {
      if (line.length() == 0 || line.at(0) != '@') break;
      std::string t = line.substr(0, 3);
      if ( t == "@HD" || t == "@RG") {
        newhdr += line +"\n";
      }
    }
    SeqLib::BamHeader bh(newhdr);
    bam_writer_.Open(path);
    bam_writer_.SetHeader(bh);
    bam_writer_.WriteHeader();
  }

  void Open(std::string path, std::string rgid, std::string sample) {
    std::string header_str = "@HD\tVN:1.5\tGO:none\n";
    header_str += "@RG\tID:" + rgid + "\tSM:" + sample + "\n";
    SeqLib::BamHeader bh(header_str);
    bam_writer_.Open(path);
    bam_writer_.SetHeader(bh);
    bam_writer_.WriteHeader();
  }

//  void Init(std::string path, std::string readgroup, std::string sample) {
//    header_str += "@RG\tID:"
//    header_str += readgroup;
//    header_str += "\tSM:";
//    header_str += sample;
//    header_str += "\n";
//
//  }

  ~UnMappedBamWriter() {
    bam_writer_.Close();
  }
  bool IsOpen() {
    return bam_writer_.IsOpen();
  }


  void WriteRecord(const SeqLib::BamRecord & R1, const SeqLib::BamRecord& R2, std::string seq1, std::string seq2, std::string qual1, std::string qual2) {
    // Write paired end records
    // R1, R2 must be strictly First in pair and, Second in pair
    if (not cpputil::ProperPair(R1) || not cpputil::ProperPair(R2)) {
      throw std::runtime_error("not a proper pair");
    }
    auto r1 = CreateUBamRecord(R1, seq1, qual1, false);
    auto r2 = CreateUBamRecord(R2, seq2, qual2, false);
    //std::cout << out << std::endl;
    bool status1 = bam_writer_.WriteRecord(r1);
    bool status2 = bam_writer_.WriteRecord(r2);
    if (not status1 or not status2) {
      std::cerr << "cannot write bam record " << R1.Qname() << std::endl;
    }
  }

  //simply strip off mapping information and output ubam
  //this is used for intermolecular bams
  void WriteRecord(const SeqLib::BamRecord & R1, const SeqLib::BamRecord& R2) {
    auto r1 = CreateUBamRecord(R1, R1.Sequence(), R1.QualitySequence(), false);
    auto r2 = CreateUBamRecord(R2, R2.Sequence(), R2.QualitySequence(), false);
    bool status1 = bam_writer_.WriteRecord(r1);
    bool status2 = bam_writer_.WriteRecord(r2);
    if (not status1 or not status2) {
      std::cerr << "cannot write bam record " << R1.Qname() << std::endl;
    }
  }

  void WriteRecord(const SeqLib::BamRecord & R1, std::string seq, std::string qual) {
    // Write Single end record
    // R1, R2 must be strictly First in pair and, Second in pair
    if (not cpputil::ProperPair(R1)){
      throw std::runtime_error("not a proper pair");
    }
    auto out = CreateUBamRecord(R1, seq, qual, true);
    bool status = bam_writer_.WriteRecord(out);
    if (not status) {
      std::cerr << "cannot write bam record " << R1.Qname() << std::endl;
    }
  }

  void WriteRecord(ExtFastxRecord& fxr, bool first_read) {
    auto out = CreateUBamRecord(fxr, first_read);
    bool status = bam_writer_.WriteRecord(out);
    if (not status) {
      std::cerr << "cannot write bam record " << fxr.broad_id() << std::endl;
    }
  }
};

}
#endif //ADAPTERTRIM_INCLUDE_BAMIO_H_
