#ifndef CPPUTIL_INCLUDE_VARIANT_H_
#define CPPUTIL_INCLUDE_VARIANT_H_

#include "SeqLib/BamRecord.h"
#include "SeqLib/RefGenome.h"
#include "DNAUtils.h"
#include "BamRecordExt.h"

namespace cpputil {

struct Bases {
  std::string bases;
  std::vector<int8_t> quals;
  Bases() = default;
  Bases(string bb, string qq): bases(bb), quals(bb.size()){
    std::transform(qq.begin(), qq.end(), quals.begin(), [](char c) -> int8_t {return c - 33;});
  }
  bool empty() const {
    return bases.empty();
  }
  bool operator==(const Bases& other) const {
    return this->bases == other.bases;
  }
  bool operator!=(const Bases& other) const {
    return !(*this == other);
  }
  void clear() {
    bases.clear();
    quals.clear();
  }
};

struct Variant {
  std::string contig;
  int32_t contig_start;
  std::string contig_seq;
  std::string read_id;
  int32_t alt_start;
  std::string alt_seq;
  std::string alt_qual;
  int family_size;
  bool first_of_pair;
  bool rev_strand;
  int var_qual;
  int32_t read_count;
  int32_t dist_to_fragend;
  int32_t r1_start, r2_start;
  std::string alt_qual_str;

  Variant(std::string ctg,
          int32_t cstart,
          std::string cseq,
          std::string rid,
          int32_t rstart,
          std::string rseq,
          std::string rqual,
          int fs,
          bool fop,
          bool isrev):

      contig(ctg),
      contig_start(cstart),
      contig_seq(cseq),
      read_id(rid),
      alt_start(rstart),
      alt_seq(rseq),
      alt_qual(rqual),
      family_size(fs),
      first_of_pair(fop),
      rev_strand(isrev),
      var_qual(-1),
      read_count(1),
      dist_to_fragend(-1),
      r1_start(-1),
      r2_start(-1),
      alt_qual_str("")
  {
    if (Type() == "SNV") {
      int minqual = std::numeric_limits<int>::max();
      for(char x : alt_qual) {
        if ((int)x < minqual) minqual = x;
      }
      var_qual = minqual - 33;
      alt_qual_str = std::to_string(var_qual);
    }
  }

//  Variant():
//        contig_start(-1),
//        alt_start(-1),
//        family_size(1),
//        first_of_pair(true),
//        rev_strand(false),
//        var_qual(-1),
//        read_count(0),
//        dist_to_fragend(-1)
//        {}

  std::string Type() const {
    if (alt_seq.size() == contig_seq.size()) {
      return "SNV";
    }
    else if (alt_seq.size() > contig_seq.size()) {
      return "INS";
    } else if (alt_seq.size() < contig_seq.size()) {
      return "DEL";
    } else {
      throw std::runtime_error("Unknown Var Type\n");
    }
  }

  std::string MutType() {
    std::string res;
    if (Type() == "SNV" && not isMNV()) {
      if (contig_seq == "A" or contig_seq == "G") {
        res = complementString(contig_seq) + ">" + complementString(alt_seq);
      } else {
        res = contig_seq + ">" + alt_seq;
      }
    }
    return res;
  }

  bool isIndel() const {
    if (Type() == "INS" or Type() == "DEL") return true;
    else return false;
  }

  bool IndelLen() const {
    if (alt_seq.size() > contig_seq.size())
      return alt_seq.size() - contig_seq.size();
    else
      return contig_seq.size() - alt_seq.size();
  }

  bool isMNV() const {
    if (alt_seq.size() == contig_seq.size() && alt_seq.size() > 1) return true;
    else return false;
  }

  int32_t NM() const {
  // edit distance
    if (alt_seq.size() == 1 && contig_seq.size() == 1) {
      return 1;
    } else {
      return abs(alt_seq.size() - contig_seq.size());
    }
  }

  bool operator==(const Variant& other) const {
    if (other.contig != this->contig) return false;
    if (other.contig_start != this->contig_start) return false;
    if (other.alt_seq != this->alt_seq) return false;
    return true;
  }

  bool operator!=(const Variant& other) const {
    return !(*this == other);
  }

  bool operator<(const Variant& other) const {
    if (other.contig > this->contig) return false;
    else if (other.contig == this->contig) {
      if (other.contig_start > this->contig_start) return false;
      else if (other.contig_start == this->contig_start) {
        if (other.alt_seq >= this->alt_seq) return false;
      }
    }
    return true;
  }
};

template<typename Stream>
Stream &operator<<(Stream &os, const Variant &v) {
  int32_t fs;
  if (v.dist_to_fragend == -1) {
    fs = v.alt_start < 0 ? abs(v.alt_start) : v.alt_start + 1;
  } else {
    fs = v.dist_to_fragend;
  }
  std::string res = v.contig + "\t" +
                    std::to_string(v.contig_start + 1) + "\t" +
                    v.contig_seq + "\t" +
                    v.alt_seq + "\t" +
                    v.Type() + "\t" +
                    std::to_string(fs) + "\t" +
                    v.alt_qual_str + "\t" +
                    std::to_string(v.read_count) + "\t" +
                    v.read_id + "\t" +
                    std::to_string(v.family_size);
      os << res;
  return os;
}

std::vector<Variant> var_atomize(const Variant& var) {
  std::vector<Variant> res;
  if (var.isMNV()) {
    for (unsigned i = 0; i< var.alt_seq.size(); ++ i) {
      res.emplace_back(var.contig,
                       var.contig_start + i,
                       var.contig_seq.substr(i, 1),
                       var.read_id,
                       var.alt_start + i,
                       var.alt_seq.substr(i, 1),
                       var.alt_qual.substr(i, 1),
                       var.family_size,
                       var.first_of_pair,
                       var.rev_strand);
      res.back().read_count = var.read_count;
      res.back().r1_start = var.r1_start;
      res.back().r2_start = var.r2_start;
    }
  } else {
    res.push_back(var);
  }
  return res;
}

Variant squash_vars(const std::vector<Variant>& vars) {
  //squash varaints from the two reads of a read-pair into one
  if (vars.size() == 1) return vars[0];
  assert(vars.size() == 2);
  assert(vars[0].contig == vars[1].contig &&
          vars[0].contig_start == vars[1].contig_start &&
          vars[0].read_id == vars[1].read_id &&
          vars[0].alt_seq == vars[1].alt_seq);
  Variant ret = vars[0];
  ret.read_count = 2;
  ret.alt_qual_str = vars[0].alt_qual_str + "," + vars[1].alt_qual_str;
  ret.dist_to_fragend = std::min(abs(vars[0].alt_start), abs(vars[1].alt_start));
  if (vars[0].first_of_pair) {
    ret.r1_start = vars[0].r1_start;
    ret.r2_start = vars[1].r2_start;
  } else {
    ret.r2_start = vars[0].r2_start;
    ret.r1_start = vars[1].r1_start;
  }
  return ret;
}

std::vector<Variant> breakmnv(const Variant& var, int minq) {
  assert(var.isMNV());
}

std::vector<std::vector<Variant>> snppair_consolidate(const Variant& var1, const Variant& var2, int minq) {
  //var1 and var2 belongs to read1 and read2 respectively
  //filter or keep as a pair
  //break MNV to SNV if need
  assert(var1 == var2);
  assert(var1.Type() == "SNV");
  std::vector<std::vector<Variant>> res;
  if (var1.var_qual >= minq && var2.var_qual >= minq) {
    res.push_back({var1, var2});
  } else {
    if (var1.isMNV()) {
      auto var1_atoms = var_atomize(var1);
      auto var2_atoms = var_atomize(var2);
      for (unsigned i = 0; i < var1.alt_seq.size(); ++i) {
        if (var1_atoms[i].var_qual < minq || var2_atoms[i].var_qual < minq) {
          var1_atoms[i].var_qual = std::min(var1_atoms[i].var_qual, var2_atoms[i].var_qual);
          var2_atoms[i].var_qual = std::min(var1_atoms[i].var_qual, var2_atoms[i].var_qual);
        }
        res.push_back({var1_atoms[i], var2_atoms[i]});
      }
    } else {
      auto cpvar1 = var1;
      auto cpvar2 = var2;
      cpvar1.var_qual = std::min(var1.var_qual, var2.var_qual);
      cpvar2.var_qual = std::min(var1.var_qual, var2.var_qual);
      res.push_back({cpvar1, cpvar2});
    }
  }
  return res;
}

std::vector<Variant> GetVar(const SeqLib::BamRecord &rec, const SeqLib::BamHeader& header,
                                     const SeqLib::RefGenome& refgenome) {
  std::vector<cpputil::Variant> vars;
  if (GetNM(rec) == 0) return vars;
  const string rname = header.IDtoName(rec.ChrID());
  const int32_t refstart = rec.Position();
  const int32_t refend = rec.PositionEnd();
  const int32_t readstart = rec.AlignmentPosition();
  const int32_t readend = rec.AlignmentEndPosition();

  const auto cigar = rec.GetCigar();
  const std::string seq = rec.Sequence();
  if (refstart > refend - 1) {
    std::cerr << rec << std::endl;
    throw std::runtime_error("invalid cigar string");
  }
  std::string refstr = refgenome.QueryRegion(rname, refstart, refend - 1); // QueryRegion use closed interval
  std::string readstr = seq.substr(readstart, readend - readstart);
  std::string qualstr = rec.Qualities().substr(readstart, readend - readstart);

  std::string refgapstr, readgapstr, qualgapstr;
  int refpos = 0, readpos = 0;
  for (auto cit = cigar.begin(); cit != cigar.end(); ++cit) {
    if (cit->Type() == 'M' or cit->Type() == '=' or cit->Type() == 'X') {
      refgapstr = refstr.substr(refpos, cit->Length());
      readgapstr = readstr.substr(readpos, cit->Length());
      qualgapstr = qualstr.substr(readpos, cit->Length());

      int vlen = 0;
      for (int i = 0, previ = -1; i < (int) refgapstr.size(); ++i) {
        if (refgapstr[i] != readgapstr[i] && readgapstr[i] != 'N' &&
            (refgapstr[i] == 'A' || refgapstr[i] == 'T' || refgapstr[i] == 'G' || refgapstr[i] == 'C' )) {
          vlen = previ + 1 == i ? vlen + 1 : 1;
          previ = i;
        } else {
          if (vlen > 0) { // output mismatch
            int ss = i - vlen;
            int rstart = rec.ReverseFlag() ?   readstart + readpos + ss - seq.size(): readstart + readpos + ss;
            int r1_start = 0, r2_start = 0;
            cpputil::Variant v(rname,
                               refstart + refpos + ss,
                               refgapstr.substr(ss, vlen),
                               rec.Qname(),
                               rstart,
                               readgapstr.substr(ss, vlen),
                               qualgapstr.substr(ss, vlen),
                               GetFamilySize(rec),
                               rec.FirstFlag(),
                               rec.ReverseFlag());
            if (rec.FirstFlag()) v.r1_start = readstart + readpos + ss;
            else v.r2_start = readstart + readpos + ss;
            vars.push_back(v);
            vlen = 0;
          }
        }
      }
      if (vlen >0) {
        int ss = refgapstr.size() - vlen;
        int rstart = rec.ReverseFlag() ?   readstart + readpos + ss - seq.size(): readstart + readpos + ss;
        cpputil::Variant v(rname,
                           refstart + refpos + ss,
                           refgapstr.substr(ss, vlen),
                           rec.Qname(),
                           rstart,
                           readgapstr.substr(ss, vlen),
                           qualgapstr.substr(ss, vlen),
                           GetFamilySize(rec),
                           rec.FirstFlag(),
                           rec.ReverseFlag());
        if (rec.FirstFlag()) v.r1_start = readstart + readpos + ss;
        else v.r2_start = readstart + readpos + ss;
        vars.push_back(v);
        vlen = 0;
      }

//      if (cit->Type() == 'X') {
//        cpputil::Variant v(rname,
//                           refstart + refpos,
//                           refgapstr,
//                           rec.Qname(),
//                           rec.ReverseFlag()? readstart + readpos - seq.size() : readstart + readpos,
//                           readgapstr,
//                           qualgapstr,
//                           rec.FirstFlag());
//        vars.push_back(v);
//      }

      refpos += cit->Length();
      readpos += cit->Length();

    } else if (cit->Type() == 'D') {
      string altseq;
      if (refpos > 0) {
        altseq =  std::string(1, readstr[readpos - 1]);
      }
      if (altseq != "N") {
        int rstart = rec.ReverseFlag() ?   readstart + readpos -1 - seq.size(): readstart + readpos - 1;
        if (refpos == 0)  { // del from the begining, not using 1bp anchor
          cpputil::Variant v(rname,
                             refstart + refpos ,
                             refstr.substr(refpos, cit->Length()),
                             rec.Qname(),
                             rstart,
                             altseq,
                             "",
                             GetFamilySize(rec),
                             rec.FirstFlag(),
                             rec.ReverseFlag());
          if (rec.FirstFlag()) v.r1_start = readstart + readpos -1;
          else v.r2_start = readstart + readpos -1;
          vars.push_back(v);

        } else {
          cpputil::Variant v(rname,
                             refstart + refpos - 1,
                             refstr.substr(refpos - 1, cit->Length() + 1),
                             rec.Qname(),
                             rstart,
                             altseq,
                             std::string(1, qualstr[readpos - 1]),
                             GetFamilySize(rec),
                             rec.FirstFlag(),
                             rec.ReverseFlag());
          if (rec.FirstFlag()) v.r1_start = readstart + readpos -1;
          else v.r2_start = readstart + readpos -1;
          vars.push_back(v);
        }
      }
      refpos += cit->Length();
    } else if (cit->Type() == 'I') {
//      refgapstr += std::string(cit->Length(), '-');
//      readgapstr += readstr.substr(readpos, cit->Length());
      std::string insseq = readstr.substr(readpos, cit->Length());
      if (insseq.find('N') == std::string::npos) {
        std::string refseq = refgenome.QueryRegion(rname, refstart + refpos - 1, refstart + refpos - 1);
        std::string altseq = refseq + insseq;
        std::string qualseq = readpos == 0 ? "I" + qualstr.substr(readpos, cit->Length()) : qualstr.substr(readpos - 1,
                                                                                                           cit->Length()
                                                                                                               + 1);
        int rstart = rec.ReverseFlag() ? readstart + readpos - 1 - seq.size() : readstart + readpos - 1;
        cpputil::Variant v(rname,
                           refstart + refpos - 1,
                           refseq,
                           rec.Qname(),
                           rstart,
                           altseq,
                           qualseq,
                           GetFamilySize(rec),
                           rec.FirstFlag(),
                           rec.ReverseFlag());
        if (rec.FirstFlag()) v.r1_start = readstart + readpos -1;
        else v.r2_start = readstart + readpos -1;
        vars.push_back(v);
      }
      readpos += cit->Length();
    }
  }
  return vars;
}

}
#endif
