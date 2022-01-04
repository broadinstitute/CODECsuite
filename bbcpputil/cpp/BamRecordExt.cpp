//
// Created by Ruolin Liu on 3/3/20.
//

#include <string>
#include <BamRecordExt.h>
#include "Algo.h"

// Get the number of alignment N bases. Excluding soft clipping.
namespace cpputil {

bool ProperPair(const SeqLib::BamRecord& bam) {
  if (bam.Interchromosomal()) return false;
  return bam.PairOrientation() == 0 ;
}

std::pair<int32_t, int32_t> MatePositionAndPositionEndWithSoftClip(const SeqLib::BamRecord & bam) {
  int32_t mp = bam.MatePosition();
  int32_t mate_position_end_with_sclip = mp;
  std::string  mcigar_str;
  auto stat = bam.GetZTag("MC", mcigar_str);
  if (!stat) {
    std::cerr << bam.Qname() + " does not have MC tag; consider run with the -m/--mate_pair_missing option.\n";
    return std::make_pair(0, 0);
  }
  SeqLib::Cigar mcigar(mcigar_str);
  if (mcigar.front().Type() == 'S')
    mp -= mcigar.front().Length();
  for (auto c = mcigar.begin(); c != mcigar.end(); ++c) {
    if (c->ConsumesReference()) {
      mate_position_end_with_sclip += c->Length();
    }
  }
  if (mcigar.back().Type() == 'S') {
    mate_position_end_with_sclip += mcigar.back().Length();
  }
  return std::make_pair(mp, mate_position_end_with_sclip);
}

int32_t GetUnclippedFramgentLength(const SeqLib::BamRecord &b) {
  if (!b.PairedFlag()) {
    return b.Sequence().size();
  } else {
    int32_t fl = abs(b.InsertSize());
    if (fl == 0 || b.PairOrientation() != 0) {
      return b.Sequence().size();
    }
    std::string  mcigar_str;
    auto stat = b.GetZTag("MC", mcigar_str);
    if (!stat) { // return insert size if no MC tag exist
      return fl;
    } else {
      int32_t mbegin, mend;
      std::tie(mbegin, mend) = MatePositionAndPositionEndWithSoftClip(b);
      if (b.ReverseFlag()) {
        return b.PositionEndWithSClips() - mbegin;
      } else {
        return mend - b.PositionWithSClips();
      }
    }
  }
}


std::pair<int32_t, int32_t> MatePositionAndPositionEnd(const SeqLib::BamRecord & bam) {
  int32_t mp = bam.MatePosition();
  int32_t mate_position_end = mp;
  std::string  mcigar_str;
  bam.GetZTag("MC", mcigar_str);
  SeqLib::Cigar mcigar(mcigar_str);
  for (auto c = mcigar.begin(); c != mcigar.end(); ++c) {
    if (c->ConsumesReference()) {
      mate_position_end += c->Length();
    }
  }
  return std::make_pair(mp, mate_position_end);
}

//overlap len of paired end reads in the reference coordinate, excluding soft clipping
int32_t InsertSize(const SeqLib::BamRecord & read1, const SeqLib::BamRecord& read2) {
  assert(read1.Qname() == read2.Qname());
  if (read1.ChrID() != read2.ChrID()) return 0;
  if (read1.ReverseFlag() ^ read2.ReverseFlag() == 0) return 0;
  if (not read1.MappedFlag() || not read2.MappedFlag()) return 0;
  int right, left;
  if (read1.ReverseFlag()) {
    right = read1.PositionEnd();
    left = read2.Position();
  } else {
    left = read1.Position();
    right = read2.PositionEnd();
  }
  if (right < left) return 0;
  return right - left;
}

std::pair<uint64_t, uint64_t> ProperPairFramgentEndsWithSclip(const SeqLib::BamRecord &b) {
  if (b.PairOrientation() == 0) { // proper FR pair read
    int32_t mbegin, mend;
    std::tie(mbegin, mend) = MatePositionAndPositionEndWithSoftClip(b);
    if (b.ReverseFlag()) {
      return std::make_pair( (uint64_t) b.MateChrID() << 32 | (uint32_t) mbegin,
                             (uint64_t) b.ChrID() << 32 | (uint32_t) b.PositionEndWithSClips());
    } else {
      return std::make_pair( (uint64_t) b.ChrID() << 32 | (uint32_t) b.PositionWithSClips(),
                             (uint64_t) b.MateChrID() << 32 | (uint32_t) mend);
    }
  }
  return std::make_pair((uint64_t)0, (uint64_t)0);
}

void PrintQual(const SeqLib::BamRecord &b) {
  const std::string qual = b.Qualities();
  const std::string seq = b.Sequence();
  std::string display1, display2;
  std::string marker(qual.size(), ' ');
  for (unsigned ii = 0; ii < qual.size(); ++ii) {
    int q = (int) qual[ii] - 33;
    if (q < 30) marker[ii] = '*';
    int tens = q / 10;
    int units = q % 10;
    display1 += std::to_string(tens);
    display2 += std::to_string(units);
  }
  std::cerr << seq << std::endl;
  std::cerr << display1 << std::endl;
  std::cerr << display2 << std::endl;
  std::cerr << marker << std::endl;
}

//int32_t CountNOrLowQInMatchedBases(const SeqLib::BamRecord &b, const int qcutoff) {
//  //Valid bases are not N and baseq >= qcutoff
//  uint32_t *c = bam_get_cigar(b.raw());
//  int32_t readpos = 0; //left clipping
//  size_t i = 0;
//  for (; i < b.raw()->core.n_cigar; ++i) {
//    if (bam_cigar_opchr(c[i]) == 'S')
//      readpos += bam_cigar_oplen(c[i]);
//    else if (bam_cigar_opchr(c[i]) != 'H')
//      break;
//  }
//
//  int n = 0;
//  uint8_t *p = bam_get_seq(b.raw());
//  uint8_t *bq = bam_get_qual(b.raw());
//
//  for (; i < b.raw()->core.n_cigar; ++i) {
//    if (bam_cigar_opchr(c[i]) == 'M') {
//      for(int ww = 0; ww < bam_cigar_oplen(c[i]); ++ww) {
//        if (bam_seqi(p, ww + readpos) == 15 || bq[ww + readpos] < qcutoff)
//          ++n;
//      }
//      readpos += bam_cigar_oplen(c[i]);
//    } else if(bam_cigar_opchr(c[i]) == 'I') {
//      readpos += bam_cigar_oplen(c[i]);
//    }
//  }
//  return n;
//}

void AddMatchedBasesToCycleCount( const SeqLib::BamRecord& b,
    std::vector<int64_t>& q0_cycle_count,
    std::vector<int64_t>& q30_cycle_count,
    int start,
    int end) {
  //Valid bases are not N and baseq >= qcutoff
  if (end <= start) return;
  if (end > b.AlignmentEndPosition()) {
    end = b.AlignmentEndPosition();
  }
  if (start < b.AlignmentPosition()) {
    start = b.AlignmentPosition();
  }
  uint32_t *c = bam_get_cigar(b.raw());
  int32_t readpos = 0; //left clipping
  size_t i = 0;
  for (; i < b.raw()->core.n_cigar; ++i) {
    if (bam_cigar_opchr(c[i]) == 'S')
      readpos += bam_cigar_oplen(c[i]);
    else if (bam_cigar_opchr(c[i]) != 'H')
      break;
  }

  const uint8_t *p = bam_get_seq(b.raw());
  const uint8_t *bq = bam_get_qual(b.raw());
  const int32_t rl = b.raw()->core.l_qseq;
  const bool is_reverse = b.ReverseFlag();
  for (; i < b.raw()->core.n_cigar; ++i) {
    char cigar = bam_cigar_opchr(c[i]);
    if (cigar == 'M' or cigar == 'X' or cigar == '=') {
      for(int ww = 0; ww < bam_cigar_oplen(c[i]); ++ww) {
        if (ww + readpos <  start) continue;
        if (ww + readpos >= end) break;
        if (bam_seqi(p, ww + readpos) != 15) {
          int cycle = is_reverse ? rl - readpos - ww - 1: readpos + ww;
          ++q0_cycle_count[cycle];
          if (bq[ww + readpos] >= 30)
            ++q30_cycle_count[cycle];
        }
      }
      readpos += bam_cigar_oplen(c[i]);
    } else if(bam_cigar_opchr(c[i]) == 'I') {
      readpos += bam_cigar_oplen(c[i]);
    }
  }
  return;
}

int GetFamilySize(const SeqLib::BamRecord &bam) {
  uint8_t *p1 = bam_aux_get(bam.raw(), "cD");
  if (!p1) {
    return 1; // cD tag not exist
  }
  int32_t nm = bam_aux2i(p1);
  return nm;
}

int32_t CountNBasesInAlignment(const SeqLib::BamRecord &b) {
  uint32_t *c = bam_get_cigar(b.raw());
  int32_t lc = 0; //left clipping
  for (size_t i = 0; i < b.raw()->core.n_cigar; ++i) {
    if (bam_cigar_opchr(c[i]) == 'S')
      lc += bam_cigar_oplen(c[i]);
    else if (bam_cigar_opchr(c[i]) != 'H')
      break;
  }

  int32_t rc = 0; // right clipping
  for (int32_t i = b.raw()->core.n_cigar - 1; i >= 0; --i) { // loop from the end
    if ((bam_cigar_opchr(c[i]) == 'S'))
      rc += bam_cigar_oplen(c[i]);
    else if (bam_cigar_opchr(c[i]) != 'H')// not a clip, so stop counting
      break;
  }

  uint8_t *p = bam_get_seq(b.raw());
  int32_t n = 0;
  for (int ww = lc; ww < b.raw()->core.l_qseq - rc; ww++)
    if (bam_seqi(p, ww) == 15)
      ++n;
  return n;
}

int32_t GetTotalIndelLen(const SeqLib::BamRecord &bam) {
  uint32_t *c = bam_get_cigar(bam.raw());
  int32_t l = 0;
  for (size_t i = 0; i < bam.raw()->core.n_cigar; ++i) {
    if (bam_cigar_opchr(c[i]) == 'D' or bam_cigar_opchr(c[i]) == 'I')
      l += bam_cigar_oplen(c[i]);
  }
  return l;
}

int32_t GetNM(const SeqLib::BamRecord &bam) {
  uint8_t *p1 = bam_aux_get(bam.raw(), "NM");
  if (!p1) {
    return -1; // NM tag not exist
  }
  int32_t nm = bam_aux2i(p1);
  return nm;
}

int32_t GetNMismatchX(const SeqLib::BamRecord &bam) {
  //only support "=/X" type of cigar now
  uint32_t *c = bam_get_cigar(bam.raw());
  bool newtype_cigar = false;
  int nm = 0;
  for (int i = 0; i < bam.raw()->core.n_cigar; ++i) {
    char cigar = bam_cigar_opchr(c[i]);
    if ( cigar == '=') {
      newtype_cigar = true;
    } else if (cigar == 'X') {
      nm += bam_cigar_oplen(c[i]);
    }
  }
  if (!newtype_cigar) return -1;
  else return nm;
}

int32_t IndelLen(const SeqLib::BamRecord &bam) {
  uint32_t *c = bam_get_cigar(bam.raw());
  int il = 0;
  for (int i = 0; i < bam.raw()->core.n_cigar; ++i) {
    char cigar = bam_cigar_opchr(c[i]);
    if ( cigar == 'I' or cigar == 'D') {
      il += bam_cigar_oplen(c[i]);
    }
  }
  return il;
}

int32_t GetNMismatch(const SeqLib::BamRecord &bam, bool NisMM) {
  //only works if MD tag exists
  std::string  mdstr;
  int nm = 0;
  auto status = bam.GetZTag("MD", mdstr);
  if (!status) {
    std::cerr << "MD tag not exist. use NM " << bam.Qname() << std::endl;
    int nm = 0;
    bam.GetIntTag("NM", nm);
    return nm;
    //throw std::runtime_error("MD tag not exist");
  }
  bool del = false;
  for (size_t i = 0; i < mdstr.size(); ++i) {
    if (mdstr[i] == '^') {
      del = true;
      continue;
    }
    if (del) {
      del = mdstr[i] < 65 ? false : true;
    } else {
      if (mdstr[i] >= 65 and mdstr[i] <= 90) { // A-Z
        ++nm;
      }
    }
  }
  if (not NisMM) {
    int numN = CountNBasesInAlignment(bam);
    nm -= numN;
  }
  return nm;
}


bool HasClusteredMuts(const SeqLib::BamRecord &rec, const SeqLib::BamHeader& header,
                      const SeqLib::RefGenome& refgenome, const int cutoff) {
  /*
   * Filter clustered mutations near the end of alignment.
   * if nmut within a window = cutoff * dist  is equal or larger than cutoff,
   */
  const int dist = 15;
  if (GetNMismatch(rec) < cutoff) return false;
  const std::string rname = header.IDtoName(rec.ChrID());
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
  std::vector<int> mutpos;
  int refpos = 0, readpos = 0;
  for (auto cit = cigar.begin(); cit != cigar.end(); ++cit) {
    if (cit->Type() == 'M' or cit->Type() == '=' or cit->Type() == 'X') {
      refgapstr = refstr.substr(refpos, cit->Length());
      readgapstr = readstr.substr(readpos, cit->Length());
      qualgapstr = qualstr.substr(readpos, cit->Length());

      for (int i = 0; i < (int) refgapstr.size(); ++i) {
        if (refgapstr[i] != readgapstr[i] && readgapstr[i] != 'N' &&
            (refgapstr[i] == 'A' || refgapstr[i] == 'T' || refgapstr[i] == 'G' || refgapstr[i] == 'C' )) {
          mutpos.push_back(readstart + readpos + i);
        }
      }

      refpos += cit->Length();
      readpos += cit->Length();

    } else if (cit->Type() == 'D') {
      refpos += cit->Length();
    } else if (cit->Type() == 'I') {
      readpos += cit->Length();
    }
  }
  int beg, end, n;
  n = largest_cluster(mutpos, cutoff * dist, beg, end);
  if (n >= cutoff) return true;
  else return false;
}

int32_t NumSoftClip5End(const SeqLib::BamRecord &bam) {
  int32_t p = 0;
  uint32_t* c = bam_get_cigar(bam.raw());
  if (bam.ReverseFlag()) {
    for (int i = bam.raw()->core.n_cigar - 1; i >=0; --i) {
      if (bam_cigar_opchr(c[i]) == 'S' ) {
        p = bam_cigar_oplen(c[i]);
      } else if (bam_cigar_opchr(c[i]) != 'H') {
        break;
      }
    }
  } else {
    for (int i = 0; i < bam.raw()->core.n_cigar ; ++i) {
      if (bam_cigar_opchr(c[i]) == 'S') {
        p = bam_cigar_oplen(c[i]);
      }else if (bam_cigar_opchr(c[i]) != 'H') {
        break;
      }
    }
  }
  return p;
}

uint32_t GetNumNonIndelAlignedBases(const SeqLib::BamRecord &bam) {
  uint32_t* c = bam_get_cigar(bam.raw());
  uint32_t dmax = 0;
  for (size_t i = 0; i < bam.raw()->core.n_cigar; i++) {
    char cigar = bam_cigar_opchr(c[i]);
    if (cigar == 'M' or cigar == 'X' or cigar == '=')
      dmax += bam_cigar_oplen(c[i]);
  }
  return dmax;
}

int SoftClip3end(SeqLib::BamRecord &bam) {
  auto cigar = bam.GetCigar();
  if (cigar.size() == 0) return 0;
  if (bam.ReverseFlag()) {
    if (cigar.front().Type() == 'S') {
      int front_trim = cigar.front().Length();
      std::string newseq = bam.Sequence().substr(front_trim, bam.Sequence().size() - front_trim);
      std::string newqual = bam.Qualities().substr(front_trim, bam.Qualities().size() - front_trim);
      SeqLib::Cigar newc;
      for (unsigned ii = 1; ii < cigar.size(); ++ii) {
        newc.add(cigar[ii]);
      }
      bam.SetCigar(newc);
      bam.SetSequence(newseq);
      bam.SetQualities(newqual, 33);
      return front_trim;
    }
  } else {
    if (cigar.back().Type() == 'S') {
      int back_trim = cigar.back().Length();

      std::string newseq = bam.Sequence().substr(0, bam.Sequence().size() -  back_trim);
      std::string newqual = bam.Qualities().substr(0, bam.Qualities().size() - back_trim);

      SeqLib::Cigar newc;
      for (unsigned ii = 0; ii < cigar.size() - 1; ++ii) {
        newc.add(cigar[ii]);
      }
      bam.SetCigar(newc);
      bam.SetSequence(newseq);
      bam.SetQualities(newqual, 33);
      return back_trim;
    }
  }
  return 0;
}

bool SoftClipBamRecord(SeqLib::BamRecord &bam) {
  auto cigar = bam.GetCigar();
  if (cigar.front().Type() != 'S' && cigar.back().Type() != 'S') {
    return false;
  }
  SeqLib::Cigar newc;
  for (auto c = cigar.begin(); c != cigar.end(); ++c) {
    if (c->Type() != 'S') {
      newc.add(*c);
    }
  }
  int front_trim = 0;
  int back_trim = 0;
  if (cigar.front().Type() == 'S') {
    front_trim = cigar.front().Length();
  }
  if (cigar.back().Type() == 'S') {
    back_trim = cigar.back().Length();
  }
  std::string newseq = bam.Sequence().substr(front_trim, bam.Sequence().size() - front_trim - back_trim);
  std::string newqual = bam.Qualities().substr(front_trim, bam.Qualities().size() - front_trim - back_trim);

  bam.SetCigar(newc);
  bam.SetSequence(newseq);
  bam.SetQualities(newqual, 33);
  return true;
}

void MaskBaseBelowMinBq(SeqLib::BamRecord &bam, int32_t mbp) {
  //TODO: update the NM tag
  //TODO: even better the MD tag
  if (mbp == 0) return;
  const int offset = 33;
  std::string quals = bam.Qualities();
  std::string dnas = bam.Sequence();
  for (unsigned i = 0; i < quals.size(); ++i) {
    if ((int) quals[i] < mbp + offset) {
      dnas[i] = 'N';
    }
  }
  bam.SetSequence(dnas);
  bam.SetQualities(quals, offset);
}

void TrimBamFromFragEnd(SeqLib::BamRecord &bam, int32_t mp, int32_t mate_position_end_with_sclip, int32_t end5, int32_t end3) {
  if (end5 == 0 and end3 == 0) return; // no trim
  std::string seq = bam.Sequence();
  std::string qual = bam.Qualities();
  //if (bam.SupplementaryFlag()) return;
  int32_t replace5, replace3;
  if (!bam.Interchromosomal()
      && bam.PairOrientation() == 0
      && mp != 0
      && mate_position_end_with_sclip != 0)
  { // only FR paired read
    if (bam.ReverseFlag()) {
      replace5  = std::max(mp + (int32_t) end5 - bam.PositionWithSClips(), 0);
      replace5 = std::min(replace5, (int32_t) seq.size());
      replace3 = std::min(end3, (int32_t) seq.size());
    } else {
      replace3 = std::max(bam.PositionEndWithSClips() - mate_position_end_with_sclip + end3, 0);
      replace3 = std::min((int32_t) seq.size(), replace3);
      replace5 = std::min(end5, (int32_t) seq.size());
    }
  } else {
    replace5 = std::min(end5, (int32_t) seq.size());
    replace3 = std::min(end3, (int32_t) seq.size());
  }

  seq.replace(0, replace5, std::string(replace5, 'N'));
  seq.replace(seq.size() - replace3, replace3, std::string(replace3, 'N'));
  bam.SetSequence(seq);
  bam.SetQualities(qual, 33);
}

void TrimPairFromFragEnd(SeqLib::BamRecord &left, SeqLib::BamRecord&right, int32_t n_trim) {
  if (n_trim <= 0) return;
  int32_t mp = right.PositionWithSClips();
  int32_t mep = right.PositionEndWithSClips();
  TrimBamFromFragEnd(left, mp, mep, n_trim, n_trim);
  mp = left.PositionWithSClips();
  mep = left.PositionEndWithSClips();
  TrimBamFromFragEnd(right, mp, mep, n_trim, n_trim);
}

void TrimSingleFromFragEnd(SeqLib::BamRecord &bam, int32_t n_trim) {
  if (n_trim <= 0) return;
  int32_t mp = bam.PositionWithSClips();
  int32_t mep = bam.PositionEndWithSClips();
  TrimBamFromFragEnd(bam, mp, mep, n_trim, n_trim);
}

int RefPosToQueryPos(const SeqLib::BamRecord &bam, const int refpos) {
  if (refpos < bam.PositionWithSClips() || refpos >= bam.PositionEndWithSClips()) {
    return -1;
  }
  if (refpos == bam.PositionEndWithSClips()) return 0;
  int32_t refscan = bam.PositionWithSClips();
  int32_t readscan = 0;
  uint32_t *c = bam_get_cigar(bam.raw());
  for (int i = 0; i < bam.raw()->core.n_cigar; ++i) {
    char cigar = bam_cigar_opchr(c[i]);
    if ( cigar == 'M' || cigar  == 'S' || cigar == 'X' || cigar == '=') {
      refscan += bam_cigar_oplen(c[i]);
      readscan += bam_cigar_oplen(c[i]);
    } else if(bam_cigar_opchr(c[i]) == 'I') {
      readscan += bam_cigar_oplen(c[i]);
    } else if (bam_cigar_opchr(c[i]) == 'D') {
      refscan += bam_cigar_oplen(c[i]);
    }
    if (refscan > refpos) {
      if (cigar == 'M' || cigar == 'S' || cigar == 'X' || cigar == '=') {
        readscan = readscan - refscan + refpos;
      }
      break;
    }
  }
  return readscan;
}

bool IsPairOverlap(const SeqLib::BamRecord& one, const SeqLib::BamRecord& two) {
  assert(one.Qname() == two.Qname());
  if (one.Interchromosomal() || one.Position() >= two.PositionEnd() || one.PositionEnd() <= two.Position()) {
    return false;
  } else {
    return true;
  }
}

std::pair<int, int> GetPairOverlapRStartAndRStop(const SeqLib::BamRecord& fwd, const SeqLib::BamRecord& rev) {
  if (!IsPairOverlap(fwd, rev)) {
    return std::make_pair(0,0);
  }
  int left = std::max(fwd.Position(), rev.Position());
  int right = std::min(fwd.PositionEnd(), rev.PositionEnd());
  return std::make_pair(left, right);
}

int OverlapLenInRef(const std::vector<SeqLib::BamRecord>& seg) {
  if (seg.size() == 1) {
    return seg[0].PositionEnd() - seg[0].Position();
  } else {
    auto res = GetPairOverlapRStartAndRStop(seg[0], seg[1]);
    return res.second - res.first;
  }
}

std::pair<std::pair<int,int>,std::pair<int,int>>
GetPairOverlapQStartAndQStop(const SeqLib::BamRecord& fwd, const SeqLib::BamRecord& rev) {
  std::pair<int,int> fwdpos(0,0);
  std::pair<int,int> revpos(0,0);
  if (!IsPairOverlap(fwd, rev)) {
    return std::make_pair(fwdpos, revpos);
  }
  int left = std::max(fwd.Position(), rev.Position());
  int right = std::min(fwd.PositionEnd(), rev.PositionEnd()) - 1;
  fwdpos = std::make_pair(RefPosToQueryPos(fwd, left), RefPosToQueryPos(fwd, right) + 1);
  revpos = std::make_pair(RefPosToQueryPos(rev, left), RefPosToQueryPos(rev, right) + 1);
  return make_pair(fwdpos, revpos);
}

std::pair<int,int>
GetBamOverlapQStartAndQStop(const SeqLib::BamRecord& record, const SeqLib::GenomicRegion& gr) {
  std::pair<int,int> pos(0,0);
  if (gr.GetOverlap(record.AsGenomicRegion()) ==  0) return pos;
  int left = std::max(record.Position(), gr.pos1);
  int right = std::min(record.PositionEnd(), gr.pos2) - 1;
  return std::make_pair(RefPosToQueryPos(record, left), RefPosToQueryPos(record, right) + 1);
}


} //end namespace
