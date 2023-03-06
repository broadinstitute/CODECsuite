//
// Created by Ruolin Liu on 3/3/20.
//
//
// Created by Ruolin Liu on 2/18/20.
//

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include <chrono>
#include <ctime>

//#include "SeqLib/BWAWrapper.h"
#include "SeqLib/BamWriter.h"
#ifdef BGZF_MAX_BLOCK_SIZE
#pragma push_macro("BGZF_MAX_BLOCK_SIZE")
#undef BGZF_MAX_BLOCK_SIZE
#define BGZF_MAX_BLOCK_SIZE_BAK
#endif

#ifdef BGZF_BLOCK_SIZE
#pragma push_macro("BGZF_BLOCK_SIZE")
#undef BGZF_BLOCK_SIZE
#define BGZF_BLOCK_SIZE_BAK
#endif

#include "InsertSeqFactory.h"
#include "ReadVCF.h"
#include "Variant.h"
#include "BamRecordExt.h"
#include "Alignment.h"
#include "StringUtils.h"
#include "MAF.h"
#include "TargetLayout.h"
#include "Stats.h"

using std::string;
using std::vector;
int MYINT_MAX = std::numeric_limits<int>::max();

struct Options {
  string tumor_bam;
  string normal_bam;
  string outprefix;
  int mapq = 10;
  int bqual_min = 0;
  bool load_supplementary = false;
  bool load_secondary = false;
  bool load_unpair = false;
  bool load_duplicate = false;
  bool allow_normal_nocov = false;
  int max_mismatch_filter = MYINT_MAX;
  int verbose = 0;
  int anchor = 8;
  string vcf;
  //string sample = "";
  string reference;
  int fragend_dist_filter = 0;
  bool overlap_only = false;
  string bed_file;
  string uid_tag_name;
};


static struct option  long_options[] = {
    {"tumor_bam",                required_argument,      0,        't'},
    {"normal_bam",               required_argument,      0,        'n'},
    {"bed",                      required_argument,      0,        'L'},
    {"outprefix",                required_argument,      0,        'o'},
    {"unique_molecular_id",      required_argument,      0,        'U'},
    {"load_unpair",              no_argument,            0,        'u'},
    {"load_supplementary",       no_argument,            0,        'S'},
    {"load_secondary",           no_argument,            0,        '2'},
    {"load_duplicate",           no_argument,            0,        'D'},
    {"allow_normal_nocov",       no_argument,            0,        'e'},
    {"vcf",                      required_argument ,     0,        'V'},
    {"maf",                      required_argument ,     0,        'M'},
    //{"sample",                   required_argument,      0,        's'},
    {"reference",                required_argument,      0,        'r'},
    {"bqual_min",                required_argument,      0,        'q'},
    {"anchor",                   required_argument,      0,        'a'},
    {"overlap_only",             no_argument,            0,        'O'},
    {"fragend_dist_filter",      required_argument,      0,        'd'},
    {"max_mismatch_filter",      required_argument,      0,        'x'},
    {"verbose",                  required_argument,      0,        'v'},
    {0,0,0,0}
};

//const char* msi_short_options = "b:m:v:S2us:r:q:Od:n:x:V:L:t:";
const char* msi_short_options = "b:m:v:S2Dur:q:Od:n:x:L:t:U:a:o:V:e";

void msi_print_help()
{
  std::cerr<< "---------------------------------------------------\n";
  std::cerr<< "Usage: msi -t tumor.bam -n normal -L region.bed [options]\n";
  std::cerr<< "General Options:\n";
  std::cerr<< "-v/--verbose,                          [default 0]\n";
  std::cerr<< "-t/--tumor_bam,                        input bam\n";
  std::cerr<< "-n/--normal_bam,                       input bam\n";
  std::cerr<< "-L/--bed,                              targeted region\n";
  std::cerr<< "-o/--outprefix,                        output prefix. outprefix.all, outprefix.msi\n";
  std::cerr<< "-m/--mapq,                             min mapping quality [10].\n";
  std::cerr<< "-S/--load_supplementary,               include supplementary alignment [false].\n";
  std::cerr<< "-2/--load_secondary,                   include secondary alignment [false].\n";
  std::cerr<< "-D/--load_duplicate,                   include duplicated alignment [false].\n";
  std::cerr<< "-u/--load_unpair,                      include unpaired alignment [false].\n";
  std::cerr<< "-e/--allow_normal_nocov,               Allow sites where normal sample has 0 coverage [false].\n";
  std::cerr<< "-U/--unique_molecular_id,              UID tag name [""].\n";
  std::cerr<< "-V/--vcf,                              .\n";
  //std::cerr<< "-s/--sample,                           sample from the VCF file [null].\n";
  std::cerr<< "-r/--reference,                        reference sequence in fasta format [null].\n";
  std::cerr<< "\nFiltering Options:\n";
  std::cerr<< "-q/--bqual_min,                        Not counting a base if its quality is less than this value [0].\n";
  std::cerr<< "-x/--max_nonNedit_filter,              Skip a read if the number of non-N bases edits is larger than this value [INT_MAX].\n";
  std::cerr<< "-d/--fragend_dist_filter,              Consider a variant if its distance to the fragment end is at least this value [0].\n";
  std::cerr<< "-a/--anchor,                           number of anchor bases required for spanning repeats [8].\n";
  std::cerr<< "-O/--overlap_only,                     Count only overlapped region of a read pair. Default [false].\n";
}

int msi_parse_options(int argc, char* argv[], Options& opt) {
  int option_index;
  int next_option = 0;
  do {
    next_option = getopt_long(argc, argv, msi_short_options, long_options, &option_index);
    switch (next_option) {
      case -1:break;
      case 'v':
        opt.verbose = atoi(optarg);
        break;
      case 't':
        opt.tumor_bam = optarg;
        break;
      case 'n':
        opt.normal_bam = optarg;
        break;
      case 'L':
        opt.bed_file = optarg;
        break;
      case 'o':
        opt.outprefix = optarg;
        break;
      case 'U':
        opt.uid_tag_name = optarg;
        break;
      case 'm':
        opt.mapq = atoi(optarg);
        break;
      case 'q':
        opt.bqual_min = atoi(optarg);
        break;
      case 'x':
        opt.max_mismatch_filter = atoi(optarg);
        break;
      case 'd':
        opt.fragend_dist_filter = atoi(optarg);
        break;
      case 'a':
        opt.anchor = atoi(optarg);
        break;
      case 'S':
        opt.load_supplementary = true;
        break;
      case '2':
        opt.load_secondary = true;
        break;
      case 'u':
        opt.load_unpair = true;
        break;
      case 'D':
        opt.load_duplicate = true;
        break;
      case 'e':
        opt.allow_normal_nocov = true;
        break;
      case 'O':
        opt.overlap_only = true;
        break;
//      case 's':
//        opt.sample = optarg;
//        break;
      case 'V':
        opt.vcf = optarg;
        break;
      case 'r':
        opt.reference = optarg;
        break;
      default:msi_print_help();
        return 1;
    }
  } while (next_option != -1);

  return 0;
}

struct Repeat {
  int32_t refstart;
  std::string pat;
  int32_t nunit;
  char previous_base;
  double pre_entropy;
  double post_entropy;
  Repeat(): refstart(0), pat(""), nunit(0), previous_base(0), pre_entropy(-1.0), post_entropy(-1.0){}
  Repeat(int32_t p, std::string u, int32_t n, char pb): refstart(p), pat(u), nunit(n), previous_base(pb), pre_entropy(-1.0), post_entropy(-1.0){}
  int32_t refend() const {
    return refstart + (int) pat.size() * nunit;
  }
  friend bool operator<(const Repeat&l, const Repeat& r) {
    return std::tie(l.refstart, l.pat, l.nunit) < std::tie(r.refstart, r.pat, r.nunit);
  }
  float VariantAf(const cpputil::BCFReader& bcf_reader, const std::string contig, int var) {
    // var: number of bases at INS (>0) or DEL (<0)
    std::string ref = std::string(1, previous_base);
    std::string alt = std::string(1, previous_base);
    if (var < 0) { // current del only;
      for (int i = var; i < 0; ++i) {
        ref += pat;
      }
    } else {
      for (int i = 0; i < var; ++i) {
        alt += pat;
      }
    }
    return bcf_reader.var_population_freq(contig, refstart - 1, ref, alt);
  }
};

std::ostream& operator<<(std::ostream& os, const Repeat& rp) {
  os << rp.refstart <<"\t" << rp.pat <<"\t" << rp.nunit;
  return os;
}

bool IsHomopolymer(const std::string& seq) {
  assert(seq.size() > 1);
  char c = seq[0];
  for(auto i = 1; i < seq.size(); ++i) {
    if (seq[i] != c) return false;
  }
  return true;
}

//search up to simple quadTR
std::vector<Repeat> ScanRef(const SeqLib::RefGenome& ref, const SeqLib::GenomicRegion& gr, const SeqLib::BamHeader& bh,
            const int hp_min = 10, const int tr_min = 5, const int period_max = 4) {
  auto seq = ref.QueryRegion(gr.ChrName(bh), gr.pos1 - 1, gr.pos2); // getting an extra previous base
  //make sure no 'N'
  std::vector<Repeat> res;
  for (unsigned ii = 1; ii < seq.size(); ++ii) {
    if (seq[ii] == 'N') return res;
  }
  for (int period = 1; period  <= period_max; period++) {
    for (int shift = 0; shift < period; shift++) {
      string prev = seq.substr(shift, period);
      int strike = 1;
      //bool found = false;
      for (unsigned ii = period + shift; ii < seq.size();) {
        string curr = seq.substr(ii, period);
        if (period > 1 and IsHomopolymer(prev)){
          prev = curr;
          ii += period;
          continue;
        }
        if (curr == prev) {
          strike++;
        }
        else { // reach the end of repeat
          if ((period == 1 && strike >= hp_min) or (period > 1 && strike >= tr_min)) {
            int32_t s = ii - period * strike + gr.pos1 - 1;
            char baseprev = seq[ii - period * strike - 1];
            bool overlapped = false;
            for (unsigned tt = 0; tt < res.size(); ++tt) {
              if (s < res[tt].refend()) {
                overlapped = true;
                if (strike > res[tt].nunit or (strike == res[tt].nunit && s < res[tt].refstart)) { // select a better one
                  res[tt] = Repeat(s, prev, strike, baseprev);
                }
                break;
              }
            }
            if (not overlapped) {
              res.emplace_back(s, prev, strike, baseprev);
              //found = true;
            }
          }
          strike = 1;
        }
        prev = curr;
        ii += period;
      }
      //if (found) break;
    }
  }
  return res;
}

bool SegOverlapRepeat(const cpputil::Segments& segs, const Repeat& re, const bool allow_indel_in_anchor, const bool allow_N_in_anchor, const int anchor) {
  /*
   * Return  num. of DEL(<0) or INS(>0)
   */
  for (const auto& br : segs) {
    if (br.Position() + anchor > re.refstart or br.PositionEnd() < re.refend() + anchor) {
      return 0;
    }
    // closed intervals
    int front_anchor_start = cpputil::RefPosToQueryPos(br, re.refstart - anchor);
    int front_anchor_end = cpputil::RefPosToQueryPos(br, re.refstart - 1);
    int back_anchor_start = cpputil::RefPosToQueryPos(br, re.refend());
    int back_anchor_end = cpputil::RefPosToQueryPos(br, re.refend() + anchor - 1);
    if (not allow_indel_in_anchor and (front_anchor_end - front_anchor_start + 1 != anchor || back_anchor_end - back_anchor_start + 1 != anchor)) {
      return 0;
    }
    std::string prefix =  br.Sequence().substr(front_anchor_start, front_anchor_end - front_anchor_start + 1);
    std::string suffix =  br.Sequence().substr(back_anchor_start, back_anchor_end - back_anchor_start + 1);
    if (not allow_N_in_anchor and (std::count(prefix.begin(), prefix.end(), 'N') or std::count(suffix.begin(), suffix.end(), 'N')) ) {
      return 0;
    }

  }
  return 1;
}

int ScanRead(const SeqLib::BamRecord& br, const Repeat& re, int& num_n, bool allow_disrupt, bool allow_endwith_different_repeat, const int anchor) {
  /*
   * if the repeat is not continuous return -1;
   */
  assert (br.Position() + anchor <= re.refstart and br.PositionEnd() >= re.refend() + anchor);
  int readstart = cpputil::RefPosToQueryPos(br, re.refstart);
  int readend = cpputil::RefPosToQueryPos(br, re.refend()-1);
  std::string rs = br.Sequence().substr(readstart, readend - readstart + 1);
  num_n = std::count(rs.begin(), rs.end(), 'N');
  int occ = 0;
  std::string::size_type pos =0;
  std::string obs = "";
  while((pos = rs.find(re.pat, pos)) != std::string::npos) {
    ++occ;
    obs += re.pat;
    pos += re.pat.length();
  }

  // not a single unit observed, likely a large deletion overlapped
  if (occ == 0 or rs.length() < re.pat.length()) {
    return -1;
  }
  //if pattern is disruputed
  if (not allow_disrupt and rs.find(obs) == std::string::npos) {
    return -1;
  }
  //if last unit is not the pattern
  if (not allow_endwith_different_repeat and rs.substr(rs.length() - re.pat.length(), re.pat.length()) != re.pat) {
    return -1;
  }
  return occ;
//  int indel_len = 0;
//  if (occ + num_n < re.nunit) {
//    indel_len = occ + num_n - re.nunit;
//  } else if (occ > re.nunit) {
//    indel_len = occ - re.nunit;
//  }
  //return indel_len;
}

bool HomozygousRef(const Repeat& rp, const std::vector<int>& obs) {
  /*
   * Return true if homozygous reference
   */
  bool homo = true;
  for (auto& o : obs) {
    if (o != rp.nunit) {
      return false;
    }
  }
  return true;
}

int NonRefGt(const Repeat& rp, const std::vector<int>& obs, const int num_diff_small = 2, const int num_diff_large = 3) {
  /*
   * Return lenght. <0 DEL, >0 Insert
   */
  std::vector<int> nonref;
  const int large_del = 20;
  int num_diff = rp.nunit >= large_del ? num_diff_large : num_diff_small;
  for (auto& o : obs) {
    if (o + num_diff <= rp.nunit or o >= rp.nunit + num_diff)
      nonref.push_back(o - rp.nunit);
  }
  if (nonref.empty()) {
    return 0;
  } else {
    int mostlik = cpputil::GetMode(nonref);
    return mostlik;
  }
}

std::string HomoLen2Str(const std::vector<int>& obs) {
  if(obs.empty()) {
    return "NA";
  }
  std::string str = "";
  for (const auto& d : obs) {
    str += std::to_string(d);
    str +=",";
  }
  str = str.substr(0, str.length()-1);
  return str;
}

int main(int argc, char ** argv) {
  Options opt;
  int parse_ret =  msi_parse_options(argc, argv, opt);
  if (parse_ret) return 1;
  if (argc == 1) {
    msi_print_help();
    return 1;
  }
  cpputil::BCFReader bcf_reader;
  if (!opt.vcf.empty()) {
    bcf_reader.Open(opt.vcf.c_str());
  }
  if (opt.outprefix.empty()) {
    std::cerr << "\n Error: output prefix must be specified" << std::endl;
    msi_print_help();
    return 1;
  }
  SeqLib::RefGenome ref;
  ref.LoadIndex(opt.reference);
  cpputil::InsertSeqFactory tumor_factory(opt.tumor_bam, opt.mapq, opt.load_supplementary, opt.load_secondary, opt.load_duplicate, true, false);
  cpputil::InsertSeqFactory normal_factory;
  std::ofstream all(opt.outprefix + ".all");
  std::ofstream msi(opt.outprefix + ".msi");
  if (!opt.normal_bam.empty()) {
    normal_factory = cpputil::InsertSeqFactory(opt.normal_bam, opt.mapq, opt.load_supplementary, opt.load_secondary, false, false, false);
  }
  cpputil::TargetLayout tl(tumor_factory.bamheader(), opt.bed_file);
  vector<vector<cpputil::Segments>> tumor_chunk;
  vector<vector<cpputil::Segments>> normal_chunk;
  for (int i = 0; i < tl.NumRegion(); ++i) {
    const auto& gr = tl[i];
    if (i % 1000 == 0) {
      auto timenow = std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
      std::cerr << i + 1 << " region processed. Last position: " << gr << std::ctime(&timenow);
    }
    auto repeats = ScanRef(ref, gr, tumor_factory.bamheader());
    vector<Repeat> filtered_repeat;
    int window = 10;
    for (auto& rep : repeats) {
      auto prefix = ref.QueryRegion(gr.ChrName(tumor_factory.bamheader()), rep.refstart - window, rep.refstart-1);
      auto suffix = ref.QueryRegion(gr.ChrName(tumor_factory.bamheader()), rep.refend(), rep.refend() + window - 1);
      double pe = cpputil::entropy(prefix);
      double se = cpputil::entropy(suffix);
      //std::cout << gr.ChrName(tumor_factory.bamheader()) << "\t" << rep <<"\t" << prefix << "\t" << pe << "\t" << suffix << "\t" << se << "\t" << std::endl;
      if (rep.pat.length() == 1) { // if homopolyer, check 10bp before and after if the pattern over represented
        rep.pre_entropy = (double) std::count(prefix.begin(), prefix.end(), rep.pat[0]) / prefix.length();
        rep.post_entropy = (double) std::count(suffix.begin(), suffix.end(), rep.pat[0]) / suffix.length();
      } else {
        rep.pre_entropy = 0;
        rep.post_entropy = 0;
      }
      if (pe > 1.0 and se > 1.0) {
        filtered_repeat.push_back(rep);
      }
    }
    repeats = filtered_repeat;
    const int max_indel_len = 100;
    std::vector<int> rp_hist(max_indel_len + 1);
    std::map<Repeat, std::vector<int>> tumor_obs;
    if (tumor_factory.ReadByRegion(cpputil::ArePEAlignmentOverlapAtLeastK, gr, tumor_chunk, 0, opt.uid_tag_name, opt.load_unpair)) {
      for (auto rp : repeats) {
        for (auto frag : tumor_chunk) {
          std::fill(rp_hist.begin(), rp_hist.end(), 0);
          if (opt.verbose) {
            std::cerr << "frag size: " << frag.size() << std::endl;
          }
          for (auto &seg: frag) {

            bool keep = true;
            for (const auto &s: seg) {
              if ( cpputil::GetNM(s) - cpputil::CountNBasesInAlignment(s) - cpputil::GetTotalIndelLen(s) > opt.max_mismatch_filter) {
                if (opt.verbose) {
                  std::cerr << "Discard read by mixmatch filter " << s.Qname() <<"\n";
                }
                keep = false;
                break;
              }
            }
            if (not keep) continue;

            if (!SegOverlapRepeat(seg, rp, false, false, opt.anchor)) continue;
            int num_n1, num_n2;
            auto r1 = ScanRead(seg[0], rp, num_n1, false, false, opt.anchor);
            auto r2 = ScanRead(seg[1], rp, num_n2, false, false, opt.anchor);
            if (r1 != r2 or num_n1 or num_n2 or // if disagree between pair or has N
                r1 == -1 or r2 == -1) { // or non-perfect repeat
              continue;
            }
            //debug info
            if (opt.verbose) {
              std::cerr << "T: " << gr.ChrName(tumor_factory.bamheader()) << "\t" << rp << "\t" << seg[0].Qname() << "\t" << r1 << std::endl;
            }
            ++rp_hist[r1];
          }
          int mode_pos = std::max_element(rp_hist.begin(), rp_hist.end()) - rp_hist.begin();
          if (rp_hist[mode_pos] * 2 > (int) frag.size()) {
            tumor_obs[rp].push_back(mode_pos);
          }
        }
      } // end for
    } //end if

    std::map<Repeat, std::vector<int>> normal_obs;
    if (!opt.normal_bam.empty()) {
      if (normal_factory.ReadByRegion(cpputil::ArePEAlignmentOverlapAtLeastK, gr, normal_chunk, 0, "", opt.load_unpair)) {
        for (auto rp : repeats) {
          for (auto frag : normal_chunk) {
            for (auto seg: frag) {


              if (!SegOverlapRepeat(seg, rp, true, true, opt.anchor)) continue;
              int num_n1, num_n2;
              auto r1 = ScanRead(seg[0], rp, num_n1, true, true, opt.anchor);
              auto r2 = ScanRead(seg[1], rp, num_n2, true, true, opt.anchor);
              if (r1 != r2 or r1 == -1 or r2 == -1) { // if disagree or non perfect repeat
                continue;
              }
              if (opt.verbose) {
                std::cerr << "N: " << gr.ChrName(tumor_factory.bamheader()) << "\t" << rp << "\t" << seg[0].Qname() << "\t" << std::endl;
              }
              normal_obs[rp].push_back(r1);
            }
          }
        } // end for
      } //end if
    }

    //process repeats in the region
    for (auto rp : repeats) {
      //if (!opt.normal_bam.empty() && (normal_obs[rp].empty() || !HomozygousRef(rp, normal_obs[rp])) ) {
      if (!opt.normal_bam.empty() &&
          (!HomozygousRef(rp, normal_obs[rp]) || (!opt.allow_normal_nocov && normal_obs[rp].empty()) ) ) {
        continue;
      }
      auto normal_str = opt.normal_bam.empty()? "NA" : HomoLen2Str(normal_obs[rp]);
      if (!tumor_obs[rp].empty()) {
        auto tumor_str = HomoLen2Str(tumor_obs[rp]);
        string af = "NA";
        all << gr.ChrName(tumor_factory.bamheader()) << "\t" << rp << "\t" << tumor_str << "\t" << normal_str << "\t"
            << rp.pre_entropy << "\t" << rp.post_entropy <<"\t" << af
            << std::endl;
        int n = NonRefGt(rp, tumor_obs[rp]);
        if (n < 0) { // deletion only
          if (!opt.vcf.empty()) {
            float f = rp.VariantAf(bcf_reader, gr.ChrName(tumor_factory.bamheader()), n);
            af = std::to_string(f);
          }
          msi << gr.ChrName(tumor_factory.bamheader()) << "\t" << rp << "\t" << tumor_str << "\t" << normal_str << "\t"
                    << rp.pre_entropy << "\t" << rp.post_entropy << "\t" << af
                    << std::endl;
        }
      }
    }
  } //end for
  std::cerr << "All region processed \n";
  return 0;
}
