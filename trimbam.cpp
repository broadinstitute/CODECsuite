//
// Created by Ruolin Liu on 3/13/20.
//

#include <iostream>
#include <string>
#include <getopt.h>

#include "SeqLib/BamWriter.h"
#include "SeqLib/BamReader.h"
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
#include "BamRecordExt.h"
#include "InsertSeqFactory.h"

using std::string;
struct Options {
  string bam;
  int32_t trim5 = 0;
  int32_t trim3 = 0;
  string trimbam;
  bool mate_cigar_missing = false;
  int mapq = 0;
  bool load_supplementary = true;
  bool clip3 = false;
  const bool load_unpair = true;
  int consensus_mode = 0;
};


static struct option  long_options[] = {
    {"bam",                      required_argument,      0,        'b'},
    {"trim_bam",                 required_argument,      0,        'T'},
    {"trim5",                    required_argument,      0,        '5'},
    {"trim3",                    required_argument ,     0,        '3'},
    {"mate_cigar_missing",       no_argument ,           0,        'm'},
    {"clip3",                    no_argument ,           0,        'c'},
    {0,0,0,0}
};

const char*short_options = "b:5:3:T:mc";

void print_help()
{
  std::cerr<< "---------------------------------------------------\n";
  std::cerr<< "Usage: consensus [options] -b a.bam\n";
  std::cerr<< "General Options:\n";
  std::cerr<< "-b/--bam,                              input bam\n";
  std::cerr<< "-5/--trim5,                            Trim number of bases based on this option from the 5'end of the fragment (not read), by change them to N. [0]\n";
  std::cerr<< "-3/--trim3,                            Trim number of bases based on this option from the 3'end of the fragment (not read), by change them to N. [0]\n";
  std::cerr<< "-T/--trim_bam,                         Output file for trimmed bam [null]\n";
  std::cerr<< "-m/--mate_cigar_missing,               No mate cigar in the input bam. The algorithm will close pair first and this is slower.[false]\n";
  std::cerr<< "-c/clip3,                              Trim the 3'end of the read if it is soft-clipped. It is likely to be adapter.[false]\n";
}

int parse_options(int argc, char* argv[], Options& opt) {
  int option_index;
  int next_option = 0;
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
    switch (next_option) {
      case -1:break;
      case 'b':
        opt.bam = optarg;
        break;
      case '5':
        opt.trim5 = atoi(optarg);
        break;
      case '3':
        opt.trim3 = atoi(optarg);
        break;
      case 'T':
        opt.trimbam = optarg;
        break;
      case 'm':
        opt.mate_cigar_missing = true;
        break;
      case 'c':
        opt.clip3 = true;
        break;
      default:print_help();
        return 1;
    }
  } while (next_option != -1);

  return 0;
}

int main(int argc, char ** argv) {
  Options opt;
  int parse_ret = parse_options(argc, argv, opt);
  if (parse_ret) return 1;
  if (argc == 1) {
    print_help();
    exit(0);
  }

  SeqLib::BamReader bam_reader_;
  SeqLib::BamWriter trim_bam_writer;
  bam_reader_.Open(opt.bam);
  if (opt.trimbam.empty()) {
    std::cerr << "trimbam output is missing. Use -T\n";
    exit(0);
  }
  trim_bam_writer.SetHeader(bam_reader_.Header());
  trim_bam_writer.Open(opt.trimbam);
  trim_bam_writer.WriteHeader();

//  if (opt.mate_cigar_missing) {
//    cpputil::InsertSeqFactory isf(opt.bam, opt.mapq, opt.load_supplementary, true,);
//    while (!isf.finished()) {
//      std::vector<std::vector<cpputil::Segments>> chunk;
//      while ((chunk = isf.ReadByChrom(cpputil::SegmentNotEmpty, 0, opt.clip3, opt.consensus_mode, opt.load_unpair)).size() > 0) {
//
//        for (auto frag : chunk) {
//          for (auto seg : frag) {
//            if (seg.size() == 2) {
//              cpputil::TrimPairFromFragEnd(seg.front(), seg.back(), opt.trim5, opt.trim3);
//              trim_bam_writer.WriteRecord(seg.front());
//              trim_bam_writer.WriteRecord(seg.back());
//            } else if (seg.size() == 1){
//              if (seg.front().ReverseFlag()) {
//                cpputil::TrimBamFromFragEnd(seg.front(), 0, 0, 0, opt.trim3);
//                trim_bam_writer.WriteRecord(seg.front());
//              } else {
//                cpputil::TrimBamFromFragEnd(seg.front(), 0, 0, opt.trim5, 0);
//                trim_bam_writer.WriteRecord(seg.front());
//              }
//            }
//          }
//        }
//
//      }
//    }
//
//  } else {
    SeqLib::BamRecord b;
    while(true) {
      auto ret = bam_reader_.GetNextRecord(b);
      if (!ret) break;
      if (opt.clip3) {
        cpputil::SoftClip3end(b);
      }
      int32_t mp, mep;
      std::tie(mp, mep) = cpputil::MatePositionAndPositionEndWithSoftClip(b);
      cpputil::TrimBamFromFragEnd(b, mp, mep, opt.trim5, opt.trim3);
      trim_bam_writer.WriteRecord(b);
    }
//  }
  return 0;
}
