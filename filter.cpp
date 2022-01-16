//
// Created by Ruolin Liu on 12/29/21.
//

//
// Created by Ruolin Liu on 2/18/20.
//

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>

#include "Alignment.h"
#include "BamIO.h"

//#ifdef BGZF_MAX_BLOCK_SIZE
//#pragma push_macro("BGZF_MAX_BLOCK_SIZE")
//#undef BGZF_MAX_BLOCK_SIZE
//#define BGZF_MAX_BLOCK_SIZE_BAK
//#endif
//
//#ifdef BGZF_BLOCK_SIZE
//#pragma push_macro("BGZF_BLOCK_SIZE")
//#undef BGZF_BLOCK_SIZE
//#define BGZF_BLOCK_SIZE_BAK
//#endif

#include "InsertSeqFactory.h"

//#ifdef BGZF_MAX_BLOCK_SIZE_BAK
//#undef BGZF_MAX_BLOCK_SIZE_BAK
//#pragma pop_macro("BGZF_MAX_BLOCK_SIZE")
//#endif
//
//#ifdef BGZF_BLOCK_SIZE_BAK
//#undef BGZF_BLOCK_SIZE_BAK
//#pragma pop_macro("BGZF_BLOCK_SIZE")
//#endif


using std::string;
struct FilterOptions {
  string bam;
  string outbam;
  int mapq = 10;
  int minbq = 10;
  int min_family_size = 1;
  bool load_supplementary = false;
  bool load_secondary = false;
};


static struct option  filter_long_options[] = {
    {"bam",                      required_argument,      0,        'b'},
    {"min_family_size",          required_argument,      0,        'f'},
    {"load_supplementary",       no_argument,            0,        's'},
    {"load_secondary",           no_argument,            0,        '2'},
    {"mapq",                     required_argument ,     0,        'm'},
    {"baseq",                    required_argument ,     0,        'q'},
    {"outbam",                   required_argument ,     0,        'o'},
    {0,0,0,0}
};

const char* filter_short_options = "b:s2m:b:o:f:";

void filter_print_help()
{
  std::cerr<< "---------------------------------------------------\n";
  std::cerr<< "Usage: codec consensus [options]\n";
  std::cerr<< "General Options:\n";
  std::cerr<< "-b/--bam,                              Input bam [required]\n";
  std::cerr<< "-o/--outbam,                           Output unmapped bamfile [default stdout].\n";
  std::cerr<< "-m/--mapq,                             Min mapping quality [10].\n";
  std::cerr<< "-f/--min_family_size,                  Minimum number of raw read pairs contributing to the consensus (cD tag in the bam) [1].\n";
  std::cerr<< "-q/--baseq,                            paired baseq calibration. No effect now\n";
  std::cerr<< "-s/--load_supplementary,               Include supplementary alignment [false].\n";
  std::cerr<< "-2/--load_secondary,                   Include secondary alignment [false].\n";
}

int filter_parse_options(int argc, char* argv[], FilterOptions& opt) {
  int option_index;
  int next_option = 0;
  do {
    next_option = getopt_long(argc, argv, filter_short_options, filter_long_options, &option_index);
    switch (next_option) {
      case -1:break;
      case 'b':
        opt.bam = optarg;
        break;
      case 'o':
        opt.outbam = optarg;
        break;
      case 'm':
        opt.mapq = atoi(optarg);
        break;
      case 'q':
        opt.minbq = atoi(optarg);
        break;
      case 'f':
        opt.min_family_size = atoi(optarg);
        break;
      case 's':
        opt.load_supplementary = true;
        break;
      case '2':
        opt.load_secondary = true;
        break;
    }
  } while (next_option != -1);

  return 0;
}


int codec_filter(int argc, char ** argv) {
  FilterOptions opt;
  int parse_ret =  filter_parse_options(argc, argv, opt);
  if (parse_ret) return 1;
  if (argc == 1) {
    filter_print_help();
    return 1;
  }
  if (opt.bam.empty()) {
    std::cerr << "-b/--bam is required" << std::endl;
    filter_print_help();
    return 1;
  }

  cpputil::InsertSeqFactory isf(opt.bam, opt.mapq, opt.load_supplementary, opt.load_secondary, true, false);
  SeqLib::BamWriter bam_writer;
  if (opt.outbam.empty()) {
    bam_writer.Open("-");
  } else {
    bam_writer.Open(opt.outbam);
  }
  bam_writer.SetHeader(isf.bamheader());
  bam_writer.WriteHeader();

  int64_t pass_counter = 0;
  int64_t failed_counter = 0;
  std::string dup;
  int fsize;
  while (!isf.finished()) {
    std::vector<cpputil::Segments> frag;
    while( (frag= isf.FetchReadNameSorted()).size() > 0) {
      for (auto seg : frag) {
        if (seg.size() != 2) {
          std::cerr << seg[0] << std::endl;
        }
        if (!seg.front().FirstFlag()) {
          std::cerr << "first flag first!!\n";
        }
        if (seg.front().GetZTag("du", dup)) {
          ++failed_counter;
        } else {
          if (seg.front().GetIntTag("cD", fsize)) {
            if (fsize >= opt.min_family_size) {
              ++pass_counter;
              bam_writer.WriteRecord(seg.front());
              bam_writer.WriteRecord(seg.back());
            } else {
              ++failed_counter;
            }
          } else {
            ++pass_counter;
            bam_writer.WriteRecord(seg.front());
            bam_writer.WriteRecord(seg.back());
          }
        }
      }
    }
  }
  std::cerr << "kept " << pass_counter << " consensus pairs" << std::endl;
  std::cerr << "discard " << failed_counter << " consensus pairs" << std::endl;
  return 0;
}