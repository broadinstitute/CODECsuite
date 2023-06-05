//
// Created by Ruolin Liu on 3/13/20.
//

#include <iostream>
#include <getopt.h>
#include <cassert>

#include "Index.h"
#include "Files.h"

using std::string;
struct DemuxOptions {
  string library_file;
  string fastq1;
  string fastq2;
  string outprefix = "./test";
  string reference;
  int index_begin = 3;
  int index_len = 18;
  int max_ed = 2;
  int min_readlen = 30;
  bool include_non_pf = false;
  bool verbose = false;
  bool out_unmatch = false;
  bool out_hopped = false;
  bool count_PF = false;
};


static struct option  demux_long_options[] = {
    {"library_param",            required_argument,      0,        'p'},
    {"q1",                       required_argument,      0,        '1'},
    {"q2",                       required_argument,      0,        '2'},
    {"outprefix",                required_argument ,     0,        'o'},
    {"ref",                      required_argument ,     0,        'r'},
    {"index_begin",              required_argument,      0,        'b'},
    {"index_len",                required_argument,      0,        'l'},
    {"min_read_len",             required_argument,      0,        'm'},
    {"max_ed",                   required_argument,      0,        'e'},
    {"verbose",                  no_argument ,           0,        'v'},
    {"out_unmatch",                  no_argument ,           0,        'u'},
    {"out_hopped",                  no_argument ,           0,        'h'},
    {"include_non_pf",           no_argument ,           0,        'i'},
    {"count_pf",                 no_argument ,           0,        'c'},
    {0,0,0,0}
};

const char* demux_short_options = "p:1:2:o:r:vib:l:e:cuh:m:";

void codec_demux_usage()
{
  std::cerr<< "---------------------------------------------------\n";
  std::cerr<< "Usage: codec demux [options]\n";
  std::cerr<< "General Options:\n";
  std::cerr<< "-p/--library_param,                    Sample, barcode mapping in CSV format. Header must be \"SampleName,IndexBarcode1,IndexBarcode2\"\n";
  std::cerr<< "-1/--q1,                               Input read1\n";
  std::cerr<< "-2/--q2,                               Input read2\n";
  std::cerr<< "-b/--index_begin,                      The read position where the index begins (Default: 3) \n";
  std::cerr<< "-l/--index_len,                        Index length (Default: 18)\n";
  std::cerr<< "-m/--min_read_len,                     Minimum read length (Default: 30)\n";
  std::cerr<< "-e/--max_ed,                           Maximum edit distance allowed as a match (Default: 2)\n";
  std::cerr<< "-o/--outprefix,                        Output path, e.g., /tmp/test\n";
  std::cerr<< "-r/--ref,                              Reference genome fasta file, for judging index hopping\n";
  std::cerr<< "-i/--include_non_pf,                   Include non-pass filter reads\n";
  std::cerr<< "-v/--verbose,                          Print verbose information\n";
  std::cerr<< "-c/--count_pf,                         Just count number of pass filter pairs. Do not do anything else\n";
  std::cerr<< "-u/--out_unmatch,                      Output reads having no matching barcodes\n";
  std::cerr<< "-u/--out_hopped,                       Output reads hopped\n";
}

int demux_parse_options(int argc, char* argv[], DemuxOptions& opt) {
  int option_index;
  int next_option = 0;
  do {
    next_option = getopt_long(argc, argv, demux_short_options, demux_long_options, &option_index);
    switch (next_option) {
      case -1:break;
      case 'p':
        opt.library_file = optarg;
        break;
      case '1':
        opt.fastq1 = optarg;
        break;
      case '2':
        opt.fastq2 = optarg;
        break;
      case 'b':
        opt.index_begin = atoi(optarg);
        break;
      case 'l':
        opt.index_len = atoi(optarg);
        break;
      case 'm':
        opt.min_readlen = atoi(optarg);
        break;
      case 'e':
        opt.max_ed = atoi(optarg);
        break;
      case 'o':
        opt.outprefix = optarg;
        break;
      case 'r':
        opt.reference = optarg;
        break;
      case 'i':
        opt.include_non_pf = true;
        break;
      case 'v':
        opt.verbose = true;
        break;
      case 'c':
        opt.count_PF = true;
        break;
      case 'h':
        opt.out_hopped = true;
        break;
      case 'u':
        opt.out_unmatch = true;
        break;
      default:codec_demux_usage();
        return 1;
    }
  } while (next_option != -1);

  return 0;
}

int codec_demux(int argc, char ** argv) {
  DemuxOptions opt;
  int parse_ret = demux_parse_options(argc, argv, opt);
  if (parse_ret) return 1;
  if (argc == 1) {
    codec_demux_usage();
    return 1;
  }
  //AffineGap ag("CACTGATCGTCAGCTGAC", "TGAATCTGAGGCACTGTA");
//  AffineGap ag("AGT", "TGAGTT");
//  ag.PrintAllPaths();
//  exit(0);

  CDS::IndexBarcode ibmatcher(opt.library_file, opt.outprefix, opt.max_ed, opt.out_unmatch, opt.out_hopped, opt.verbose);
  if (!opt.reference.empty()) {
    ibmatcher.LoadBwa(opt.reference);
  }
  if (not cpputil::FileExist(opt.library_file)) {
    std::cerr << opt.library_file << " does not exist\n";
    return 1;
  }
  cpputil::FastxReader R1_reader(opt.fastq1);
  cpputil::FastxReader R2_reader(opt.fastq2);
  cpputil::FastxRecord read1;
  cpputil::FastxRecord read2;
  uint64_t total_pf_reads = 0;
  uint64_t total_reads = 0;
  while (R1_reader.yield(read1)) {
    R2_reader.yield(read2);
    ++total_reads;
    if (not opt.include_non_pf and (read1.is_filtered() or read2.is_filtered())) continue;
    if (read1.seq.length() < opt.min_readlen or read2.seq.length() < opt.min_readlen) continue;
    ++total_pf_reads;
    //if (opt.count_PF) continue;
    assert(read1.name() == read2.name());
    ibmatcher.DecodePair(read1, read2, opt.index_begin, opt.index_len);
  }
//  if (opt.count_PF) {
uint64_t total_matched = ibmatcher.total_matched();
for (unsigned i = 0; i < ibmatcher.samples().size(); ++ i) {
  string s = ibmatcher.samples()[i];
  uint64_t n = ibmatcher.nmatched()[i];
  std::cout << "#sample, matched, matched%: " << s << ", " << n << ", " << (double) n / total_matched << std::endl;
}
std::cout << "#total, #PF, #matched, matched%: " << total_reads << ", " << total_pf_reads << ", " << total_matched << ", " << (double) total_matched / total_pf_reads << std::endl;
//  }
  return 0;
}
