//
// Created by Ruolin Liu on 2/18/20.
//

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>
#include "spoa/spoa.hpp"

#include "Alignment.h"
#include "AlignmentConsensus.h"
#include "BamIO.h"

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

#ifdef BGZF_MAX_BLOCK_SIZE_BAK
#undef BGZF_MAX_BLOCK_SIZE_BAK
#pragma pop_macro("BGZF_MAX_BLOCK_SIZE")
#endif

#ifdef BGZF_BLOCK_SIZE_BAK
#undef BGZF_BLOCK_SIZE_BAK
#pragma pop_macro("BGZF_BLOCK_SIZE")
#endif


using std::string;
struct CssOptions {
  string bam;
  int mapq = 10;
  bool load_supplementary = false;
  bool load_secondary = false;
  bool allow_nonoverlapping_pair = false;
  bool clip3 = false;
  int consensus_mode = 0;
  int pair_min_overlap = 1;
  bool trim_overhang = false;
  string tmpdir = "/tmp";
  //bool output_singleend = false;
  int minbq = 0;
  string outbam = "";
  int thread = 1;
};


static struct option  consensus_long_options[] = {
    {"bam",                      required_argument,      0,        'b'},
    {"load_supplementary",       no_argument,            0,        'l'},
    {"clip3",                    no_argument,            0,        'C'},
    {"mapq",                     required_argument ,     0,        'm'},
    {"baseq",                    required_argument ,     0,        'q'},
    {"outbam",                   required_argument ,     0,        'o'},
    {"pair_min_overlap",         required_argument,      0,        'p'},
    {"trim_overhang",            no_argument,            0,        't'},
    {"allow_nonoverlapping_pair",no_argument,            0,        'i'},
    //{"output_singleend",         no_argument,            0,        's'},
    //hidden parameter
    {"consensus_mode",           required_argument ,     0,        'M'},
    {"dirtmp",                   required_argument ,     0,        'd'},
    {"thread",                   required_argument,      0,        'T'},
    {0,0,0,0}
};

const char* consensus_short_options = "b:m:M:o:lCp:q:d:tT:i";

void consensus_print_help()
{
  std::cerr<< "---------------------------------------------------\n";
  std::cerr<< "Usage: codec consensus [options]\n";
  std::cerr<< "General Options:\n";
  std::cerr<< "-b/--bam,                              Input bam [required]\n";
  std::cerr<< "-o/--outbam,                           Output unmapped bamfile [required].\n";
  std::cerr<< "-m/--mapq,                             Min mapping quality [10].\n";
  std::cerr<< "-q/--baseq,                            paired baseq calibration. If only one of the baseq < cutoff, make the other one baseq = cutoff -1. [0].\n";
  std::cerr<< "-l/--load_supplementary,               Include supplementary alignment [false].\n";
  std::cerr<< "-t/--trim_overhang,                    When perform paired-end consensus, if true then only do consensus of the overlapped region [false].\n";
  std::cerr<< "-C/--clip3,                            trim the 3'end soft clipping [false].\n";
//  std::cerr<< "-s/--output_singleend,                 The R1R2 consensus will be output in a single end format [false].\n";
  std::cerr<< "-p/--pair_min_overlap,                 When using selector, the minimum overlap between the two ends of the pair [1]. "
              "-1 means two segs complete overlap excluding the soft clip part [-1]. 0 for allowing non-overlapping pair\n";
  std::cerr<< "-d/--dirtmp,                           Temporary dir for sorted bam [/tmp]\n";
  std::cerr<< "-T/--thread,                           Number of threads for sort [1]\n";
  std::cerr<< "-i/--allow_nonoverlapping_pair,        Allow output of non-overlaping pairs, usually caused by intermolecular ligation. This will simply print the original reads.  [false]\n";
  //std::cerr<< "-M/--consensus_mode,                   0 for paired end consensus. [0]\n";
}

int consensus_parse_options(int argc, char* argv[], CssOptions& opt) {
  int option_index;
  int next_option = 0;
  do {
    next_option = getopt_long(argc, argv, consensus_short_options, consensus_long_options, &option_index);
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
      case 'l':
        opt.load_supplementary = true;
        break;
      case 't':
        opt.trim_overhang = true;
        break;
      case 'i':
        opt.allow_nonoverlapping_pair = true;
        break;
      case 'C':
        opt.clip3 = true;
        break;
//      case 's':
//        opt.output_singleend = true;
//        break;
      case 'd':
        opt.tmpdir = optarg;
        break;
      case 'M':
        opt.consensus_mode = atoi(optarg);
        break;
      case 'p':
        opt.pair_min_overlap = atoi(optarg);
        break;
      case 'T':
        opt.thread = atoi(optarg);
        break;
      default:consensus_print_help();
        return 1;
    }
  } while (next_option != -1);

  return 0;
}


int codec_consensus(int argc, char ** argv) {
  CssOptions opt;
  int parse_ret =  consensus_parse_options(argc, argv, opt);
  if (parse_ret) return 1;
  if (argc == 1) {
    consensus_print_help();
    return 1;
  }
  if (opt.outbam.empty()) {
    std::cerr << "-o/--outbam is required" << std::endl;
    consensus_print_help();
    return 1;
  }
  if (opt.bam.empty()) {
    std::cerr << "-b/--bam is required" << std::endl;
    consensus_print_help();
    return 1;
  }
  if (opt.allow_nonoverlapping_pair) {
    if (opt.pair_min_overlap != 0) {
      std::cerr << "-p/--pair_min_overlap has to be 0 if -i/--allow_nonoverlapping_pair is true\n";
    }
    return 1;
  }

  char temp[100];
  strcpy(temp, opt.tmpdir.c_str());
  strcat(temp, "/tempsort.XXXXXX");
  int fd = mkstemp(temp);
  if (fd == -1) {
    std::cerr << "unable to create temp file for sorting bam in queryname order\n";
    return 1;
  }
  string samsort = "samtools sort -n " + opt.bam + " -o " + string(temp) + " -@ " + std::to_string(opt.thread);
  std::cout << "sorting bam: " << samsort << std::endl;
  std::system(samsort.c_str());
  std::cout << "sorting done" << std::endl;

  cpputil::InsertSeqFactory isf(temp, opt.mapq, opt.load_supplementary, opt.load_secondary, true, opt.clip3);
  cpputil::UnMappedBamWriter writer(opt.outbam, isf.bamheader());
  int64_t read_counter = 0;
  if (opt.consensus_mode == 0) {
    while (!isf.finished()) {
      std::vector<cpputil::Segments> frag;
      while( (frag= isf.FetchReadNameSorted()).size() > 0) {
        for (auto seg : frag) {
          assert (seg.size() == 2);
          int ol = cpputil::GetNumOverlapBasesPEAlignment(seg);
          if (ol < opt.pair_min_overlap) continue;
          ++ read_counter;
          if (ol > 0) {
            std::vector<std::string> orig_quals; // consensus output
            std::vector<std::string> seqs;
            for (auto&s : seg) {
              seqs.push_back(s.Sequence());
            }
            auto seq = cpputil::MergePair(seg, seqs, opt.trim_overhang, opt.minbq, orig_quals); // trim overhang if true
            if (seg.front().FirstFlag()) {
              writer.WriteRecord(seg.front(), seg.back(), seq.first, seq.second, orig_quals.front(), orig_quals.back());
            }
            else if (seg.back().FirstFlag()) {
              writer.WriteRecord(seg.back(), seg.front(), seq.second, seq.first, orig_quals.back(), orig_quals.front());
            }
          } else if (opt.allow_nonoverlapping_pair) {
            writer.WriteRecord(seg.front(), seg.back());
          }
        }
      }
    }
  } else if (opt.consensus_mode == -1) { // dummy code for denovo consensus
    std::vector<std::string> sequences = {
        "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
        "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
        "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
        "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
        "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
        "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
    };

    auto alignment_engine = spoa::AlignmentEngine::Create(
        spoa::AlignmentType::kOV, 3, -5, -5, -2);  // linear gaps

    spoa::Graph graph{};

    for (const auto& it : sequences) {
      auto alignment = alignment_engine->Align(it, graph);
      graph.AddAlignment(alignment, it);
    }

    auto consensus = graph.GenerateConsensus();

    std::cerr << ">Consensus LN:i:" << consensus.size() << std::endl
              << consensus << std::endl;

    auto msa = graph.GenerateMultipleSequenceAlignment(true);

    for (const auto& it : msa) {
      std::cerr << it << std::endl;
    }
  } else {
    std::cerr << "Consensus mode should be either 0 or 1\n";
    return 1;
  }
  std::cout << "generate " << read_counter << " consensus reads" << std::endl;
  return 0;
}
