//
// Created by Ruolin Liu on 2/18/20.
//

#include <iostream>
#include <fstream>
#include <string>
#include <getopt.h>

#include "Alignment.h"
#include "AlignmentConsensus.h"
#include "BamIO.h"
#include "MutCounter.h"

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
  string reference;
  int mapq = 10;
  bool output_nononverlapping_pair = false;
  bool clip3 = false;
  int consensus_mode = 0;
//  int pair_min_overlap = 1;
  bool trim_overhang = false;
  string tmpdir = "/tmp";
  //bool output_singleend = false;
  int bqual_min = 0;
  string outbam = "";
  int thread = 1;

  int min_fraglen = 30;
  int max_fraglen = std::numeric_limits<int>::max();
  float min_passQ_frac = 0;
  bool filter_5endclip = false;
  int max_N_filter = std::numeric_limits<int>::max();
  int max_snv_filter = std::numeric_limits<int>::max();
  int clustered_mut_cutoff = std::numeric_limits<int>::max();
  float max_frac_prim_AS = std::numeric_limits<float>::max();
  float max_pair_mismatch_frac = 1.0;

  //hidden options
  int count_read = 0;
  bool standard_ngs_filter = false;
  bool load_unpair = true;
  bool load_supplementary = false;
  bool load_secondary = false;
};


static struct option  filter_long_options[] = {
    {"bam",                      required_argument,      0,        'b'},
    {"reference",                required_argument,      0,        'r'},
//    {"load_supplementary",       no_argument,            0,        'l'},
//    {"load_unpair",              no_argument,            0,        'u'},
//    {"load_secondary",           no_argument,            0,        '2'},
//    {"load_duplicate",           no_argument,            0,        'D'},
    {"clip3",                    no_argument,            0,        'C'},
    {"mapq",                     required_argument ,     0,        'm'},
    {"baseq",                    required_argument ,     0,        'q'},
    {"outbam",                   required_argument ,     0,        'o'},
//    {"pair_min_overlap",         required_argument,      0,        'p'},
    {"trim_overhang",            no_argument,            0,        't'},
    {"output_nononverlapping_pair",no_argument,            0,        'i'},

    {"min_fraglen",              required_argument,      0,        'g'},
    {"max_fraglen",              required_argument,      0,        'G'},
    {"max_frac_prim_AS",         required_argument,      0,         'B'},
    {"min_passQ_frac",           required_argument,      0,        'Q'},
    {"max_N_filter",             required_argument,      0,        'y'},
    {"max_snv_filter",           required_argument,      0,        'x'},
    {"max_pair_mismatch_frac",   required_argument,      0,        'N'},
    {"clustered_mut_cutoff",     required_argument,      0,        'c'},
    {"filter_5endclip",          no_argument,            0,        '5'},

    {"consensus_mode",           required_argument ,     0,        'M'},
    {"dirtmp",                   required_argument ,     0,        'd'},
    {"thread",                   required_argument,      0,        'T'},
    {0,0,0,0}
};

const char* filter_short_options = "b:m:M:o:lCq:d:tT:ig:G:B:Q:y:x:c:5r:N:";

void filter_print_help()
{
  std::cerr<< "---------------------------------------------------\n";
  std::cerr<< "Usage: codec filter [options]\n";
  std::cerr<< "General Options:\n";
  std::cerr<< "-b/--bam,                              Input bam [required]\n";
  std::cerr<< "-r/--reference,                        reference sequence in fasta format [null].\n";
  std::cerr<< "-o/--outbam,                           Output unmapped bamfile [required].\n";
  std::cerr<< "-m/--mapq,                             Min mapping quality [10].\n";
  std::cerr<< "-q/--baseq,                            If one of the baseq < cutoff, make all baseq = 2 so that VC will ingnore them. [0].\n";
  std::cerr<< "-t/--trim_overhang,                    When perform paired-end consensus, if true then only do consensus of the overlapped region [false].\n";
  std::cerr<< "-C/--clip3,                            trim the 3'end soft clipping [false].\n";
//  std::cerr<< "-s/--output_singleend,                 The R1R2 consensus will be output in a single end format [false].\n";
//  std::cerr<< "-p/--pair_min_overlap,                 When using selector, the minimum overlap between the two ends of the pair [1].\n";
  std::cerr<< "-d/--dirtmp,                           Temporary dir for sorted bam [/tmp]\n";
  std::cerr<< "-T/--thread,                           Number of threads for sort [1]\n";
  std::cerr<< "-i/--output_nononverlapping_pair,        Allow output of non-overlaping pairs, usually caused by intermolecular ligation. This will simply print the original reads.  [false]\n";

  std::cerr<<"filtering options\n";
  std::cerr<< "-G/--max_fraglen,                      Filter out a read if its fragment length is larger than this value [INT_MAX].\n";
  std::cerr<< "-g/--min_fraglen,                      Filter out a read if its fragment length is less than this value [30].\n";
  std::cerr<< "-5/--filter_5endclip,                  Filtering out reads with 5'end soft clipping [False].\n";
  std::cerr<< "-B/--max_frac_prim_AS,                 Filter out a read if the AS of the secondary alignment is >= this fraction times the AS of the primary alignment [FLOAT_MAX].\n";
  std::cerr<< "-Q/--min_passQ_frac,                   Filter out a read if the fraction of bases passing quality threshold (together with -q) is less than this number [0].\n";
  std::cerr<< "-x/--max_snv_filter,                   Skip a read if the number of mismatch bases is larger than this value [INT_MAX].\n";
  std::cerr<< "-y/--max_N_filter,                     Skip a read if its num of N bases is larger than this value [INT_MAX].\n";
  std::cerr<< "-N/--max_pair_mismatch_frac,           Filter out a read-pair if its R1 and R2 has mismatch fraction larger than this value  [1.0].\n";
  std::cerr<< "-c/--clustered_mut_cutoff,             Filter out a read if at least this number of mutations occur in a window of 30 bp near the read end (<15 bp). [INT_MAX].\n";
  //std::cerr<< "-M/--consensus_mode,                   0 for paired end consensus. [0]\n";
}

int filter_parse_options(int argc, char* argv[], CssOptions& opt) {
  int option_index;
  int next_option = 0;
  do {
    next_option = getopt_long(argc, argv, filter_short_options, filter_long_options, &option_index);
    switch (next_option) {
      case -1:break;
      case 'r':
        opt.reference = optarg;
        break;
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
        opt.bqual_min = atoi(optarg);
        break;
      case 'l':
        opt.load_supplementary = true;
        break;
      case 't':
        opt.trim_overhang = true;
        break;
      case 'i':
        opt.output_nononverlapping_pair = true;
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
//      case 'p':
//        opt.pair_min_overlap = atoi(optarg);
//        break;
      case 'T':
        opt.thread = atoi(optarg);
        break;
      case 'g':
        opt.min_fraglen = atoi(optarg);
        break;
      case 'G':
        opt.max_fraglen = atoi(optarg);
        break;
      case '5':
        opt.filter_5endclip = true;
        break;
      case 'Q':
        opt.min_passQ_frac = atof(optarg);
        break;
      case 'B':
        opt.max_frac_prim_AS = atof(optarg);
        break;
      case 'y':
        opt.max_N_filter = atoi(optarg);
        break;
      case 'x':
        opt.max_snv_filter = atoi(optarg);
        break;
      case 'N':
        opt.max_pair_mismatch_frac = atof(optarg);
        break;
      case 'c':
        opt.clustered_mut_cutoff = atoi(optarg);
        break;
      default:filter_print_help();
        return 1;
    }
  } while (next_option != -1);

  return 0;
}


int codec_filter(int argc, char ** argv) {
  CssOptions opt;
  int parse_ret =  filter_parse_options(argc, argv, opt);
  if (parse_ret) return 1;
  if (argc == 1) {
    filter_print_help();
    return 1;
  }
  if (opt.outbam.empty()) {
    std::cerr << "-o/--outbam is required" << std::endl;
    filter_print_help();
    return 1;
  }
  if (opt.bam.empty()) {
    std::cerr << "-b/--bam is required" << std::endl;
    filter_print_help();
    return 1;
  }
//  if (opt.output_nononverlapping_pair) {
//    if (opt.pair_min_overlap != 0) {
//      std::cerr << "-p/--pair_min_overlap has to be 0 if -i/--output_nononverlapping_pair is true\n";
//      return 1;
//    }
//  }

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

  SeqLib::RefGenome ref;
  ref.LoadIndex(opt.reference);
  const int L = 250;
  auto errorstat = cpputil::ErrorStat(L, opt.bqual_min);
  cpputil::InsertSeqFactory isf(temp,
                                0,
                                opt.load_supplementary,
                                opt.load_secondary,
                                false,
                                !opt.load_unpair,
                                opt.clip3);
  cpputil::UnMappedBamWriter writer(opt.outbam, isf.bamheader());
  int64_t read_counter = 0;
  int64_t duplex_counter = 0;
  int64_t pass_counter = 0;
  if (opt.consensus_mode == 0) {
    while (!isf.finished()) {
      std::vector<cpputil::Segments> frag;
      while( (frag= isf.FetchReadNameSorted(opt.load_unpair)).size() > 0) {
        for (auto seg : frag) {
          assert (seg.size() == 2);
          ++ read_counter;
          int ol = cpputil::GetNumOverlapBasesPEAlignment(seg);
          if (ol > 0) {
            ++ duplex_counter;
            std::vector<std::string> orig_quals; // consensus output
            std::vector<std::string> seqs;
            for (auto&s : seg) {
              seqs.push_back(s.Sequence());
            }
            int pair_nmismatch = 0, olen = 0;
            float nqpass;
            int fail = cpputil::FailFilter(frag, isf.bamheader(), ref, opt, errorstat, true, pair_nmismatch, nqpass, olen);
            if (fail) {
              continue;
            }
            ++ pass_counter;
            auto seq = cpputil::PairConsensus(seg, seqs, opt.trim_overhang, opt.bqual_min, orig_quals); // trim overhang if true
            if (seg.front().FirstFlag()) {
              writer.WriteRecord(seg.front(), seg.back(), seq.first, seq.second, orig_quals.front(), orig_quals.back());
            }
            else if (seg.back().FirstFlag()) {
              writer.WriteRecord(seg.back(), seg.front(), seq.second, seq.first, orig_quals.back(), orig_quals.front());
            }
          } else if (opt.output_nononverlapping_pair) {
            writer.WriteRecord(seg.front(), seg.back());
          }
        }
      }
    }
  } else if (opt.consensus_mode == -1) { // dummy code for denovo consensus
    std::cerr << "denovo consensus yet to be build\n"; 
  } else {
    std::cerr << "Consensus mode should be either 0 or 1\n";
    return 1;
  }
  std::cout << "total molecule" << "\t" << "duplex molecule" << "\t" << "pass filter duplex" << "\t" << "filtered by mapq" \
            <<  "\t" << "filtered by pairmismatch rate"  << "\t" << "filtered by q30 rate" << "\t" << "filtered by edit" \
            << "\t" << "filtered by clustered" << "\t" << "filtered by sclip" << "\t" << "filtered by largefrag" << std::endl;

  std::cout << read_counter << "\t" << duplex_counter << "\t" << pass_counter << "\t" << errorstat.n_filtered_by_mapq << \
            "\t" << errorstat.n_filtered_pairmismatch_rate << "\t" << errorstat.n_filtered_q30rate << "\t" << \
            errorstat.n_filtered_edit << "\t" << errorstat.n_filtered_clustered << "\t" << errorstat.n_filtered_sclip << \
            "\t" << errorstat.n_filtered_largefrag  << std::endl;
  return 0;
}
