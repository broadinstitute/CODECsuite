#include <iostream>
#include <string>
#include <getopt.h>
#include <cassert>
#include <algorithm>
#include <seqan/align.h>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include "Adapter.h"
#include "BamIO.h"

using FastxRecord = cpputil::ExtFastxRecord;
struct Options {
  std::string r1;
  std::string r2;
  std::string prefix = "/tmp/cds";
  std::string linker_type = "adapter_v2";
  std::string rgid = "A";
  std::string rgsm = "A";
  bool split_bam_output = false;
  int debug = 0;
  const unsigned MIN_READ_LEN = 15; // this has to be large than TRIM_5 + TRIM_3
  int nbase_skip = 0;
  unsigned JUNCTION_ADAPTER_MIN_LEN= 30;
  unsigned THREE_END_ADAPTER_MIN_LEN = 3;
  int POST_TRIM5 = 0;
  int POST_TRIM3 = 0;
  int R1_UMI_LEN = 0;
  int R2_UMI_LEN = 0;
  int INTERMOL_SIZE_DIFF = 5;
  //Shouldn't be below -50
  int GAP_SCORE = -8;
  int MISMATCH_SCORE = -4;
  int MATCH_SCORE = 1;
  int MIN_BQ_FROM_BACK = 3;
};

int HIGH_CONF = 0;
int LOW_CONF = 0;
int READ1_UNTRIMMED = 0;
int READ2_UNTRIMMED = 0;
int BOTH_UNTRIMMED = 0;
int LOST_BOTH = 0;
int LOST_READ1 = 0;
int LOST_READ2 = 0;
int SINGLE_INSERT_HIGHCONF = 0;
int SINGLE_INSERT_LOWCONF = 0;
int DOUBLE_LIGATION = 0;

static struct option  long_options[] = {
  {"R1",                       required_argument,      0,        '1'},
  {"R2",                       required_argument,      0,        '2'},
  {"prefix",                   required_argument,      0,        'o'},
  {"rgid",                     required_argument,      0,        'i'},
  {"rgsm",                     required_argument,      0,        's'},
  {"trim_full_linker_len",     required_argument,      0,        'F'},
  {"trim_prefix_linker_len",   required_argument,      0,        'P'},
  {"match",                    required_argument,      0,        'A'},
  {"nbase_skip",               required_argument,      0,        'n'},
  {"mismatch",                 required_argument,      0,        'B'},
  {"gap",                      required_argument,      0,        'G'},
  {"min_bq_from_back",         required_argument,      0,        'T'},
  {"r2_umi_len",               required_argument,      0,        'u'},
  {"r2_umi_len",               required_argument,      0,        'U'},
  {"post_trim5",               required_argument,      0,        'f'},
  {"post_trim3",               required_argument,      0,        't'},
  {"debug",                    required_argument,      0,        'd'},
  {"split_bam_output",         no_argument,      0,        'S'},
  //hidden parameter
  //{"linker_type",              required_argument,      0,        'l'},
  {0,0,0,0}
};

const char* trim_short_options = "F:P:1:2:A:B:G:o:t:f:d:l:u:U:T:i:s:n:S";


void codec_trim_usage()
{
  std::cerr<< "---------------------------------------------------\n";
  std::cerr<< "Usage: codec trim [options]\n";
  std::cerr<< "General Options:\n";
  std::cerr<< "-1/--R1,                               R1.fastq. [required]\n";
  std::cerr<< "-2/--R2,                               R1.fastq. [required]\n";
  std::cerr<< "-o/--prefix,                           Output prefix. [/tmp/cds.]\n";
  std::cerr<< "-n/--nbase_skip,                       Do not trim in the first n base of the read. [n=0]\n";
  std::cerr<< "-i/--rgid,                             Read group id when output bam [A]\n";
  std::cerr<< "-s/--rgsm,                             Read group sample when output bam [A]\n";
  //std::cerr<< "-l/--linkertype,                     One of those, medium_v1, medium_v2, long_v1, adapter_v2. [hidden]\n";
  std::cerr<< "-T/--min_bq_from_back,                 After adpters being removed, trim bases with bq less than this value from the back [3].\n";
  std::cerr<< "-u/--r1_umi_len,                       num of umi bases to be trimmed from 5'end of read1 [0].\n";
  std::cerr<< "-U/--r2_umi_len,                       num of umi bases to be trimmed from 5'end of read2 [0].\n";
  std::cerr<< "-F/--trim_full_linker_len,             Minimum trimming length if linker is found in the middle of a read [30].\n";
  std::cerr<< "-P/--trim_prefix_linker_len,           Minimum trimming length if linker is found at the 3'end of a read [3]\n";
  std::cerr<< "-t/--post_trim3,                       Num of bases to be trimmed at the 3'end of a read after adapter removed [0]\n";
  std::cerr<< "-f/--post_trim5,                       Num of bases to be trimmed at the 5'end of a read after adapter removed [0]\n";
  std::cerr<< "-A/--match                             Score for a sequence match [1]\n";
  std::cerr<< "-B/--mismatch                          Penalty for a mismatch [4]\n";
  std::cerr<< "-G/--gap                               Penalty for open or extend a gap [5]\n";
  std::cerr<< "-d/--debug,                            1: detail pairwise alignment plot, 2: Qscore plot, default no debug[0]\n";
  std::cerr<< "-S/--split_bam_output,                     Output byproduct to separate bams [false]\n";
}

int trim_parse_options(int argc, char* argv[], Options& opt) {
  int option_index;
  int next_option = 0;
  do {
      next_option = getopt_long(argc, argv, trim_short_options, long_options, &option_index);
      switch(next_option) {
        case -1:
          break;
        case '1':
          opt.r1 = optarg;
          break;
        case '2':
          opt.r2 = optarg;
          break;
        case 'o':
          opt.prefix = optarg;
          break;
        case 'i':
          opt.rgid = optarg;
          break;
        case 's':
          opt.rgsm = optarg;
          break;
        case 'd':
          opt.debug = atoi(optarg);
          break;
        case 'n':
          opt.nbase_skip = atoi(optarg);
          break;
        case 'P':
          opt.THREE_END_ADAPTER_MIN_LEN = atoi(optarg);
          break;
        case 'F':
          opt.JUNCTION_ADAPTER_MIN_LEN = atoi(optarg);
          break;
        case 'A':
          opt.MATCH_SCORE = atoi(optarg);
          break;
        case 'B':
          opt.MISMATCH_SCORE = -atoi(optarg);
          break;
        case 'G':
          opt.GAP_SCORE = -atoi(optarg);
          break;
        case 'l': // hidden paramter
          opt.linker_type = optarg;
          break;
        case 'f':
          opt.POST_TRIM5 = atoi(optarg);
          break;
        case 't':
          opt.POST_TRIM3 = atoi(optarg);
          break;
        case 'T':
          opt.MIN_BQ_FROM_BACK = atoi(optarg);
          break;
        case 'u':
          opt.R1_UMI_LEN = atoi(optarg);
          break;
        case 'U':
          opt.R2_UMI_LEN = atoi(optarg);
          break;
        case 'S':
          opt.split_bam_output = true;
          break;
        default:
          codec_trim_usage();
          return 1;
      }
   }while(next_option != -1); 
    
   return 0;
}



int codec_trim(int argc, char** argv)
{
  Options opt;
  int parse_ret =  trim_parse_options(argc, argv, opt);
  if (parse_ret) return 1;
  if (argc == 1) {
     codec_trim_usage();
     exit(0);
  }

  //Input
  cpputil::FastxReader R1_reader(opt.r1);
  cpputil::FastxReader R2_reader(opt.r2);
  FastxRecord read1;
  FastxRecord read2;

  //Output
  cpputil::UnMappedBamWriter highconf(opt.prefix + ".trim.bam", opt.rgid, opt.rgsm);

  cpputil::UnMappedBamWriter lost;
  cpputil::UnMappedBamWriter singleton;
  // obsolete
  //cpputil::UnMappedBamWriter lowconf;
  //cpputil::UnMappedBamWriter untrimboth;
  if (opt.split_bam_output) {
    lost.Open(opt.prefix + ".lost.bam", opt.rgid, opt.rgsm);
    singleton.Open(opt.prefix + ".singleton.bam", opt.rgid, opt.rgsm);
    //lowconf.Open(opt.prefix + ".lowconf.bam", opt.rgid, opt.rgsm);
    //untrimboth.Open(opt.prefix + ".untrimboth.bam", opt.rgid, opt.rgsm);
  }
  //cpputil::UnMappedBamWriter trimone(opt.prefix + ".trimone.bam", opt.rgid, opt.rgsm);
  //cpputil::UnMappedBamWriter single_highconf(opt.prefix + ".singleinsert.bam", opt.rgid, opt.rgsm);
  //cpputil::UnMappedBamWriter single_lowconf(opt.prefix + ".single_lowconf.bam", opt.rgid, opt.rgsm);

  int num_debug_output = 0;
  while (R1_reader.yield(read1)) {
    R2_reader.yield(read2);
    if (read1.name() != read2.name()) {
      std::cerr << "read 1 name " << read1.name() << " and read 2 name " << read2.name() << " do not match!\n";
    }
    CDS::FragTrimS fts = CDS::TrimPairedRead(opt, read1, read2);
    read1.tm = fts.first.tm;
    read2.tm = fts.second.tm;
    // output fastq and stat
    if (opt.debug == 1) {
      std::cerr << CDS::ToString(fts.first.tm) << ", " << CDS::ToString(fts.second.tm) << std::endl;
    }
    if (read1.seq.size() < opt.MIN_READ_LEN  && read2.seq.size() < opt.MIN_READ_LEN) {
      LOST_BOTH ++;
      if (opt.split_bam_output) {
        lost.WriteRecord(read1, true);
        lost.WriteRecord(read2, false);
      }
    }
    else if (read1.seq.size() < opt.MIN_READ_LEN) {
      LOST_READ1 ++;
      if (opt.split_bam_output) {
        singleton.WriteRecord(read1, true);
        singleton.WriteRecord(read2, false);
      } else {
        if (read1.seq.size() == 0) {
          read1.seq = "ATG"; // fake a read
          read1.qual = "!!!";
        }
        highconf.WriteRecord(read1, true);
        highconf.WriteRecord(read2, false);
      }
    } else if (read2.seq.size() < opt.MIN_READ_LEN) {
      LOST_READ2 ++;
      if (opt.split_bam_output) {
        singleton.WriteRecord(read1, true);
        singleton.WriteRecord(read2, false);
      } else {
        if (read2.seq.size() == 0) {
          read2.seq = "ATG"; // fake a read
          read2.qual = "!!!";
        }
        highconf.WriteRecord(read1, true);
        highconf.WriteRecord(read2, false);
      }
    }
    else {
      if (opt.linker_type == "adapter_v2")  {
        //assert(fts.first.tm == CDS::TRIM && fts.second.tm == CDS::TRIM);
        if (fts.first.tm == CDS::TRIM && fts.second.tm == CDS::TRIM) {
          HIGH_CONF++;
          highconf.WriteRecord(read1, true);
          highconf.WriteRecord(read2, false);
        } else {
          LOW_CONF++;
          //No longer have low conf
//          if (opt.split_bam_output) {
//            lowconf.WriteRecord(read1, true);
//            lowconf.WriteRecord(read1, false);
//          }
        }
      }
      else { // not used, keep for legacy reason
        if (fts.second.tm == CDS::TRIM_INSUF &&
            fts.first.tm == CDS::TRIM_INSUF ) {
          BOTH_UNTRIMMED++;
          //untrimboth.WriteRecord(read1, true);
          //untrimboth.WriteRecord(read2, false);
        }
        else if (fts.first.tm == CDS::DOUBLE_LIGATION || fts.second.tm == CDS::DOUBLE_LIGATION) {
          ++DOUBLE_LIGATION;
          //single_highconf.WriteRecord(read1, true);
          //single_highconf.WriteRecord(read2, false);
        }
        else if (fts.first.tm == CDS::TRIM && fts.second.tm == CDS::TRIM)
        {
          HIGH_CONF++;
          //highconf.WriteRecord(read1, true);
          //highconf.WriteRecord(read2, false);
        }
        else if (fts.first.tm == CDS::SINGLE && fts.second.tm == CDS::SINGLE) {
          SINGLE_INSERT_HIGHCONF++;
          //single_highconf.WriteRecord(read1, true);
          //single_highconf.WriteRecord(read2, false);
        }
        else if (fts.first.tm == CDS::SINGLE || fts.second.tm == CDS::SINGLE) {
          SINGLE_INSERT_LOWCONF++;
          //single_lowconf.WriteRecord(read1, true);
          //single_lowconf.WriteRecord(read2, false);
        }
        else if (fts.first.tm == CDS::TRIM || fts.first.tm == CDS::TRIM13
            || fts.second.tm == CDS::TRIM || fts.second.tm == CDS::TRIM13) {
          //else if (fts.first.tm == TRIM || fts.second.tm == TRIM) {
          LOW_CONF++;
          //lowconf.WriteRecord(read1, true);
          //lowconf.WriteRecord(read2, false);
        }
      }
    }
  }

  std::cout << "LOST_BOTH: " << LOST_BOTH << std::endl;
  std::cout << "LOST_READ1: " << LOST_READ1 << std::endl;
  std::cout << "LOST_READ2: " << LOST_READ2 << std::endl;
  std::cout << "HIGH_CONF: " << HIGH_CONF << std::endl;
  std::cout << "LOW_CONF: " << LOW_CONF << std::endl;
  std::cout << "BOTH_UNTRIMMED: " << BOTH_UNTRIMMED << std::endl;
  std::cout << "READ1_UNTRIMMED: " << READ1_UNTRIMMED << std::endl;
  std::cout << "READ2_UNTRIMMED: " << READ2_UNTRIMMED << std::endl;
  std::cout << "SINGLE_INSERT: " << SINGLE_INSERT_HIGHCONF << std::endl;
  std::cout << "SINGLE_INSERT_LOWCONF: " << SINGLE_INSERT_LOWCONF << std::endl;
  std::cout << "DOUBLE_LIGATION: " << DOUBLE_LIGATION << std::endl;
  std::cout << "TOTAL: " << LOST_BOTH + LOST_READ1 + LOST_READ2 + DOUBLE_LIGATION +\
                         HIGH_CONF + LOW_CONF + BOTH_UNTRIMMED + \
                         READ2_UNTRIMMED + READ1_UNTRIMMED + \
                         SINGLE_INSERT_LOWCONF + SINGLE_INSERT_HIGHCONF << std::endl;

  return 0;
}
