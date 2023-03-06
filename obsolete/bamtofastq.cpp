//
// Created by Ruolin Liu on 10/6/20.
//

#include <iostream>
#include <getopt.h>
#include <cassert>
#include <cstdlib>
#include <unistd.h>
#include <string.h>
#include <tuple>
#include <SeqLib/BamReader.h>
#include <SeqLib/BamRecord.h>

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
#include "FastxIO.h"

using std::string;
struct Options {
  string bam;
  string fastq1;
  string fastq2;
  string tmpdir = "/tmp";
  int thread = 1;
};


static struct option  long_options[] = {
    {"input1",                  required_argument,      0,        '1'},
    {"input2",                  required_argument,      0,        '2'},
    {"bam",                     required_argument ,     0,        'b'},
    {"tmpdir",                  required_argument ,     0,        't'},
    {"thread",                  required_argument,      0,        'p'},
    {0,0,0,0}
};

const char*short_options = "1:2:b:t:p:";

void print_help()
{
  std::cerr<< "---------------------------------------------------\n";
  std::cerr<< "Usage: consensus [options]\n";
  std::cerr<< "General Options:\n";
  std::cerr<< "-b/--bam,                            Bam input\n";
  std::cerr<< "-1/--fastq1,                         Output Fastq1\n";
  std::cerr<< "-2/--fastq2,                         Output Fastq2\n";
  std::cerr<< "-t/--tmpdir,                         Temporary dir for sorted bam [/tmp]\n";
  std::cerr<< "-p/--thread,                         Number of threads for sort [1]\n";
}

int parse_options(int argc, char* argv[], Options& opt) {
  int option_index;
  int next_option = 0;
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
    switch (next_option) {
      case -1:break;
      case '1':
        opt.fastq1 = optarg;
        break;
      case '2':
        opt.fastq2 = optarg;
        break;
      case 'b':
        opt.bam = optarg;
        break;
      case 't':
        opt.tmpdir = optarg;
        break;
      case 'p':
        opt.thread = atoi(optarg);
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
  char temp[100];
  strcpy(temp, opt.tmpdir.c_str());
  strcat(temp, "/tempsort.XXXXXX");
  int fd = mkstemp(temp);
  if (fd == -1) {
    std::cerr << "unable to create temp file for sorting bam in queryname order\n";
    return 1;
  }
  string samsort = "samtools sort -n " + opt.bam + " -o " + string(temp) + " -@ " + std::to_string(opt.thread);
  std::cout << samsort << std::endl;
  std::system(samsort.c_str());
  SeqLib::BamRecord read;
  SeqLib::BamReader input;
  input.Open(temp);
  cpputil::FastqWriter R1(opt.fastq1);
  cpputil::FastqWriter R2(opt.fastq2);
  std::vector<SeqLib::BamRecord> pair(2);
  // After sorted by name, using yield like approach
  while(input.GetNextRecord(read)) {
    if (read.FirstFlag() && !read.SupplementaryFlag() && !read.SecondaryFlag()) {
      if (pair[0].isEmpty()) {
        pair[0] = read;
        if (!pair[1].isEmpty()) {
          if (pair[1].Qname() != pair[0].Qname()) {
            throw std::runtime_error("Bam file must be query name sorted! Exit at read " + read.Qname());
          } else {
            cpputil::FastxRecord fx1(pair[0], true);
            R1.Write(fx1);
            cpputil::FastxRecord fx2(pair[1], true);
            R2.Write(fx2);
            pair.clear();
            pair.resize(2);
          }
        }
      } else {
        throw std::runtime_error("Duplicated read name, " + read.Qname());
      }
    }

    if (!read.FirstFlag() && !read.SupplementaryFlag() && !read.SecondaryFlag()) {
      if (pair[1].isEmpty()) {
        pair[1] = read;
        if (!pair[0].isEmpty()) {
          if (pair[1].Qname() != pair[0].Qname()) {
            throw std::runtime_error("Bam file must be query name sorted! Exit at read " + read.Qname());
          } else {
            cpputil::FastxRecord fx1(pair[0], true);
            R1.Write(fx1);
            cpputil::FastxRecord fx2(pair[1], true);
            R2.Write(fx2);
            pair.clear();
            pair.resize(2);
          }
        }
      } else {
        throw std::runtime_error("Duplicated read name, " + read.Qname());
      }
    }
  }
  close(fd);
  unlink(temp);
  return 0;
}
