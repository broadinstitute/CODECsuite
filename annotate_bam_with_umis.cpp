//
// Created by Ruolin Liu on 7/2/20.
// In order to avoid sorting, the fastq is saved in the memory
// Therefore this program is designed to be fast and mem efficient for large hash (1B keys)
// The nucleotide is sotre as 3bits integer
//
#include <iostream>
#include <string>
#include <getopt.h>
#include <cassert>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <unordered_map>
#include "SeqLib/BamReader.h"
#include "SeqLib/BamWriter.h"

#include "SeqLib/FermiAssembler.h"
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
#include "DNAUtils.h"

using std::string;
struct Options {
  string i1;
  string i2;
  string input_bam;
  string output;
};

static struct option  long_options[] = {
    {"index1",                       required_argument,      0,        '1'},
    {"index2",                       required_argument,      0,        '2'},
    {"ibam",                         required_argument,      0,        'i'},
    {"obam",                         required_argument,      0,        'o'},
    {0,0,0,0}
};

const char*short_options = "1:2:o:i:";

void print_help()
{
  std::cerr<< "---------------------------------------------------\n";
  std::cerr<< "Usage: annotate_bam_with_umis -1 index.1.fq.gz -2 index.2.fq.gz -i input.bam -o output.bam\n";
  std::cerr<< "Important!! The program relies on the read_id appearing in the same order in the bam as in the index file.\n";
  std::cerr<< "Hyphenate index1 and index2 to a single UMI structure as AAA-BBB.\n";
  std::cerr<< "Then add the UMI to the unaligned bam as tag RX\n";
  std::cerr<< "The UMI fastq is at index1.fastq.gz and index2.fastq.gz.\n";
  std::cerr<< "The reads in the bam can be a strict subset of that in the index fastq, but not vise-versa.\n";
  std::cerr<< "-1/--index1,                               read1 UMI\n";
  std::cerr<< "-2/--index2,                               read2 UMI\n";
  std::cerr<< "-i/--ibam,                                 input input bam with UMI annotated.\n";
  std::cerr<< "-o/--obam,                                 output unmapped bam with UMI annotated.\n";
}

int parse_options(int argc, char* argv[], Options& opt) {
  int option_index;
  int next_option = 0;
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
    switch(next_option) {
      case -1:
        break;
      case '1':
        opt.i1 = optarg;
        break;
      case '2':
        opt.i2 = optarg;
        break;
      case 'o':
        opt.output = optarg;
        break;
      case 'i':
        opt.input_bam = optarg;
        break;
      default:
        print_help();
        return 1;
    }
  }while(next_option != -1);

  return 0;
}

int main(int argc, char** argv) {
  Options opt;
  int parse_ret = parse_options(argc, argv, opt);
  if (parse_ret) return 1;
  if (argc == 1) {
    print_help();
    exit(1);
  }

  SeqLib::BamReader br;
  SeqLib::BamRecordVector brv;
  SeqLib::BamRecord r;
  SeqLib::FermiAssembler fermi;
  br.Open(opt.input_bam);
  size_t count = 0;
  while (br.GetNextRecord(r) && count++ < 20000) {
    brv.push_back(r);
  }
  fermi.AddReads(brv);
  std::cout << "size: " << fermi.NumSequences() << std::endl;
  fermi.PerformAssembly();
  std::vector<std::string> contigs = fermi.GetContigs();
  for (size_t i = 0; i < contigs.size(); ++i) {
    std::cout << ">contig" << i << std::endl << contigs[i] << std::endl;
  }
  return 0;

  cpputil::FastxReader R1_reader(opt.i1);
  cpputil::FastxReader R2_reader(opt.i2);
  SeqLib::BamReader bam_reader;
  SeqLib::BamWriter bam_writer;
  bam_reader.Open(opt.input_bam);
  bam_writer.SetHeader(bam_reader.Header());
  bam_writer.Open(opt.output);
  bam_writer.WriteHeader();

  cpputil::FastxRecord read1;
  cpputil::FastxRecord read2;
  SeqLib::BamRecord b;
  uint32_t read_counter = 0;
  std::string seq;
  while (bam_reader.GetNextRecord(b)) {
    if (read_counter++ % 2 == 0) {
      while(R1_reader.yield(read1)) {
        R2_reader.yield(read2);
        if (read1.name() != read2.name()) {
          string errmsg = "read1 query name " + read1.id + " and read2 query name " + read2.id + " do not match!";
          throw std::runtime_error(errmsg);
        }
        if (read1.name() == b.Qname()) {
          seq = read1.seq + "-" + read2.seq;
          break;
        }
      };
    }
    if (read1.name() != b.Qname()) {
      string errmsg = "read query name " + read1.name() + " and bam record query name " + b.Qname() + " do not match!";
      throw std::runtime_error(errmsg);
    }
    b.AddZTag("RX", seq);
    bam_writer.WriteRecord(b);
    if (read_counter % 1000000 == 0) {
      printf("Processed %d reads\n", read_counter);
    }
  }
  printf("Finished processing %d reads\n", read_counter);

}
