//
// Created by Ruolin Liu on 10/6/20.
//

#include <iostream>
#include <getopt.h>
#include <cassert>
#include <FastxIO.h>

using std::string;
struct Options {
  string fastq1;
  string fastq2;
  string index1;
  string index2;
  string out1;
  string out2;
};


static struct option  long_options[] = {
    {"input1",                  required_argument,      0,        '1'},
    {"input2",                  required_argument,      0,        '2'},
    {"index1",                  required_argument ,     0,        'i'},
    {"index2",                  required_argument ,     0,        'I'},
    {"out1",                    required_argument ,     0,        'o'},
    {"out2",                    required_argument ,     0,        'O'},
    {0,0,0,0}
};

const char*short_options = "1:2:o:O:i:I:";

void print_help()
{
  std::cerr<< "---------------------------------------------------\n";
  std::cerr<< "Usage: consensus [options]\n";
  std::cerr<< "General Options:\n";
  std::cerr<< "-1/--input1,                         Input Read 1\n";
  std::cerr<< "-2/--iunput2,                        Input Read 2\n";
  std::cerr<< "-i/--index1,                         index/UMI file for Read 1\n";
  std::cerr<< "-I/--index2,                         index/UMI file for Read 2\n";
  std::cerr<< "-o/--output1,                        Output Read 1\n";
  std::cerr<< "-o/--output2,                        Output Read 2\n";
}

int parse_options(int argc, char* argv[], Options& opt) {
  int option_index;
  int next_option = 0;
  do {
    next_option = getopt_long(argc, argv, short_options, long_options, &option_index);
    switch (next_option) {
      case -1:break;
      case 'i':
        opt.index1 = optarg;
        break;
      case 'I':
        opt.index2 = optarg;
        break;
      case '1':
        opt.fastq1 = optarg;
        break;
      case '2':
        opt.fastq2 = optarg;
        break;
      case 'o':
        opt.out1 = optarg;
        break;
      case 'O':
        opt.out2= optarg;
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

  cpputil::FastxReader R1_reader(opt.fastq1);
  cpputil::FastxReader R2_reader(opt.fastq2);
  cpputil::FastxReader I1_reader(opt.index1);
  cpputil::FastxReader I2_reader(opt.index2);
  cpputil::FastqWriter O1_writer(opt.out1);
  cpputil::FastqWriter O2_writer(opt.out2);
  cpputil::FastxRecord read1;
  cpputil::FastxRecord read2;
  cpputil::FastxRecord index1;
  cpputil::FastxRecord index2;
  cpputil::FastxRecord outread1;
  cpputil::FastxRecord outread2;
  while (R1_reader.yield(read1)) {
    R2_reader.yield(read2);
    I1_reader.yield(index1);
    I2_reader.yield(index2);
    assert(read1.name() == read2.name());
    assert(index1.name() == index2.name());
    assert(index1.name() == read1.name());
    outread1.id = read1.id;
    outread1.seq = index1.seq + read1.seq;
    outread1.qual = index1.qual + read1.qual;
    outread2.id = read2.id;
    outread2.seq = index2.seq + read2.seq;
    outread2.qual = index2.qual + read2.qual;
    O1_writer.Write(outread1);
    O2_writer.Write(outread2);
  }
  return 0;
}
