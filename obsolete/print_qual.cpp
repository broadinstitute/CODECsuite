//
// Created by Ruolin Liu on 7/12/20.
//

#include <iostream>
#include "DNAUtils.h"

int main(int argc, char** argv) {
  std::string qual = argv[1];
  int bq = std::stoi(argv[2]);
  cpputil::PrintQualString(qual, bq);
  return 0;
}