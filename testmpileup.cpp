//
// Created by Ruolin Liu on 10/28/21.
//
#include "pileup.h"
#include "Algo.h"

int main(int argc, char **argv) {
//  cpputil::PileHandler ph;
//  ph = cpputil::PileHandler(argv[1], 0);
//  if (cpputil::ScanAllele(&ph, argv[2], std::stoi(argv[3]), argv[4][0], false) < 0)
//    return EXIT_FAILURE;
//  return EXIT_SUCCESS;
  int a, b;
  std::cout << cpputil::largest_cluster({10,13,30, 60, 90, 100, 103, 119}, 30, a, b) << std::endl;
  std::cout << a <<", " << b << std::endl;

  std::cout << cpputil::largest_cluster({10}, 30, a, b) << std::endl;
  std::cout << a <<", " << b << std::endl;

  std::cout << cpputil::largest_cluster({10, 60}, 30, a, b) << std::endl;
  std::cout << a <<", " << b << std::endl;
}

