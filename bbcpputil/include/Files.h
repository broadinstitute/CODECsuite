//
// Created by Ruolin Liu on 5/11/21.
//
#ifndef CPPUTIL_INCLUDE_FILES_H_
#define CPPUTIL_INCLUDE_FILES_H_
#include <string>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>

namespace cpputil {

inline bool FileExist(const std::string &name) {
  std::ifstream f(name.c_str());
  return f.good();
}
}

#endif //CPPUTIL_INCLUDE_FILES_H_
