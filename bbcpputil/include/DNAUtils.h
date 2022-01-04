//
// Created by Ruolin Liu on 12/20/19.
//

#ifndef REALIGN_INCLUDE_DNAUTILS_H_
#define REALIGN_INCLUDE_DNAUTILS_H_
#include <string>
#include <algorithm>
#include <cassert>

namespace cpputil {

inline char complement(char n) {
  switch (n) {
    case 'A':return 'T';
    case 'T':return 'A';
    case 'G':return 'C';
    case 'C':return 'G';
    case 'N':return 'N';
    case 'n':return 'n';
    case 'a':return 't';
    case 't':return 'a';
    case 'c':return 'g';
    case 'g':return 'c';
  }
  assert(false);
  return ' ';
}

inline std::string complementString(std::string x) {
  std::transform(std::begin(x), std::end(x), std::begin(x), complement);
  return x;
}

inline void reverse_complement(std::string &seq) {
  std::transform(std::begin(seq), std::end(seq), std::begin(seq), complement);
  for (int i = 0, j = seq.size() - 1; i < j; i++, j--) {
    std::swap(seq[i], seq[j]);
  }
}

inline void reverse(std::string &seq) {
  for (int i = 0, j = seq.size() - 1; i < j; i++, j--) {
    std::swap(seq[i], seq[j]);
  }
}

inline void PrintQualString(const std::string& qual, int min_bq = 20, int offset = 33) {
  std::string line1;
  std::string line2;
  std::string stat;
  for (unsigned i = 0; i < qual.size(); ++i) {
    int q = (int) qual[i] - offset;
    int div = q / 10;
    int reminder = q % 10;
    line1 += std::to_string(div);
    line2 += std::to_string(reminder);
    stat += q >= min_bq ? "*" : " ";
  }
  std::cout << line1 << std::endl;
  std::cout << line2 << std::endl;
  std::cout << stat << std::endl;
}

inline int TrimLowBQfromBack(const std::string& qual, char bq) {
  int i = qual.size();
  while (i > 0 && qual[i-1] < bq) {i--;}
  return i;
}

inline int LastNfromBack(const std::string& seq) {
  //one pass the first not N from the back
  int i = seq.size();
  while (i > 0 && seq[i-1] == 'N') {i--;}
  return i;
}

inline int FirstNotNfromFront(const std::string& seq) {
  int i = 0;
  while (i < (int) seq.size() && seq[i] == 'N') {i++;}
  return i;
}

}
#endif //REALIGN_INCLUDE_DNAUTILS_H_
