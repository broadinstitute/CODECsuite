//
// Created by Ruolin Liu on 3/26/20.
//

#ifndef CPPUTIL_INCLUDE_STRINTUTILS_H_
#define CPPUTIL_INCLUDE_STRINTUTILS_H_

#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <numeric>


namespace cpputil {
inline void split_by_char(const std::string &s, char c,
           std::vector<std::string> &v) {
  int i = 0;
  int j = s.find(c);

  while (j >= 0) {
    v.push_back(s.substr(i, j - i));
    i = ++j;
    j = s.find(c, j);

    if (j < 0) {
      v.push_back(s.substr(i, s.length()));
    }
  }
}

inline std::vector<std::string> split(const std::string& s, const std::string& delims)
{
  std::vector<std::string> result;
  std::string::size_type lastPos = s.find_first_not_of(delims, 0);
  std::string::size_type pos = s.find_first_of(delims, lastPos);
  while (std::string::npos != pos || std::string::npos != lastPos) {
    result.push_back(s.substr(lastPos, pos - lastPos));
    lastPos = s.find_first_not_of(delims, pos);
    pos = s.find_first_of(delims, lastPos);
  }
  return result;
}

/*
 * https://stackoverflow.com/questions/9277906/stdvector-to-string-with-custom-delimiter
 * By Shadow2531
 */
template <typename T>
std::string join(const T& v, const std::string& delim) {
  std::ostringstream s;
  for (const auto& i : v) {
    if (&i != &v[0]) {
      s << delim;
    }
    s << i;
  }
  return s.str();
}

inline double entropy(const std::string& dna_seq) {
  std::map<char, int> cnt;
  for (const char& d : dna_seq)  {
    cnt[d]++;
  }
  std::vector<int> vec;
  std::vector<double> p;
  int tot = 0;
  for (auto it : cnt) {
    tot += it.second;
    vec.push_back(it.second);
  }
  p.resize(vec.size());
  std::transform(vec.begin(), vec.end(), p.begin(), [&tot](double x) {return x / tot;});
  std::transform(p.begin(), p.end(), p.begin(), [](double x) {return x * log2(x);});
  return -std::accumulate(p.begin(), p.end(), 0.0);
}

}

#endif //CPPUTIL_INCLUDE_STRINTUTILS_H_
