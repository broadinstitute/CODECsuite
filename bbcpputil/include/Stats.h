//
// Created by Ruolin Liu on 7/1/21.
//

#ifndef BBCPPUTIL_INCLUDE_STATS_H_
#define BBCPPUTIL_INCLUDE_STATS_H_
#include <vector>
#include <algorithm>
namespace cpputil {

inline int GetMode(const std::vector<int>& array) {
  /*
   * Return the mode of an integer vector. Works for short and medium size array and
   * the span of the array is not too larger
   *
   */
  assert(not array.empty());
  const auto ret = std::minmax_element(begin(array), end(array));
  std::vector<int> hist(*ret.second - *ret.first + 1);
  for (const auto& ii: array) {
    ++hist[ii - *ret.first];
  }
  return std::max_element(hist.begin(), hist.end()) - hist.begin() + *ret.first;
}


}
#endif //BBCPPUTIL_INCLUDE_STATS_H_
