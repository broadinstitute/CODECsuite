//
// Created by Ruolin Liu on 10/31/21.
//

#ifndef CODECSUITE_BBCPPUTIL_INCLUDE_ALGO_H_
#define CODECSUITE_BBCPPUTIL_INCLUDE_ALGO_H_
#include <unordered_set>
#include <queue>
#include <string>

namespace cpputil {

class UniqueQueue {
  /*
   * Pair up Queue with Set so that Queue has only unique element
   */
public:
  UniqueQueue(size_t cap) : capacity_(cap), n_added_(0) {};
  UniqueQueue() : UniqueQueue(0) {};

  bool exist(const std::string& in) const {
    //auto h = hash_string(in.c_str());
    if (s_.find(in) == s_.end()) return false;
    else return true;
  }

  void add(const std::string& in) {
    //auto h = hash_string(in.c_str());
    if (s_.find(in) == s_.end()) {
      if (q_.size() < capacity_) {
        q_.push(in);
        s_.insert(in);
      } else {
       auto key = q_.front();
       q_.pop();
       s_.erase(key);
       q_.push(in);
       s_.insert(in);
      }
      ++n_added_;
    }
  }

  void clearQueue() {
    std::queue<std::string> empty;
    std::swap(q_, empty);
    s_.clear();
  }

  uint64_t NumAdded() const {
    return n_added_;
  }

private:
  // This is FNV-1, see http://en.wikipedia.org/wiki/Fowler_Noll_Vo_hash
  inline uint64_t hash_string(const char* __s) const
  {
    uint64_t hash = 0xcbf29ce484222325ull;
    for ( ; *__s; ++__s)
    {
      hash *= 1099511628211ull;
      hash ^= *__s;
    }
    return hash;
  }

  std::unordered_set<std::string> s_;
  std::queue<std::string> q_;
  const size_t capacity_;
  uint64_t n_added_;
};

inline int largest_cluster(std::vector<int> sortedpos, int window, int &beg, int &end) {
  std::sort(sortedpos.begin(), sortedpos.end());
  if (sortedpos.size() == 0) return 0;
  if (sortedpos.size() == 1) {
    beg = sortedpos[0];
    end = sortedpos[0];
    return 1;
  }
  int max_size = 1;
  for (unsigned i= 0; i < sortedpos.size() - 1; ++i) {
    unsigned j = i+1;
    int s = 1;
    for (; j < sortedpos.size();) {
      if (sortedpos[j] - sortedpos[i] < window) {
        ++s;
        ++j;
      } else {
        if (s > max_size) {
          beg = sortedpos[i];
          end = sortedpos[j - 1];
          max_size = s;
        }
        break;
      }
      if (j == sortedpos.size() && s > max_size) {
        beg = sortedpos[i];
        end = sortedpos[j - 1];
        max_size = s;
      }
    }
  }
  return max_size;
}

}

#endif //CODECSUITE_BBCPPUTIL_INCLUDE_ALGO_H_
