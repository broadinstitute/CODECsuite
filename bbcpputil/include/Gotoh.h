//
// Created by Ruolin Liu on 9/27/20.
//

#ifndef ADAPTERTRIM_CPPUTIL_INCLUDE_GOTOH_H_
#define ADAPTERTRIM_CPPUTIL_INCLUDE_GOTOH_H_

/*
This code snippet finds all optimal paths using affine gap penalty.
So many public codes and textbooks implement or teach solutions leading to suboptimal solutions.
Inspired by paper https://www.biorxiv.org/content/10.1101/031500v1.full.pdf,
I decided to implement the correct algorithm which is documented at
https://www.researchgate.net/publication/19580571_Optimal_sequence_alignment_using_affine_gap_costs
*/

#include <iostream>
#include <string>
#include <vector>
#include <climits>
#include <array>
#include <assert.h>
#include <stack>

class Alignment {
  std::vector<std::string> display;
  int nm = 0;
public:
  struct Node {
    int i;
    int j;
    Node(int ii, int jj) : i(ii), j(jj) {}
    Node() = default;
  };

  Alignment(const std::string& row, const std::string col, const std::vector<Node>& path) {
    std::string rowgap;
    std::string colgap;
    std::string visual;
    for (unsigned ii = 0; ii < path.size() - 1; ++ii) {
      const Node & cur = path[ii];
      const Node & next = path[ii + 1];
      if (cur.i == next.i + 1 && cur.j == next.j + 1) {
        rowgap += row[next.i];
        colgap += col[next.j];
        if (row[next.i] != col[next.j]) {
          nm++;
          visual += ' ';
        } else {
          visual += '|';
        }
      } else if (cur.i == next.i && cur.j == next.j + 1) {
        nm ++;
        visual += ' ';
        rowgap += '-';
        colgap += col[next.j];
      } else if (cur.i == next.i + 1 && cur.j == next.j) {
        nm ++;
        visual += ' ';
        colgap += '-';
        rowgap += row[next.j];
      }
    }
    std::reverse(colgap.begin(), colgap.end());
    std::reverse(visual.begin(), visual.end());
    std::reverse(rowgap.begin(), rowgap.end());
    display.push_back(colgap);
    display.push_back(visual);
    display.push_back(rowgap);
  }

  int NM() const {
    return nm;
  }
  friend std::ostream& operator<<(std::ostream&, const Alignment&);
};

std::ostream& operator<<(std::ostream& os, const Alignment& align) {
  os << align.nm << '\n';
  for (const std::string& s : align.display) {
    os << s << '\n';
  }
  return os;
}

using std::vector;
class AffineGap {
  using Node =  Alignment::Node;
  /*
   * matrix as Reference at row top Query at column left
   */
  std::string rowstring_;
  std::string colstring_;
  vector<Node> align_path_;
  vector<vector<Node>> paths_;
  int nrow_;
  int ncol_;
  vector<vector<int>> R_; // match matrix
  vector<vector<int>> P_; // vertical insertion matrix
  vector<vector<int>> Q_; // horizontal deletion matrix
  vector<vector<bool>> vert_whole_; // a
  vector<vector<bool>> hori_whole_; // b
  vector<vector<bool>> diag_whole_; // c
  vector<vector<bool>> vert_top_half_; // d
  vector<vector<bool>> vert_bottom_half_; // e
  vector<vector<bool>> hori_left_half_; // f
  vector<vector<bool>> hori_right_half_; // g

  //gap score = gap_open + gap_ext * gap_len
  const static int gap_open_ = -1;
  const static int gap_ext_ = -1;
  static int DiagScore_(const char& a, const char& b) {
    return a == b? 0: -1;
  }

 public:
  AffineGap(const std::string& query, const std::string& ref): rowstring_(query), colstring_(ref), //align_path_(colstring_, rowstring_),
                                                     nrow_((int) rowstring_.length() + 1),
                                                     ncol_((int) colstring_.length() + 1),
                                                     R_(vector<vector<int>>(nrow_, vector<int>(ncol_))),
                                                     P_(vector<vector<int>>(nrow_, vector<int>(ncol_))),
                                                     Q_(vector<vector<int>>(nrow_, vector<int>(ncol_))),
                                                     vert_whole_(vector<vector<bool>>(nrow_ + 1, vector<bool>(ncol_ + 1))),
                                                     hori_whole_(vector<vector<bool>>(nrow_ + 1, vector<bool>(ncol_ + 1))),
                                                     diag_whole_(vector<vector<bool>>(nrow_ + 1, vector<bool>(ncol_ + 1))),
                                                     vert_top_half_(vector<vector<bool>>(nrow_ + 1, vector<bool>(ncol_ + 1))),
                                                     vert_bottom_half_(vector<vector<bool>>(nrow_ + 1, vector<bool>(ncol_ + 1))),
                                                     hori_left_half_(vector<vector<bool>>(nrow_ + 1, vector<bool>(ncol_ + 1))),
                                                     hori_right_half_(vector<vector<bool>>(nrow_ + 1, vector<bool>(ncol_ + 1)))
  {
    // init
    for (int j = 0; j < ncol_; ++j) {
      P_[0][j] = 2 * gap_open_ + std::max(ncol_, nrow_) * gap_ext_ - 1; // ensure a large number
      R_[0][j] = gap_open_ + j * gap_ext_;
    }
    for (int i = 0; i < nrow_; ++i) {
      Q_[i][0] = 2 * gap_open_ + std::max(ncol_, nrow_) * gap_ext_ - 1; // ensure a large number
      R_[i][0] = gap_open_ + i * gap_ext_;
    }
    R_[0][0] = 0;
    diag_whole_[nrow_][ncol_] = 1;

    //cost assignment
    for (int i = 0; i < nrow_; ++i) {
      for (int j = 0; j < ncol_; ++j) {
        if (i != 0) {
          P_[i][j] = gap_ext_ + std::max(P_[i-1][j], R_[i-1][j] + gap_open_);
          if (P_[i][j] == gap_ext_ + P_[i-1][j]) vert_top_half_[i-1][j] = 1;
          if (P_[i][j] == gap_ext_ + gap_open_ + R_[i-1][j]) vert_bottom_half_[i-1][j] = 1;
        }
        if (j != 0) {
          Q_[i][j] = gap_ext_ + std::max(Q_[i][j-1], R_[i][j-1] + gap_open_);
          if (Q_[i][j] == gap_ext_ + Q_[i][j-1]) hori_left_half_[i][j-1] = 1;
          if (Q_[i][j] == gap_ext_ + gap_open_ + R_[i][j-1]) hori_right_half_[i][j-1] = 1;
        }
        if (i != 0 && j != 0 ) {
          R_[i][j] = std::max(R_[i-1][j-1] + DiagScore_(colstring_[j-1], rowstring_[i-1]), std::max(Q_[i][j], P_[i][j]));
          if (R_[i][j]  == R_[i-1][j-1] + DiagScore_(colstring_[j-1], rowstring_[i-1])) diag_whole_[i][j] = 1;
        }
        if (R_[i][j]  == P_[i][j]) vert_whole_[i][j] = 1;
        if (R_[i][j]  == Q_[i][j]) hori_whole_[i][j] = 1;
      }
    }
//    std::cout<< "cost assignment\t";
//    std::cout << diag_whole_[nrow_-1][ncol_-1] << "\t" << vert_whole_[nrow_-1][ncol_-1] << "\t" << hori_whole_[nrow_-1][ncol_-1] << "\n";

    //edge assignment
    for (int i = nrow_ - 1; i >= 0; --i) {
      for (int j = ncol_ - 1; j >= 0; --j) {
        if ((vert_whole_[i+1][j] == 0 || vert_bottom_half_[i][j] == 0) &&
            (hori_whole_[i][j+1] == 0 || hori_right_half_[i][j] == 0) &&
            diag_whole_[i+1][j+1] == 0) {
//          if (vert_bottom_half_[i][j] == 0 && vert_whole_[i+1][j] != 0) {
//            std::cerr << "vi: " << i << " vj: " << j << std::endl;
//          }
//          if (hori_whole_[i][j+1] != 0 && hori_right_half_[i][j] == 1) {
//            std::cerr << "hi: " << i << " hj: " << j << std::endl;
//          }
          vert_whole_[i][j] = 0;
          hori_whole_[i][j] = 0;
          diag_whole_[i][j] = 0;
        }
        if (vert_whole_[i+1][j] == 0  &&
            hori_whole_[i][j+1] == 0  &&
            diag_whole_[i+1][j+1] == 0) {
          continue;
        } else {
          if ( vert_whole_[i+1][j] == 1 && vert_top_half_[i][j] == 1) {
            vert_top_half_[i+1][j] = 1 -  vert_bottom_half_[i][j];
            vert_bottom_half_[i][j] = 1 - vert_whole_[i][j];
            vert_whole_[i][j] = 1;
          } else {
            vert_top_half_[i+1][j] = 0;
            vert_bottom_half_[i][j] = 0;
          }

          if ( hori_whole_[i][j + 1] == 1 && hori_left_half_[i][j] == 1) {
            hori_left_half_[i][j+1] = 1 -  hori_right_half_[i][j];
            hori_right_half_[i][j] = 1 - hori_whole_[i][j];
            hori_whole_[i][j] = 1;
          } else {
            hori_left_half_[i][j+1] = 0;
            hori_right_half_[i][j] = 0;
          }
        }
      }
    }
    // backtrack by bit array matrics
    DFS(Node(nrow_ -1, ncol_ -1), 0);
  }

  void DFS(const Node& cn, const int must_go_dir) {
    // must_go_dir, 0: not required, 1: must go left, 2: must go above
    auto prev = align_path_.empty() ? Node(0,0) : align_path_.back();
    align_path_.push_back(cn);
    if (cn.i == 0 && cn.j == 0) {
      paths_.push_back(align_path_);
    }
    else {
      if (must_go_dir == 1) {
        int next_must_go = hori_whole_[cn.i][cn.j] && hori_left_half_[cn.i][cn.j] ? 1 : 0;
        DFS(Node(cn.i, cn.j-1), next_must_go);
      }
      else if (must_go_dir == 2){
        int next_must_go = vert_whole_[cn.i][cn.j] && vert_top_half_[cn.i][cn.j] ? 2 : 0;
        DFS(Node(cn.i-1, cn.j), next_must_go);
      }
      else {
        if (diag_whole_[cn.i][cn.j]) {
          DFS(Node(cn.i - 1, cn.j - 1), 0);
        }
        if (vert_whole_[cn.i][cn.j]) {
          if (vert_bottom_half_[cn.i][cn.j]) {
            if (cn.i + 1 != prev.i || cn.j != prev.j) return;
          }
          int next_must_go = vert_top_half_[cn.i][cn.j] ? 2 : 0;
          DFS(Node(cn.i - 1, cn.j), next_must_go);
        }
        if (hori_whole_[cn.i][cn.j]) {
          if (hori_right_half_[cn.i][cn.j]) {
            if (cn.i != prev.i || cn.j + 1 != prev.j) return;
          }
          int next_must_go = hori_left_half_[cn.i][cn.j] ? 1 : 0;
          DFS(Node(cn.i, cn.j - 1), next_must_go);
        }
      }
    }
    align_path_.pop_back();
  }

  decltype(auto) Paths() const {
    return (paths_);
  }

  decltype(auto) Path() const {
    return paths_.front();
  }

  void PrintAllPaths() const {
//    int i = 2, j=4;
//    std::cout<< "edge assignment\n";
//    std::cout << vert_whole_[i][j] << "\t" << hori_whole_[i][j] << "\t" << diag_whole_[i][j] << "\t";
//    std::cout << vert_top_half_[i][j] << "\t" << vert_bottom_half_[i][j] << "\t" << hori_left_half_[i][j] << "\t" << hori_right_half_[i][j];
//    std::cout << "\n";
    for (auto& p : paths_) {
      Alignment a(rowstring_, colstring_, p);
      std::cerr << a;
      std::cerr << '\n';
    }
  };

//  void Print() const {
//    std::cerr << "R\n";
//    std::cerr << R_;
//    std::cerr << "P\n";
//    std::cerr << P_;
//    std::cerr << "Q\n";
//    std::cerr << Q_;
//    std::cerr << "a\n";
//    std::cerr << vert_whole_; // a
//    std::cerr << "b\n";
//    std::cerr << hori_whole_; // b
//    std::cerr << "c\n";
//    std::cerr << diag_whole_; // c
//    std::cerr << "d\n";
//    std::cerr << vert_top_half_; // d
//    std::cerr << "e\n";
//    std::cerr << vert_bottom_half_; // e
//    std::cerr << "f\n";
//    std::cerr << hori_left_half_; // f
//    std::cerr << "g\n";
//    std::cerr << hori_right_half_; // g
//  }

};

#endif //ADAPTERTRIM_CPPUTIL_INCLUDE_GOTOH_H_
