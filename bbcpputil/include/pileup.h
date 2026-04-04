//
// Created by Ruolin Liu on 10/28/21.
//

#ifndef CODECSUITE_BBCPPUTIL_INCLUDE_PILEUP_H_
#define CODECSUITE_BBCPPUTIL_INCLUDE_PILEUP_H_

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <ctype.h>
#include <unistd.h>
#include <string>
#include <memory>
#include <map>

#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "SeqLib/BamWalker.h"
#include "SeqLib/RefGenome.h"
#include "BamRecordExt.h"

//adapted from samtools
//https://github.com/samtools/htslib/blob/9672589346459d675d62851d5b7b5f2e5c919076/test/pileup.c

namespace cpputil {

class PileHandler {
public:
  std::string fname;
  std::shared_ptr<htsFile> fp;
  std::shared_ptr<bam_hdr_t> fp_hdr;
  std::shared_ptr<hts_idx_t> idx;
  std::shared_ptr<hts_itr_t> iter;
  int min_mapQ;

  PileHandler(): min_mapQ(0) {}

  PileHandler(std::string f, int mq): fname(f), min_mapQ(mq) {
    fp = std::shared_ptr<htsFile>(hts_open(f.c_str(), "r"), htsFile_delete());
    if (!fp) {
      fprintf(stderr, "Couldn't open \"%s\" : %s", fname.c_str(), strerror(errno));
      exit(1);
    }

    idx = std::shared_ptr<hts_idx_t>(sam_index_load(fp.get(), fname.c_str()), idx_delete());
    if (idx == 0) {
      fprintf(stderr, "\"%s\" is not indexed : %s", fname.c_str(), strerror(errno));
      exit(1);
    }

    fp_hdr = std::shared_ptr<bam_hdr_t> (sam_hdr_read(fp.get()), bam_hdr_delete());
    if (!fp_hdr.get()) {
      fprintf(stderr, "Couldn't read header from \"%s\" : %s",
              fname.c_str(), strerror(errno));
      exit(1);
    }
  }
//  ~PileHandler() {
//    if (fp_hdr) bam_hdr_destroy(fp_hdr);
//    if (fp) sam_close(fp);
//    if (idx) hts_idx_destroy(idx);
//    if (iter) hts_itr_destroy(iter);
//  }
  static int readaln(void *data, bam1_t *b) {
    PileHandler *g = (PileHandler*)data;
    int ret;

    while (1) {
      ret = g->iter ? sam_itr_next(g->fp.get(), g->iter.get(), b) : sam_read1(g->fp.get(), g->fp_hdr.get(), b);
      if (ret < 0) break;
      if ( b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP) ) continue;
      if ((int)b->core.qual < g->min_mapQ) continue;
      break;
    }
    return ret;
  }
};

std::vector<BamPileup> GetPileupReads(PileHandler *input,
                                      std::string chrom,
                                      int tpos,
                                      bool  handle_overlap){
  bam_mplp_t mplp = NULL;
  const bam_pileup1_t *pileups[1] = { NULL };
  int n_plp[1] = { 0 };
  int tid, pos, n = 0;
  tid = bam_name2id(input->fp_hdr.get(), chrom.c_str());
  if (tid < 0) {
    fprintf(stderr, "%s not found for \"%s\"\n", chrom.c_str(), input->fname.c_str());
    exit(1);
  }

  if (input->iter) input->iter.reset();
  input->iter = std::shared_ptr<hts_itr_t>( sam_itr_queryi(input->idx.get(), tid, tpos, tpos + 1), hts_itr_delete());

  mplp = bam_mplp_init(1, PileHandler::readaln, (void **) &input);
  if (!mplp) {
    perror("bam_plp_init");
    exit(1);
  }
  if (handle_overlap) {
    bam_mplp_init_overlaps(mplp);
  }

  int non_del_cnt = 0, d = 0;
  std::vector<BamPileup> res;
  while ((n = bam_mplp_auto(mplp, &tid, &pos, n_plp, pileups)) > 0) {
    if (tid < 0) break;
    if (tid >= input->fp_hdr->n_targets) {
      fprintf(stderr,
              "bam_mplp_auto returned tid %d >= header n_targets %d\n",
              tid, input->fp_hdr->n_targets);
      exit(1);
    }
    if (pos < tpos || pos > tpos) continue;
    char cc = 0, qq = 0;
    for (int j = 0; j  < n_plp[0]; ++j) {
      const bam_pileup1_t *pi = pileups[0] + j;
      res.emplace_back(pi);
    }
    //printf("%s\t%d\t%d\t%d\t%d\n", input->fp_hdr->target_name[tid], pos+1, n_plp[0], nuc_cnt, non_del_cnt);
  }
  bam_mplp_destroy(mplp);
  if (n < 0) {
    fprintf(stderr, "bam_plp_auto failed for \"%s\"\n", input->fname.c_str());
    exit(1);
  }
  return res;
}

static int ScanAllele(PileHandler *input, std::string chrom,
                      int tpos,
                      char nuc,
                      bool  handle_overlap,
                      int &nuc_cnt,
                      int minbq = 10) {
  /*
   * find num of supports for an allele
   * Collapse read-pairs to avoid counting overlaps
   */
  nuc_cnt = 0;
  bam_mplp_t mplp = NULL;
  const bam_pileup1_t *pileups[1] = { NULL };
  int n_plp[1] = { 0 };
  int tid, pos, n = 0;
  tid = bam_name2id(input->fp_hdr.get(), chrom.c_str());
  if (tid < 0) {
    fprintf(stderr, "%s not found for \"%s\"\n", chrom.c_str(), input->fname.c_str());
    return -1;
  }

  if (input->iter) input->iter.reset();
  input->iter = std::shared_ptr<hts_itr_t>( sam_itr_queryi(input->idx.get(), tid, tpos, tpos + 1), hts_itr_delete());

  mplp = bam_mplp_init(1, PileHandler::readaln, (void **) &input);
  if (!mplp) {
    perror("bam_plp_init");
    return -1;
  }
  if (handle_overlap) {
    bam_mplp_init_overlaps(mplp);
  }

  int non_del_cnt = 0, d = 0;
  while ((n = bam_mplp_auto(mplp, &tid, &pos, n_plp, pileups)) > 0) {
    if (tid < 0) break;
    if (tid >= input->fp_hdr->n_targets) {
      fprintf(stderr,
              "bam_mplp_auto returned tid %d >= header n_targets %d\n",
              tid, input->fp_hdr->n_targets);
      return -1;
    }
    if (pos < tpos || pos > tpos) continue;
    char cc = 0;
    for (int j = 0; j  < n_plp[0]; ++j) {
      const bam_pileup1_t *pi = pileups[0] + j;
      if (pi->is_del or pi->is_refskip) ++d;
      else {
        if (pi->qpos < pi->b->core.l_qseq) {
          uint8_t *qq = bam_get_qual(pi->b);
          if (qq[pi->qpos] >= minbq) {
            ++non_del_cnt;
            cc = seq_nt16_str[bam_seqi(bam_get_seq(pi->b), pi->qpos)];
            if (cc == nuc) ++nuc_cnt;
          }
        }
      }
    }
    //printf("%s\t%d\t%d\t%d\t%d\n", input->fp_hdr->target_name[tid], pos+1, n_plp[0], nuc_cnt, non_del_cnt);
  }
  bam_mplp_destroy(mplp);
  if (n < 0) {
    fprintf(stderr, "bam_plp_auto failed for \"%s\"\n", input->fname.c_str());
    return -1;
  }
  return non_del_cnt;
}

static int ScanIndel(PileHandler *input,
                     bool diff_bam,
                     const SeqLib::RefGenome& ref,
                     const std::string& chrom,
                      int tpos,
                      int indel,
                      const std::string& seq,
                      bool  handle_overlap,
                      int wiggle,
                      float nearby_germ_min_vaf,
                      int mindepth,
                      int minbq,
                      int &exact_match,
                      int &fuzzy_match) {
  /*
   * indel <0 for DEL >0 for INS
   * inseq is the insert seq, empty for DEL or not checking ins
   */

  //std::cerr << "scan indels" <<std::endl;
  bam_mplp_t mplp = NULL;
  const bam_pileup1_t *pileups[1] = { NULL };
  int n_plp[1] = { 0 };
  int tid = bam_name2id(input->fp_hdr.get(), chrom.c_str());
  if (tid < 0) {
    fprintf(stderr, "%s not found for \"%s\"\n", chrom.c_str(), input->fname.c_str());
    bam_mplp_destroy(mplp);
    return -1;
  }
  int slen = indel < 0 ? abs(indel) : 1;
  int depth = 0;
  exact_match = 0;
  int search;
  //search from the middle and clockwise
  for (int search1 = 0; search1 < slen + 2 * wiggle; ++search1) {
    if (search1 >= slen + wiggle) search = slen + wiggle - search1 - 1;
    else search = search1;

    if (input->iter) input->iter.reset();
    input->iter = std::shared_ptr<hts_itr_t>( sam_itr_queryi(input->idx.get(), tid, tpos + search, tpos + search + 1), hts_itr_delete());

    mplp = bam_mplp_init(1, PileHandler::readaln, (void **) &input);
    if (!mplp) {
      perror("bam_plp_init");
      bam_mplp_destroy(mplp);
      return -1;
    }
    if (handle_overlap) {
      bam_mplp_init_overlaps(mplp);
    }

    fuzzy_match = 0;
    int pos, n=0, d=0;
    std::set<std::string> readnames;
    //std::cerr << "search: " << search << "\t" << "tpos: " << tpos << std::endl;
    std::map<char, int> nuc_cnt;
    std::string refbase = ref.QueryRegion(chrom, tpos+search, tpos+search);
    while ((n = bam_mplp_auto(mplp, &tid, &pos, n_plp, pileups)) > 0) {
      if (tid < 0) break;
      if (tid >= input->fp_hdr->n_targets) {
        fprintf(stderr,
                "bam_mplp_auto returned tid %d >= header n_targets %d\n",
                tid, input->fp_hdr->n_targets);
        return -1;
      }
      if (pos != tpos + search) {continue;}

      for (int j = 0; j  < n_plp[0]; ++j) {
        const bam_pileup1_t *pi = pileups[0] + j;
        uint8_t *qq = bam_get_qual(pi->b);
//        if (diff_bam) {
//          std::cerr <<(int) qq[pi->qpos] << ", " << indel << ", " << search <<", " <<pi->indel <<std::endl;
//        }
        if (qq[pi->qpos] < minbq or pi->is_refskip)  continue;
        ++d;
        if (not handle_overlap) {
          readnames.emplace(bam_get_qname(pi->b));
        }
        if (pi->indel != 0) { // is anchor base of a indel
          if (pi->qpos < pi->b->core.l_qseq) {
            if (pi->indel < 0) { // DEL
              if (indel == pi->indel && search == 0) ++exact_match;
              else {
                ++fuzzy_match;
              }
            } else { // INS
              std::string inseq(pi->indel, '.');
              for (int32_t i = 0; i < pi->indel; ++i) {
                if (pi->qpos + i + 1 == pi->b->core.l_qseq) break;
                inseq[i] = seq_nt16_str[bam_seqi(bam_get_seq(pi->b), pi->qpos + i + 1)];
              }
              if (indel == pi->indel && inseq == seq) ++exact_match;
              else ++fuzzy_match;
            }
          }
        }
        if (!pi->is_del) { //if not DEL
          char cc = seq_nt16_str[bam_seqi(bam_get_seq(pi->b), pi->qpos)];
          ++nuc_cnt[cc];
        }
        else {
          ++nuc_cnt['-'];
        }
      }
      //printf("%s\t%d\t%d\t%d\t%d\n", input->fp_hdr->target_name[tid], pos+1, n_plp[0], nuc_cnt, non_del_cnt);
    }
    if (not handle_overlap) {
      d = (int) readnames.size();
    }
    if (search == 0) depth = d;
    if (not handle_overlap) {
      depth = (int) readnames.size();
    }
    if (n < 0) {
      fprintf(stderr, "bam_plp_auto failed for \"%s\"\n", input->fname.c_str());
      return -1;
    }
    // Exact match
    if (diff_bam && exact_match > 0) {
      bam_mplp_destroy(mplp);
      return depth;
    }
    if ((fuzzy_match > 1 and d <= 150) or (fuzzy_match > 2 and d > 150)) {
      //std::cerr << "found " << fuzzy_match << "\t" << tpos << "," << search << std::endl;
      bam_mplp_destroy(mplp);
      return depth;
    }
    //nearby SNV
    int max=0, secmax=0;
    char maxbase = 0, secmaxbase = 0;
    for (auto const& it: nuc_cnt) {
      if (it.second > max) {
        secmax = max;
        max = it.second;
        secmaxbase = maxbase;
        maxbase = it.first;
      } else if (it.second > secmax) {
        secmax = it.second;
        secmaxbase = it.first;
      }
    }
    if ((depth >= mindepth && max > 0 && maxbase != refbase[0] && max > d * nearby_germ_min_vaf) ||
      (depth >= mindepth && secmax > 0 && secmaxbase != refbase[0] && secmax > d * nearby_germ_min_vaf)) {
    //  std::cerr << search << std::endl;
      fuzzy_match = 2;
      bam_mplp_destroy(mplp);
      return depth;
    }
    bam_mplp_destroy(mplp);
  }
  //std::cerr << "depth " << non_del_cnt << std::endl;
  return depth;
}

int GenotypeVariant(PileHandler *inputbam,
                           const SeqLib::RefGenome& ref,
                           const std::string& chrom,
                           int tpos,
                           int indel,
                           const std::string& seq,
                           bool  handle_overlap,
                           int wiggle,
                           float germ_min_vaf,
                           int mindepth,
                           int minbq) {
  /* INPUT
   * 0 SNP, >0 INS, <0 DEL
   * seq: Inserted seq or SNV base
   * RETURN
   * -2: Error
   * 0 : germline
   * >0 : depth
   * other: error
   */
  //std::cerr << "scan indels" <<std::endl;
  int depth = 0;
  if (indel == 0) { // SNV
    int alt_count = 0;
    depth = ScanAllele(inputbam, chrom, tpos, seq[0],handle_overlap, alt_count, minbq);
    if (depth == -1)
      return -2;
    if (depth >= mindepth && (float) alt_count >= germ_min_vaf * depth)
      return 0;
  } else {
    int exact_match = 0, fuzzy_match = 0;
    depth = ScanIndel(inputbam, false, ref, chrom, tpos, indel, seq, handle_overlap, wiggle,
                      0.25, mindepth, minbq, exact_match, fuzzy_match);
    if (depth == -1)
      return -2;
    if (depth >= mindepth && (fuzzy_match > 1 || (float) exact_match >= germ_min_vaf * depth))
      return 0;
#if 0
    // search nearby germline SNP when calling INDELs
    bam_mplp_t mplp = NULL;
    const bam_pileup1_t *pileups[1] = { NULL };
    int n_plp[1] = { 0 };
    int tid = bam_name2id(germbam->fp_hdr.get(), chrom.c_str());
    if (tid < 0) {
      fprintf(stderr, "%s not found for \"%s\"\n", chrom.c_str(), germbam->fname.c_str());
      bam_mplp_destroy(mplp);
      return -2;
    }
    int slen = indel < 0 ? abs(indel) : 1;
    //std::cerr << "slen: " << slen << std::endl;
    for (int search = -wiggle; search < slen + wiggle; ++search) {
      if (germbam->iter) germbam->iter.reset();
      germbam->iter = std::shared_ptr<hts_itr_t>( sam_itr_queryi(germbam->idx.get(), tid, tpos + search, tpos + search + 1), hts_itr_delete());

      mplp = bam_mplp_init(1, PileHandler::readaln, (void **) &germbam);
      if (!mplp) {
        perror("bam_plp_init");
        bam_mplp_destroy(mplp);
        return -2;
      }
      if (handle_overlap) {
        bam_mplp_init_overlaps(mplp);
      }

      int non_skip_cnt = 0;
      int pos, n=0, exact_match=0;
      std::map<int, int> nuc_cnt;
      //std::cerr << "search: " << search << "\t" << "tpos: " << tpos << std::endl;
      std::string refbase = ref.QueryRegion(chrom, tpos+search, tpos+search);
      while ((n = bam_mplp_auto(mplp, &tid, &pos, n_plp, pileups)) > 0) {
        if (tid < 0) break;
        if (tid >= germbam->fp_hdr->n_targets) {
          fprintf(stderr,
                  "bam_mplp_auto returned tid %d >= header n_targets %d\n",
                  tid, germbam->fp_hdr->n_targets);
          return -2;
        }
        if (pos != tpos + search) {continue;}
//        std::cerr << "pos: " << pos << std::endl;
        for (int j = 0; j  < n_plp[0]; ++j) {
          const bam_pileup1_t *pi = pileups[0] + j;
          if (pi->is_refskip) continue;
          else {
            if (pi->qpos < pi->b->core.l_qseq) {
              uint8_t *qq = bam_get_qual(pi->b);
              if (qq[pi->qpos] < 10) {
                continue;
              }
              ++non_skip_cnt;
              if (!pi->is_del) {
                char cc = seq_nt16_str[bam_seqi(bam_get_seq(pi->b), pi->qpos)];
                if (pi->indel != 0) {
                  ++nuc_cnt[pi->indel];
                  if (cc != refbase[0]) ++nuc_cnt[(int)cc + 128];
                } else {
                  ++nuc_cnt[(int)cc + 128]; // avoid overlap between nuc ascii and indels
                }
              }
            }
          }
        }
      }

      if (n < 0) {
        fprintf(stderr, "bam_plp_auto failed for \"%s\"\n", germbam->fname.c_str());
        bam_mplp_destroy(mplp);
        return -2;
      }
      int max=0, secmax=0;
      int maxbase = 0, secmaxbase = 0;
      for (auto const& it: nuc_cnt) {
        if (it.second > max) {
          secmax = max;
          max = it.second;
          secmaxbase = maxbase;
          maxbase = it.first;
        } else if (it.second > secmax) {
          secmax = it.second;
          secmaxbase = it.first;
        }
      }
//      std::cerr << "max: " << maxbase <<", " << max << std::endl;
//      std::cerr << "sec max: " << secmaxbase <<", " << secmax << std::endl;
//      std::cerr << "non_skip_count " << non_skip_cnt << std::endl;
      int refint = (int) refbase[0] + 128;
      if ((non_skip_cnt >= mindepth && max > 0 && maxbase != refint && max > non_skip_cnt * germ_min_vaf) ||
          (non_skip_cnt >= mindepth && secmax > 0 && secmaxbase != refint && secmax > non_skip_cnt * germ_min_vaf)) {
        bam_mplp_destroy(mplp);
        return -1;
      }
      bam_mplp_destroy(mplp);
    }
#endif
  }
  return depth;
}


}

#endif //CODECSUITE_BBCPPUTIL_INCLUDE_PILEUP_H_
