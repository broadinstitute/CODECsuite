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

#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "SeqLib/BamWalker.h"
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
      fprintf(stderr, "%s is not indexed", fname.c_str(), strerror(errno));
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
    char cc = 0, qq = 0;
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
                      const std::string& chrom,
                      int tpos,
                      int indel,
                      const std::string& seq,
                      bool  handle_overlap,
                      int &cnt_found,
                      int minbq = 10) {
  /*
   * indel <0 for DEL >0 for INS
   * inseq is the insert seq, empty for DEL or not checking ins
   */

  cnt_found = 0;
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
    char qq = 0;
    for (int j = 0; j  < n_plp[0]; ++j) {
      const bam_pileup1_t *pi = pileups[0] + j;
      if (pi->indel != 0) {
        if (pi->qpos < pi->b->core.l_qseq) {
          if (pi->indel < 0) {
            if (indel == pi->indel) ++cnt_found;
          } else {
            std::string inseq(pi->indel + 1, '.');
            for (int32_t i = 0; i < pi->indel + 1; ++i) {
              if (pi->qpos + i == pi->b->core.l_qseq) break;
              inseq[i] = seq_nt16_str[bam_seqi(bam_get_seq(pi->b), pi->qpos + i)];
            }
            if (indel == pi->indel && inseq == seq) ++cnt_found;
          }
        }
      }
      if (pi->is_del or pi->is_refskip) ++d;
      else {
        if (pi->qpos < pi->b->core.l_qseq) {
          uint8_t *qq = bam_get_qual(pi->b);
          if (qq[pi->qpos] >= minbq) {
            ++non_del_cnt;
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

}

#endif //CODECSUITE_BBCPPUTIL_INCLUDE_PILEUP_H_
