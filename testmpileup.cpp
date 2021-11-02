//
// Created by Ruolin Liu on 10/28/21.
//
#include "pileup.h"

int main(int argc, char **argv) {
  cpputil::ptest_t g = { NULL, NULL, NULL };
  int use_mpileup = 0, opt;

  while ((opt = getopt(argc, argv, "m")) != -1) {
    switch (opt) {
      case 'm':
        use_mpileup = 1;
        break;
        default:
          fprintf(stderr, "Usage: %s [-m] <sorted.sam>\n", argv[0]);
          return EXIT_FAILURE;
    }
  }

  if (optind >= argc) {
    fprintf(stderr, "Usage: %s [-m] <sorted.sam>\n", argv[0]);
    return EXIT_FAILURE;
  }
  hts_idx_t *idx;

  g.fname = argv[optind];
  g.fp = sam_open(g.fname, "r");
  idx = sam_index_load(g.fp, g.fname);
  if (!g.fp) {
    fprintf(stderr, "Couldn't open \"%s\" : %s", g.fname, strerror(errno));
    goto fail;
  }
  if (idx == 0) {
    fprintf(stderr, "%s is not indexed", g.fname, strerror(errno));
  }
  g.fp_hdr = sam_hdr_read(g.fp);
  if (!g.fp_hdr) {
    fprintf(stderr, "Couldn't read header from \"%s\" : %s",
            g.fname, strerror(errno));
    goto fail;
  }
  int tid, beg, end, pos;
  tid = bam_name2id(g.fp_hdr, "3");
  if (tid < 0)
    goto fail;
  beg = 10345306;
  end = 10345307;

  g.min_mapQ = 0;
  if (g.iter) hts_itr_destroy(g.iter);
  g.iter = sam_itr_queryi(idx, tid, beg, end);

  if (test_mpileup(&g) < 0)
    goto fail;

  bam_hdr_destroy(g.fp_hdr);
  sam_close(g.fp);
  return EXIT_SUCCESS;

  fail:
  fprintf(stderr, "failed\n");
  if (g.fp_hdr) bam_hdr_destroy(g.fp_hdr);
  if (g.fp) sam_close(g.fp);
  return EXIT_FAILURE;
}

