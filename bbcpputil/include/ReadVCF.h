//
// Created by Ruolin Liu on 1/2/20.
//

#ifndef REALIGN_INCLUDE_READVCF_H_
#define REALIGN_INCLUDE_READVCF_H_

#include <fstream>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>
#include "Variant.h"

namespace cpputil {

using BiAllele = std::pair<cpputil::Variant, cpputil::Variant>;

class BCFReader {
  bcf_hdr_t *hdr_;
  bcf_srs_t *fp_;
  std::string bcf_file_;

  void _Load(const char *bcf_file) {
    bcf_file_ = bcf_file;
    fp_ = bcf_sr_init();
    bcf_sr_set_opt(fp_, BCF_SR_REQUIRE_IDX);
    if (!bcf_sr_add_reader(fp_, bcf_file)) {
      std::cerr << "Failed to read from " << bcf_file_ <<std::endl;
      std::cerr << "File indexed? "<<std::endl;
      exit(0);
    }
    hdr_ = fp_->readers[0].header;
  }

//  int Gettid_(const string &contig) const {
//    int32_t tid = bcf_hdr_name2id(hdr_, contig.c_str());
//    if (tid == -1) {
//      if (contig.substr(0, 3) == "chr") {
//        tid = bcf_hdr_name2id(hdr_, (contig.substr(3)).c_str());
//      } else {
//        tid = bcf_hdr_name2id(hdr_, ("chr" + contig).c_str());
//      }
//    }
//    return tid;
//  }

 public:
  BCFReader() : hdr_(NULL), fp_(NULL) {}

//  BCFReader(const char* bcf_file) : bcf_file_(bcf_file){
//    fp_ = hts_open(bcf_file, "rb");
//    hdr_ = bcf_hdr_read(fp_);
//    idx_ = bcf_index_load(bcf_file);
//    assert (idx_ != NULL && fp_ != NULL && hdr_ != NULL);
//    enum htsExactFormat format = hts_get_format(fp_)->format;
//    assert(format == bcf);
//    _Load(bcf_file);
//  }

  void Open(const char *bcf_file, const std::string sample = "") {
    _Load(bcf_file);
    if (!sample.empty()) {
      int ret = SubsetSample(sample.c_str());
      if (ret != 0) {
        throw std::runtime_error("sample name " + sample + " not found in BCF file\n");
      }
    } else {
      if (bcf_hdr_nsamples(hdr_) != 0 && bcf_hdr_nsamples(hdr_) != 1) {
        throw std::runtime_error("Must give a sample name for muti-sample BCF file\n");
      }
    }
  }

  bool IsOpen() const {
    return fp_ != NULL;
  }

  decltype(auto) header() const {
    return (hdr_);
  }

  int SubsetSample(const char *sample) {
    int ret = bcf_hdr_set_samples(hdr_, sample, 0);
    return ret;
  }

  ~BCFReader() {
    if (fp_) bcf_sr_destroy(fp_);
  }

  std::string VariantStr(bcf1_t *rec) {
    int32_t *gt = NULL, ngt_arr = 0;
    int nsmpl = bcf_hdr_nsamples(hdr_);
    int ngt = bcf_get_format_int32(hdr_, rec, "GT", &gt, &ngt_arr);
    assert(ngt > 0);
    std::string ret;
    ret = std::string(bcf_hdr_id2name(hdr_, rec->rid)) + ":";
    ret += std::to_string(rec->pos) + "\t";
    for (int i = 0; i < nsmpl; ++i) {
      for (int j = 0; j < 2; ++j) {
        int allele_index = bcf_gt_allele(gt[i * 2 + j]);
        //std::cerr << allele_index << std::endl;
        //std::cerr << rec->d.allele[allele_index] << std::endl;
        ret += std::string(rec->d.allele[allele_index]) + "/";
      }
      ret += "\t";
    }
    return ret;
  }

//  std::map<int32_t, BiAllele> get_variants(const string &contig, const int32_t &start, const int32_t &end,
//                                           bool het_only = true, bool phased_only = true) const {
//    /*
//     * Get small variants over a interval
//     * start, end are 0-base half open interval
//     * maternal allele has lower index than paternal allele
//     */
//    //Must has sample column
//    if (bcf_hdr_nsamples(hdr_) == 0) {
//      std::cerr << "Error: vcf does not have sample column. Thus genoytpe unclear\n";
//      exit(0);
//    }
//    std::map<int32_t, BiAllele> res;
//    int tid = Gettid_(contig);
//    if (tid == -1) {
//      return res;
//    }
//    hts_itr_t *itr = bcf_itr_queryi(idx_, tid, start, end);
//    while (bcf_itr_next(fp_, itr, rec) >= 0) {
//      bcf_unpack(rec, BCF_UN_ALL);
//      bcf_subset_format(hdr_, rec);
//
//      int32_t nsvt = 0;
//      char *svt = NULL;
//      int ret = bcf_get_info_string(hdr_, rec, "SVTYPE", &svt, &nsvt);
//      free(svt);
//      if (ret >= 0) continue;
//      int32_t *gt = NULL, ngt_arr = 0;
//      int ngt = bcf_get_format_int32(hdr_, rec, "GT", &gt, &ngt_arr);
//      if (ngt != 2) {
//        std::cerr << "Error: variant is not diploid " << rec->rid << "\t" << rec->pos + 1 << std::endl;
//        bcf_itr_destroy(itr);
//        exit(0);
//      }
//      int ai0 = bcf_gt_allele(gt[0]);
//      int ai1 = bcf_gt_allele(gt[1]);
//      if (ai0 >= rec->n_allele or ai0 < 0 or ai1 >= rec->n_allele or ai1 < 0) {
//        std::cerr << "Error: allele index not correct " << rec->rid << "\t" << rec->pos + 1 << std::endl;
//        bcf_itr_destroy(itr);
//        exit(0);
//      }
//      if (phased_only and not bcf_gt_is_phased(ai0) and not bcf_gt_is_phased(ai1)) {
//        continue;
//      }
//      std::string mat_allele(rec->d.allele[ai0]);
//      std::string pat_allele(rec->d.allele[ai1]);
//      if (het_only and mat_allele == pat_allele) {
//        continue;
//      }
//      auto v0 = cpputil::Variant(contig, rec->pos, mat_allele);
//      auto v1 = cpputil::Variant(contig, rec->pos, pat_allele);
//      res[rec->pos] = std::make_pair(v0, v1);
//      free(gt);
//    } //end while
//    bcf_itr_destroy(itr);
//    return res;
//  }

  double var_population_freq(const string contig, const int64_t pos, const string ref, const string alt, const double lowb=1e-6) const {
    //std::cerr << contig << "\t" << pos << "\t" << ref << "\t" << alt << "\n";
//    int tid = Gettid_(contig);
//    if (tid == -1) {
//      return lowb;
//    }
    string region_str = contig + ":" + std::to_string(pos+1); //bcf_sr take 1 base inclusive coordinates
    fp_->regions = bcf_sr_regions_init(region_str.c_str(), 0 ,0,1,-2);
    if (!fp_->regions) {
      std::cerr << "Failed to read the regions: " << region_str << std::endl;
      exit(0);
    }
    fp_->explicit_regs = 1;

    while (bcf_sr_next_line(fp_)) {
      bcf1_t* rec = fp_->readers[0].buffer[0];
      //bcf_subset_format(hdr_, rec); // no sample column usually
      if (pos != rec->pos or std::string(rec->d.allele[0]) != ref) {
        //std::cerr << "rec->pos: " << rec->pos << " \t" << "rec->d.allele[0] " << rec->d.allele[0] << std::endl;
        continue;
      }
      for (int ii = 1; ii < rec->n_allele; ++ii) {
        //std::cerr << "rec->pos: " << rec->pos << " \t" << "rec->d.allele[0] " << rec->d.allele[0] << " \t" << " alt allele " << rec->d.allele[ii] << std::endl;
        if (std::string(rec->d.allele[ii]) == alt) {
          int32_t naf = 0;
          float* afp = NULL;
          int ret = bcf_get_info_float(hdr_, rec, "AF", &afp, &naf);
          float af = *afp;
          free(afp);
          return af;
        }
      }
    }
    return lowb;
  }


  bool var_exist(const string &contig, const int32_t &pos, string allele = "") const {
    /*
     * If allele is empty, assume a SNP mask and return true if just position match.
     */
    string region_str = contig + ":" + std::to_string(pos+1); //bcf_sr take 1 base inclusive coordinates
    if (fp_->regions) bcf_sr_regions_destroy(fp_->regions);
    fp_->regions = bcf_sr_regions_init(region_str.c_str(), 0 ,0,1,-2);
    if (!fp_->regions) {
      std::cerr << "Failed to read the regions: " << region_str << std::endl;
      exit(0);
    }
    while (bcf_sr_next_line(fp_)) {
      bcf1_t* rec = fp_->readers[0].buffer[0];
      int32_t nsvt = 0;
      char *svt = NULL;
      int ret = bcf_get_info_string(hdr_, rec, "SVTYPE", &svt, &nsvt);
      free(svt);
      if (ret >= 0) continue;
      rec->d.var_type = -1;
      if (bcf_hdr_nsamples(hdr_) == 0) {
        //no sample column, if any of the alt match
        for (int ii = 1; ii < rec->n_allele; ++ii) {
          int t = bcf_get_variant_type(rec, ii);
          if (allele.empty() && pos == rec->pos) return true;
          std::string vcfallele(rec->d.allele[ii]);
          if (pos == rec->pos && vcfallele == allele) {
            return true;
          }
        }
      } else {
        //otherwise select the first sample. If multisample is given, should subset by sample first name
        int32_t *gt = NULL;
        int ngt_arr = 0;
        int ngt = bcf_get_format_int32(hdr_, rec, "GT", &gt, &ngt_arr);
        if (ngt < 0) { // no genotype
          for (int ii = 1; ii < rec->n_allele; ++ii) {
            int t = bcf_get_variant_type(rec, ii);
            if (allele.empty() &&  pos == rec->pos) {
              free(gt);
              return true;
            }
            std::string vcfallele(rec->d.allele[ii]);
            if (pos == rec->pos && vcfallele == allele) {
              free(gt);
              return true;
            }
          }
        } else {
          assert(ngt == 2); // assume single diploid sample
          for (int j = 0; j < 2; ++j) { // 2 alleles, maternal, paternal
            int allele_index = bcf_gt_allele(gt[j]);
            if (allele_index == 0) continue;
            if (allele_index < 0 or allele_index >= rec->n_allele) {
              std::cerr << "Error: allele index not exist " << rec->rid << "\t" << rec->pos
                        << "\n";
              exit(0);
            }
            int t = bcf_get_variant_type(rec, allele_index);
            if (allele.empty() && (t == VCF_MNP || t == VCF_SNP) && pos == rec->pos) {
              free(gt);
              return true;
            }
            std::string vcfallele(rec->d.allele[allele_index]);
            if (pos == rec->pos && vcfallele== allele) {
              free(gt);
              return true;
            }
          }
        }
        free(gt);
      } //end else
    } //end while
    return false;
  }

  void vcf_to_blacklist_snv(const SeqLib::GenomicRegion& gr, const SeqLib::BamHeader bh, std::set<int32_t>& res){
    if (fp_->regions) bcf_sr_regions_destroy(fp_->regions);
    std::string region_str = gr.ChrName(bh) + ":" + std::to_string(gr.pos1+1) + "-" + std::to_string(gr.pos2); //bcftools take 1-base inclusive interval
    fp_->regions = bcf_sr_regions_init(region_str.c_str(), 0 ,0,1,-2);
    if (!fp_->regions) {
      std::cerr << "Failed to read the regions: " << region_str << std::endl;
      exit(0);
    }
    while (bcf_sr_next_line(fp_)) {
      bcf1_t* rec = fp_->readers[0].buffer[0];
      rec->d.var_type = -1;
      for (int i=1; i<rec->n_allele; i++)
      {
        int t = bcf_get_variant_type(rec, i);
        if (t == VCF_MNP or t == VCF_SNP) {
          for (unsigned j = 0; j< strlen(rec->d.allele[i]); ++j) {
            res.insert(rec->pos + j);
          }
        }
      }
    }
  }

};

bool all_true_mut(vector<bool> v) {
  return count(v.begin(), v.end(), true) == (int) v.size();
}

template<class VariantsReader>
std::vector<bool> search_var_in_database(const VariantsReader& variant_reader, const cpputil::Variant& var, std::ofstream& outf,
                            std::string anno, bool fall_back_pos_only, int min_baseq, bool out_all_snps) {
  std::vector<bool> real_muts(var.alt_seq.size());
  if (variant_reader.var_exist(var.contig, var.contig_start, var.alt_seq)) { // known var
    if (var.isIndel()) {
      if (outf.is_open()) outf << var << '\t' << anno << '\n';
      return {true};
    }
    else {
      std::fill(real_muts.begin(), real_muts.end(), true);
      if (outf.is_open() && (out_all_snps || var.var_qual >= min_baseq)) outf << var << '\t' << anno << '\n';
    }
  } else {
    if (var.isMNV()) { // rescue if MNV
      auto avars = cpputil::var_atomize(var);
      for (unsigned i =0 ; i < avars.size(); ++i) {
        if (variant_reader.var_exist(avars[i].contig, avars[i].contig_start, avars[i].alt_seq)) {
          real_muts[i] = true;
          if (outf.is_open() and (out_all_snps || avars[i].var_qual >= min_baseq)) outf << avars[i] << '\t' << anno << '\n';
        } else if (fall_back_pos_only && variant_reader.var_exist(avars[i].contig, avars[i].contig_start)) {
          real_muts[i] = true;
        }
      }
    } else {
      if (fall_back_pos_only && variant_reader.var_exist(var.contig, var.contig_start)) {
        return {true};
      }
    }
  }
  return real_muts;
}


}

#endif //REALIGN_INCLUDE_READVCF_H_
