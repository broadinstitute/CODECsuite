library(data.table)
library(tidyverse)
library(binom)
library(scales)
library(ggplot2)
library(deconstructSigs)
library(stringi)
library(GenomicRanges)
library(ggpubr)
library(sigfit)
data("cosmic_signatures_v3")
data("cosmic_signatures_v3.2")
data("cosmic_signatures_v2")
library(deconstructSigs)
sessionInfo()
library(insect)
library("cowplot")
library(BSgenome)

to_miredas_errordf <- function(df) {
  df = df %>% select(sample, lane, group, q30_n_bases_eval, q30_n_errors) %>% mutate(sample_id = paste(lane, sample, sep="."))
  df = df %>% pivot_longer(c(q30_n_bases_eval, q30_n_errors))
  df = df %>% separate(name, sep = "_n_", into = c("qcut", "type")) %>% pivot_wider(id_cols = c(sample, lane, group, sample_id, qcut), names_from = type, values_from = value)
  df = df %>% mutate(group = paste(group, qcut, sep="_")) %>% select(-qcut)
  colnames(df)[5] = "total_bases_eval"
  colnames(df)[6] = "total_base_errors"
  df$total_bases_eval = df$total_bases_eval/2
  df$total_base_errors = df$total_base_errors/2
  df
}

prepare_error <- function(folder, group, pattern = "mutation_metrics.txt$", sep='\\.') {
  dfall = data.frame()
  i = 0
  for (f in list.files(folder, pattern=pattern)) {
    print(f)
    df = fread(file.path(folder, f))
    fields = unlist(str_split(f, sep))
    df$sample = fields[2]
    df$lane = fields[1]
    dfall = rbind(dfall, df)
    i = i + 1
  
    }
  dfall$group = group 
  dfall
}

parse_monomer_from_c <- function(error_met_df, mut_fam_df) {
  if ("qcutoff" %in% colnames(error_met_df) & str_detect(error_met_df$qcutoff, "/2")) {
    perbase = canno_base(break_mnv(mut_fam_df)) %>% group_by(ref, alt) %>% summarise(n=sum(read_count))
  } else{
    perbase = canno_base(break_mnv(mut_fam_df)) %>% group_by(ref, alt) %>% summarise(n=n())
  }
  perbase$den = 0
  perbase[perbase$ref == 'C',"den"] = as.numeric(sum(error_met_df$n_C_eval)) + as.numeric(sum(error_met_df$n_G_eval))
  perbase[perbase$ref == 'T',"den"] = as.numeric(sum(error_met_df$n_A_eval)) + as.numeric(sum(error_met_df$n_T_eval))
  perbase$er = perbase$n / perbase$den
  perbase %>% ungroup()
}

monomer_fromc <- function(folder, group, pattern="mutation_metrics.txt$", sep='\\.', idx=1, mutect2_folder = NULL, group_by_sample = FALSE) {
  dfall = data.frame()
  for (f in list.files(folder, pattern=pattern)) {
    print(f)
    met = f
    fam = str_replace(f, "mutation_metrics", "variants_called")
    metdf = fread(file.path(folder, met), nrows=1)
    famdf = fread(file.path(folder, fam))
    fields = unlist(str_split(f, sep))
    if(!is.null(mutect2_folder)) {
      maf = read_maf(file.path(mutect2_folder, paste0(fields[idx], "-filtered.annotated.maf")))
      validated_id =  maf %>% mutate(id = paste(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep="_")) %>% pull(id)
      famdf = famdf %>% mutate(id = paste(chrom,ref_pos ,ref, alt, sep="_"))
      before = nrow(famdf)
      famdf = famdf %>% filter(!(id %in% validated_id))
      print(paste("#total", before,"#validated", before - nrow(famdf)))
    }
    df = parse_monomer_from_c(metdf, famdf)
    df$sample = fields[idx]
    df$lane = fields[2]
    dfall = rbind(dfall, df)
  }
  if (group_by_sample) {
    tmp = dfall %>% group_by(sample) %>% summarise(n = sum(n), den = sum(unique(den)), er = n/den, ref="X", alt="Y", muttype="Overall") %>% ungroup()
    dfall = dfall %>% group_by(alt, ref, sample) %>% summarise(n = sum(n), den = sum(unique(den)), er = n/den) %>% ungroup()
  } else {
    tmp = dfall %>% group_by(sample, lane) %>% summarise(n = sum(n), den = sum(unique(den)), er = n/den, ref="X", alt="Y", muttype="Overall") %>% ungroup()
  }
  
  dfall = dfall %>% mutate ( cano_alt = ifelse(ref == 'A', rc(alt), alt), cano_ref = ifelse(ref == 'A', 'T', ref)) %>% mutate ( cano_alt = ifelse(ref == 'G', rc(cano_alt), cano_alt), cano_ref = ifelse(ref == 'G', 'C', cano_ref))
  dfall = dfall %>% mutate(muttype = paste(cano_ref, cano_alt, sep=">"))
  dfall = dfall %>% select(-c("cano_alt", "cano_ref"))
  dfall = rbind(dfall, tmp)
  dfall$group = group 
  dfall = arrange(dfall, sample)
  compute_ci(dfall)
}

MonomerErrorRate <- function(folder, group, pattern = "mutation_metrics.txt$", sep='\\.', idx=1, validation_folder = NULL, group_by_sample=FALSE) {
  mono = monomer_fromc(folder, group, pattern = pattern, sep=sep, idx=idx, mutect2_folder=validation_folder, group_by_sample = group_by_sample)
  colnames(mono)[1:4] = c("ref_allele", "alt_allele", "total_base_errors", "total_bases_eval")
  mono
}

compute_ci <- function(df, cfl=0.95) {
  if ("total_base_errors" %in% names(df)) {
    dfconf = binom.confint(as.numeric(df$total_base_errors), as.numeric(df$total_bases_eval), methods="wilson", conf.level=cfl)
  }
  else {
    if ("n_snv" %in% names(df)) {
      dfconf = binom.confint(n=as.numeric(df$n_bases_eval), x=c(df$n_snv), method='wilson')
    } else {
      dfconf = binom.confint(as.numeric(df$n), as.numeric(df$den), methods="wilson", conf.level=cfl)
    }
  }
  dfnew = cbind(df, dfconf[,4:6])
  dfnew
}

rc <- function(base) {
  vec = c('A', 'T', 'C', 'G')
  names(vec) = c('T', 'A', 'G', 'C')
  vec[base]
}

canno_base <- function(df) {
 df %>% mutate ( alt = ifelse(ref == 'A', rc(alt), alt), ref = ifelse(ref == 'A', 'T', ref)) %>% mutate ( alt = ifelse(ref == 'G', rc(alt), alt), ref = ifelse(ref == 'G', 'C', ref))
}

add_muttype <- function(df) {
 if ("ref" %in% colnames(df)) {
  indels = df %>% filter(type != "SNV")
  indels = indels %>% mutate(cano_alt = alt, cano_ref = ref, muttype = type)
  df = break_mnv(df)
  df = df %>% mutate ( cano_alt = ifelse(ref == 'A', rc(alt), alt), cano_ref = ifelse(ref == 'A', 'T', ref)) %>% mutate ( cano_alt = ifelse(ref == 'G', rc(cano_alt), cano_alt), cano_ref = ifelse(ref == 'G', 'C', cano_ref))
 } else if ("ref_allele" %in% colnames(df)) {
  df = df %>% mutate ( cano_alt = ifelse(ref_allele == 'A', rc(alt_allele), alt_allele), cano_ref = ifelse(ref_allele == 'A', 'T', ref_allele)) %>% mutate ( cano_alt = ifelse(ref_allele == 'G', rc(cano_alt), cano_alt), cano_ref = ifelse(ref_allele == 'G', 'C', cano_ref))
 } else if("Reference_Allele" %in% colnames(df)) {
  indels = df %>% filter(Variant_Type == "DEL" | Variant_Type == "INS")
  indels = indels %>% mutate(cano_alt = Tumor_Seq_Allele2, cano_ref = Reference_Allele, muttype = Variant_Type)
  df = df %>% filter(Variant_Type == "SNP")
  df = df %>% mutate ( cano_alt = ifelse(Reference_Allele == 'A', rc(Tumor_Seq_Allele2), Tumor_Seq_Allele2), cano_ref = ifelse(Reference_Allele == 'A', 'T', Reference_Allele)) %>% mutate ( cano_alt = ifelse(Reference_Allele == 'G', rc(cano_alt), cano_alt), cano_ref = ifelse(Reference_Allele == 'G', 'C', cano_ref))
 }
 df = df %>% mutate(muttype = paste(cano_ref, cano_alt, sep=">"))
 df = rbind(df, indels)
 df %>% select(-c("cano_alt", "cano_ref"))
}

complement_mut <- function(mut) {
  nucs = unlist(strsplit(mut, ">"))
  paste0(rc(nucs[1]), ">", rc(nucs[2]))
}
complement_mut = Vectorize(complement_mut)

add_varid <- function(df) {
  if ("Chromosome" %in% colnames(df)) {
    df = df %>% mutate(id = paste(Chromosome, Start_Position, Reference_Allele, Tumor_Seq_Allele2, sep="_"))
  } else if ("chrom" %in% colnames(df)) {
    df = df %>% mutate(id = paste(chrom, ref_pos, ref, alt, sep='_'))
  } else if ("contig" %in% colnames(df)) {
    df = df %>% mutate(id = paste(contig, position, ref_allele, alt_allele, sep="_")) 
  }
  df
}

break_mnv <- function(mut_fam_df) {
  mut_fam_df = mut_fam_df %>% filter(type == "SNV")
  snv = mut_fam_df %>% filter(nchar(ref) == 1)
  mnv = mut_fam_df %>% filter(nchar(ref) > 1)
  mnv2snv = data.frame()
  if (nrow(mnv) > 0) {
    for (i in 1:nrow(mnv)) {
      for (j in 1:nchar(mnv[i,]$ref)) {
        tmp = mnv[i,]
        tmp$ref = substr(mnv[i,]$ref, j, j)
        tmp$alt = substr(mnv[i,]$alt, j, j)
        tmp$ref_pos = mnv[i,]$ref_pos + j - 1
        mnv2snv = rbind(mnv2snv, tmp)
      }  
    }
  }
  rbind(snv, mnv2snv) 
}


read_maf <- function(f, snp_only = F, nrows=Inf) {
  fields = c("Chromosome", "Start_Position", "Start_position" ,"End_Position", 
              "Variant_Type","Tumor_Sample_UUID", "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele2", "t_depth" ,"t_alt_count", "t_ref_count", "n_alt_count", "n_ref_count", "ccf_hat")
  # if (include_filter) {
  #   fields = c(fields, "FILTER")
  # }
  stdngs = fread(f, 
   skip="Chromosome", 
   quote="",
   nrows=nrows,
   select = fields)
  if ("ccf_hat" %in% colnames(stdngs)) {
    print("subset to clonal mutations with ccf 0.9")
    stdngs = stdngs %>% filter(ccf_hat > 0.9)
  }
  if (snp_only && "Variant_Type" %in% colnames(stdngs)) {
    stdngs %>% filter(Variant_Type == "SNP" | Variant_Type == "DNP")
  } else {
    stdngs 
  }
}

read_accu_files <- function(folder, pattern = "mutation_metrics.txt$", sep='\\.', idx=1, secidx=2, picard=F) {
  dfall = data.frame()
  for (f in list.files(folder, pattern=pattern)) {
    print(f)
    fp = file.path(folder,f)
    if(file.size(fp) == 0L) {
      print("skip empty file") 
      next
    }
    fields = unlist(str_split(f, sep))
    samp = fields[idx]
    if (picard) {
      df =fread(fp, skip=6)
    } else {
      df = fread(fp)
    }
    df$sample = samp
    df$panel = fields[secidx]
    dfall = rbind(dfall, df)
  }
  arrange(dfall, sample)
}


check_against_mutect2 <- function(codec, mutect2, atbase=TRUE) {
  if (is.character(codec) & str_detect(codec, ".maf$")) {
    codec = read_maf(codec)
    if (!("id" %in% colnames(codec))) {
      codec = add_varid(codec)
    }
    if ("FILTER" %in% colnames(codec)) {
      codec = codec %>% filter(!str_detect(FILTER, "alignment"))
    }
  } 
  maf = read_maf(mutect2, snp_only = F)  
  maf = add_varid(maf)
  if ("Variant_Type" %in% colnames(codec)) {
    codec  = codec %>% rename("Variant_Type"="type", "Start_Position"="ref_pos", "Chromosome"="chrom", "Reference_Allele"="ref", "Tumor_Seq_Allele2"="alt")
    codec = codec %>% mutate(type = ifelse(type == "SNP", "SNV", type))
    #codec  = codec %>% rename("type"="Variant_Type", "ref_pos"="Start_Position", "chrom"="Chromosome","ref"="Reference_Allele", "alt"="Tumor_Seq_Allele2")
  }
  if (atbase) {
    codec = codec %>% mutate(ref = ifelse(type == "INS" & ref !="-", "-", ref), alt = ifelse(type == "INS" & ref !="_", str_sub(alt, 2), alt))
    codec = codec %>% mutate(ref = ifelse(type == "DEL" & alt != "_", str_sub(ref, 2), ref), alt = ifelse(type == "DEL" & alt !="_", "-", alt), ref_pos = ifelse(type == "DEL" & alt != "_", ref_pos + 1, ref_pos))
  }
  codec = add_varid(codec)
  codec = codec %>% mutate(concordant = ifelse(id %in% (maf %>% pull(id)), "high", "low"))
  #codec_indels = codec %>% filter(type == "DEL" | type == "INS")
  maf = maf %>% mutate(concordant = ifelse(id %in% codec$id, TRUE, FALSE))
  list(codec=codec, mutect2=maf)
}


recall_prec_of_codec <- function(dfall, sample, include_filter = F, cosmic = cosmic_signatures_v3) {
  ## GET stat
  precdf = data.frame()
  muttype = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G", "DEL", "INS")
  confs = c("high", "low")
  if (include_filter) {
    confs = c(confs, "medium")
  }
  dfall$codec = add_muttype(dfall$codec)
  for (conf in confs) {
    dfc = data.frame(s=sample, mut = muttype, nerror = 0, concordant=conf)
    tmpdf = dfall$codec %>% filter(concordant == conf)
    if (nrow(tmpdf) > 1) {
   #   val = get_signa(tmpdf, sample, cosmic)
      for (mut in muttype) {
   #     dfc[dfc$mut == mut, ]$nerror = sum(val$mut_cnt[,str_detect(names(val$mut_cnt), mut)])
        dfc[dfc$mut == mut, ]$nerror = nrow(tmpdf %>% filter(muttype==mut))
      }
    }
    precdf = rbind(precdf, dfc)
  }
  if (include_filter) {
    dfall$mutect2 = dfall$mutect2 %>% filter(FILTER == "PASS")
  }
  
  recalldf= data.frame()
  for (conf in c(TRUE, FALSE)) {
    dfc = data.frame(s=sample, mut = muttype, nerror = 0, concordant=conf)
    tmpdf = dfall$mutect2 %>% filter(concordant == conf)
    tmpdf = add_muttype(tmpdf)
    if (nrow(tmpdf) > 1) {
      for (mut in muttype) {
        dfc[dfc$mut == mut, ]$nerror = nrow(tmpdf %>% filter(muttype==mut))
      }
    }
    recalldf = rbind(recalldf, dfc)
  }
  precdf = precdf %>% left_join( precdf %>% group_by(mut, s) %>% summarise(ntotal = sum(nerror)) %>% ungroup())
  precdf = rbind(precdf, precdf %>% filter(str_detect(mut, ">")) %>% group_by(s, concordant) %>% summarise(ntotal = sum(ntotal), nerror=sum(nerror), mut="SNV") %>% ungroup())
  precdf = precdf %>% mutate(perc = nerror/ntotal)
  
  recalldf = recalldf %>% left_join( recalldf %>% group_by(mut, s) %>% summarise(ntotal = sum(nerror)) %>% ungroup())
  recalldf = rbind(recalldf, recalldf %>% filter(str_detect(mut, ">")) %>% group_by(s, concordant) %>% summarise(ntotal = sum(ntotal), nerror=sum(nerror), mut="SNV") %>% ungroup())
  recalldf = recalldf %>% mutate(perc = nerror/ntotal)
  list(recall=recalldf, prec=precdf)
}


process_recall_prec <- function(codec_dir, 
                                mutect2_dir, 
                                codec_pat = "variants_called.txt$", 
                                mutect2_pat = "-filtered.annotated.maf",
                                #q60rate = 0.5,
                                atbase = TRUE,
                                sep='\\.', idx=1, subsample=1) {
  include_filter = F
  if (str_detect(mutect2_pat, "unfiltered.maf")) {
    include_filter = T
  }
  recall = data.frame() 
  prec = data.frame() 
  codecraw = data.frame()
  mutect2raw = data.frame()
  i = 0
  print(codec_dir)
  for (f in list.files(codec_dir, pattern=codec_pat)) {
    print(f)
    i = i + 1
    fields = unlist(str_split(f, sep))
    maf = file.path(mutect2_dir, paste0(fields[idx], mutect2_pat))
    print(paste(f, maf))
    if(str_detect(f, ".maf$")) {
      dat = file.path(codec_dir, f)
    } else {
      dat = fread(file.path(codec_dir, f))
      #dat %>% filter(numQ60/olen > q60rate)
      if ("filter" %in% colnames(dat)) {
        dat = dat %>% filter(filter == "PASS")
      }
      indels = dat %>% filter(type != "SNV")
      #snvs = break_mnv(dat)
      snvs = dat %>% filter(type == "SNV")
      dat = rbind(snvs, indels)
      dat = add_varid(dat) 
      print(nrow(dat))
      dat = dat[!duplicated(dat$id), ]
      print(nrow(dat))
      if (subsample < 1) {
        print(nrow(dat))
        dat = sample_n(dat, as.integer(nrow(dat) * subsample))
        print(nrow(dat))
      }
    }
    #sname = paste(fields[idx],fields[idx+1], sep=sep)
    sname = fields[idx]
    tmp = check_against_mutect2(dat, maf, atbase)
    tmp$codec$sample=sname
    tmp$mutect2$sample=fields[idx]
    print(paste(fields[idx],fields[idx+1])) 
    df = recall_prec_of_codec(tmp, sname, include_filter) 
    codecraw = rbind(codecraw, tmp$codec)
    mutect2raw = rbind(mutect2raw, tmp$mutect2)
    recall = rbind(recall, df$recall)
    prec = rbind(prec, df$prec)
  } 
  list(recall=recall, prec=prec, codecraw=codecraw, mutect2raw=mutect2raw)
}


isDNA = function(xx) {
  res = c()
  for (x in xx) {
    flag = TRUE
    for (i in 1:nchar(x)) {
      if (substring(x, i, i) %in% c("A", "C", "G", "T")) 
        next
      else {
        res = c(res, FALSE)
        flag = FALSE
        break 
      }
    }
    if(flag) res = c(res, TRUE)
  }
  res
}


CountMutInTriNuc <- function(mutfam, sample) {
  if (class(mutfam) == "character") {
    df = fread(mutfam) 
  } else {
    df = mutfam
  }
  if ("Chromosome" %in% colnames(df)) {
    if ("Variant_Type" %in% colnames(df)) {
      stdngs_snp = df %>% filter(Variant_Type == "SNP")
    } else {
      stdngs_snp = df
    }
    if ("filters" %in% colnames(stdngs_snp)) {
      print("remove sites failed filters")
      stdngs_snp = stdngs_snp %>% filter(filter == "PASS")
    }
    stdngs_snp$Chromosome = paste0("chr", stdngs_snp$Chromosome)
    if (!("Sample" %in% colnames(f))) {
      stopifnot(!is.null(sample))
      stdngs_snp$Sample = sample
    }
    sig_input <- mut.to.sigs.input(mut.ref = stdngs_snp, 
                                  sample.id = "Sample",
                                  chr = "Chromosome",
                                  pos = "Start_Position", 
                                  ref = "Reference_Allele",
                                  alt = "Tumor_Seq_Allele2")
    
  } else if("chrom" %in% colnames(df)) {
    df$Sample = sample
    df = df %>% filter(type == "SNV" & nchar(ref) == 1)
    if (!"id" %in% colnames(df)) {
      df = add_varid(df)
    }
    df = df[!duplicated(df$id), ] 
    #print(nrow(df))
    #print(df)
    #df$chrom = paste0("chr", df$chrom)
    sig_input <- mut.to.sigs.input(mut.ref = as.data.frame(df), 
                                  sample.id = "Sample",
                                  chr = "chrom",
                                  pos = "ref_pos", 
                                  ref = "ref",
                                  alt = "alt")
  } else if ("contig" %in% colnames(df)) {
    df$Sample = sample
    df = df %>% filter(pass_filter == 1)
    df = df %>% filter(nchar(alt_allele) == 1 & nchar(ref_allele) == 1)
    #df$chrom = paste0("chr", df$contig)
    sig_input <- mut.to.sigs.input(mut.ref = df, 
                                  sample.id = "Sample",
                                  chr = "contig",
                                  pos = "position", 
                                  ref = "ref_allele",
                                  alt = "alt_allele")
    
  } else {
    stopifnot(FALSE)
  }
  sig_input
}

SBS_errorrate <- function(mutfam, context, sample, base_count=T) {
  #dividing the occurrences of SBSs at a trinucleotide context by the total observations of that trinucleotide within the callable genome 
  tot_tri_counts = NULL
  ctxt = fread(context)
   
  if (!("site_count" %in% colnames(ctxt))) {
    colnames(ctxt) = c("context", "base_count")
  }
  ctxt = ctxt %>% drop_na(context)
  tri = ctxt %>% filter(isDNA(context)) %>% filter(nchar(context) == 3)
  strand_agno_tri = sapply(tri$context, function(x){
    if (substring(x, 2, 2) == 'A' || substring(x, 2, 2) == 'G') 
      insect::rc(x)
    else
      x
  })
  tri$strand_agno = strand_agno_tri[tri$context]
  if (base_count) {
    tritmp = tri %>% group_by(strand_agno) %>% summarise(n=sum(base_count))
  } else {
    tritmp = tri %>% group_by(strand_agno) %>% summarise(n=sum(site_count))
  }
  tricounts = as.numeric(tritmp$n)
  names(tricounts) = tritmp$strand_agno
 
  sig_input = CountMutInTriNuc(mutfam, sample) 
  print(sum(sig_input))
  from = sapply(names(sig_input), function(x){
    a = substring(x, c(1,3,7), c(1,3,7))
    paste0(a[1], a[2], a[3])
  })
  list(mut_cnt = sig_input, opportunities = tricounts[from], error_rate=sig_input/tricounts[from])
}


trinuc_b2nonb <- function(s) {
  if (str_detect(s, '\\[')) {
    tmp = unlist(str_split(s, ">"))
    a = paste0(str_replace(tmp[1], "\\[", ""), substr(tmp[2], 3, 3))
    b = paste0(substr(tmp[1], 1, 1), str_replace(tmp[2], "\\]", ""))
    return (paste(a, b, sep=">")) 
  }
  return (s)
}
trinuc_b2nonb = Vectorize(trinuc_b2nonb)

trinuc_nonb2b <- function(s) {
  if (str_detect(s, '\\[')) {
    return (s)
  }  else {
    mid = paste0(unlist(sapply(c(2,4,6), function(x) { substr(s, x, x)})), collapse = "")
    return (paste0(substr(s, 1,1), "[", mid, "]", substr(s, 3,3)))
  }
}
trinuc_nonb2b = Vectorize(trinuc_nonb2b)

get_signa <- function(f, sample = NULL, cosmic, context=NULL, method="deconstructSigs", mask=NULL, sig_cutoff=0.05) {
  rownames(cosmic) = str_replace(rownames(cosmic), " ", "_")
  decsig_cosmic = cosmic
  colnames(decsig_cosmic) = trinuc_nonb2b(colnames(decsig_cosmic))
  if (is.null(context)) {
    if (!("data.frame" %in% class(f)) && class(f) == "character" && str_detect(f, ".maf$")) {
      stdngs = fread(f, 
                     skip="Hugo_Symbol", 
                     quote="",
                     select = c("Chromosome", "Start_Position", "End_Position", 
                                "Variant_Type", "Tumor_Sample_Barcode", "Reference_Allele", "Tumor_Seq_Allele2", 
                                "failed_filters", "FILTER"))
      if ("Variant_Type" %in% colnames(stdngs)) {
        stdngs_snp = stdngs %>% filter(Variant_Type == "SNP")
      } else {
        stdngs_snp = stdngs
      }
      if ("failed_filters" %in% colnames(stdngs_snp)) {
        print("remove sites failed filters")
        stdngs_snp = stdngs_snp %>% filter(failed_filters == "none")
      }
      stdngs_snp$Chromosome = paste0("chr", stdngs_snp$Chromosome)
      stdngs_snp$sample = sample
      sig_input <- mut.to.sigs.input(mut.ref = stdngs_snp,
                                      sample.id = "sample", 
                                      chr = "Chromosome",
                                      pos = "Start_Position",
                                      ref = "Reference_Allele",
                                      alt = "Tumor_Seq_Allele2")
  
    } else if("data.frame" %in% class(f) && !"chrom" %in% colnames(f)) {
      if ("Chromosome" %in% colnames(f)) {
        stdngs = f 
        if ("Variant_Type" %in% colnames(stdngs)) {
          stdngs_snp = stdngs %>% filter(Variant_Type == "SNP")
        } else {
          stdngs_snp = stdngs
        }
        if ("failed_filters" %in% colnames(stdngs_snp)) {
          print("remove sites failed filters")
          stdngs_snp = stdngs_snp %>% filter(failed_filters == "none")
        }
        #stdngs_snp$Chromosome = paste0("chr", stdngs_snp$Chromosome)
        #stdngs_snp$Sample = sample
        sig_input <- mut.to.sigs.input(mut.ref = stdngs_snp,
                                        sample.id = "Tumor_Sample_Barcode", 
                                        chr = "Chromosome",
                                        pos = "Start_Position",
                                        ref = "Reference_Allele",
                                        alt = "Tumor_Seq_Allele2")
      } else if("contig" %in% colnames(f)) {
        f$Sample = sample
        print(f)
        sig_input <- mut.to.sigs.input(mut.ref =  as.data.frame(f),
                                        sample.id = "Sample", 
                                        chr = "contig",
                                        pos = "position",
                                        ref = "ref_allele",
                                        alt = "alt_allele")
        
      }
    } 
    else {
      sig_input <- CountMutInTriNuc(f, sample)
    }
    if (!is.null(mask)) {
      sig_input[,mask] = 0
    }
    mut_cnt = sig_input 
    print(paste("software default WGS TNC rate normalizaiton",sum(sig_input)))
    opportunities = NULL
    if (method=="deconstructSigs") {
      sigs = whichSignatures(tumor.ref = sig_input,
                            signatures.ref = as.data.frame(decsig_cosmic),
                            contexts.needed = TRUE,
                            signature.cutoff = sig_cutoff,	
                            tri.counts.method = "genome")
    } else if(method=="sigfit") {
      tmp = whichSignatures(tumor.ref = sig_input,
                            signatures.ref = as.data.frame(decsig_cosmic),
                            contexts.needed = TRUE,
                            tri.counts.method = "genome")
      mcmc_samples_refit <- fit_signatures(counts = sig_input,
                                           opportunities = "human-genome", 
                                           signatures = cosmic,
                                           iter = 10000,
                                           chain=1,
                                           warmup = 2000)
      exposures <- retrieve_pars(mcmc_samples_refit, "exposures")
      sigs = list(tumor = tmp$tumor, weights = exposures$mean)
    } else {
      stopifnot(FALSE, "Unknown method")
    }
      
  } else {
    #stopifnot(file.exists(context))
    sbs = SBS_errorrate(f, context, sample)
    if (!is.null(mask)) {
      sbs$error_rate[, mask] = 0
      sbs$mut_cnt[, mask] = 0
    }
    print("normalize by coverage")
    if (method == "deconstructSigs") {
      sigs = whichSignatures(tumor.ref = sbs$error_rate / sum(sbs$error_rate),
                              signatures.ref = as.data.frame(decsig_cosmic),
                              signature.cutoff = sig_cutoff,
                              contexts.needed = FALSE)
    } else if (method == "sigfit") {
      mcmc_samples_refit <- fit_signatures(counts = sbs$mut_cnt,
                                           opportunities = sbs$opportunities, 
                                           signatures = cosmic,
                                           iter = 10000,
                                           chain=1,
                                           warmup = 2000)
      exposures <- retrieve_pars(mcmc_samples_refit, "exposures")
      sigs = list(weights = exposures$mean)
    }
    mut_cnt = sbs$mut_cnt
    opportunities = as.data.frame(matrix(sbs$opportunities,nrow=1))
    rownames(opportunities) = rownames(mut_cnt)
    colnames(opportunities) = names(sbs$opportunities)
  }
  list(sigs=sigs, mut_cnt=mut_cnt, opportunities=opportunities) 
}
  
pairwise_cossim <- function(mat, mat2=NULL) {
  if (is.null(mat2)) {
    result = matrix(NA, nrow=nrow(mat), ncol=nrow(mat))
    for (i in 1:nrow(mat)) {
      for (j in i:nrow(mat)) {
        result[i,j] = sigfit::cosine_sim(as.numeric(mat[i,]), as.numeric(mat[j,]))
      } 
    }
    result
  } else {
      res = data.frame(left=rownames(mat), right=rownames(mat2), cossim=NA)
    if (nrow(mat2) == 1) {
      mat2 = mat2[,colnames(mat)]
      for (i in 1:nrow(mat))  {
        res[i,]$cossim = sigfit::cosine_sim(as.numeric(mat[i,]), as.numeric(mat2)) 
      }
       
    } else {
      stopifnot(dim(mat) == dim(mat2))
      mat2 = mat2[,colnames(mat)]
      for (i in 1:nrow(mat))  {
        res[i,]$cossim = sigfit::cosine_sim(as.numeric(mat[i,]), as.numeric(mat2[i,])) 
      }
    }
    res
  }
}

plot_weights <- function(df, good_sigs, color_lab, legend=T, min_weight = 0.05, constitute_sigs=c("SBS1", "SBS3", "SBS4")) {
  if (legend) {
    plotlegend = theme(legend.position="bottom", legend.box="vertical", legend.margin=margin())
  } else {
    plotlegend = theme(legend.position = "none")
  }
  df$sample = rownames(df)
  plotdf = df %>% pivot_longer(setdiff(colnames(df), c("sample")), names_to = c("signature"), values_to = c("weight"))
  plotdf[plotdf$weight < min_weight, ]$weight = 0
  plotdf$signature = factor(plotdf$signature, levels=colnames(df))
  plotdf$good = FALSE
  plotdf[plotdf$signature %in% good_sigs,]$good = TRUE
  plotdf %>% filter(weight >0 | signature %in% c("SBS1", "SBS3", "SBS4")) %>%
  ggplot(aes(y=sample, x = signature, color = good, size = weight)) + 
  geom_point() +
  cowplot::theme_cowplot(22) + 
  theme(axis.line  = element_blank()) +
  theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  ylab('') + xlab('') +
  scale_size_continuous(limits=c(min_weight, 1)) +
  scale_color_manual(name=color_lab, values = c("gray76", "firebrick2")) +
  theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_line( size=.1, color="black" )) +
  plotlegend
}

print_weights <- function(df) {
  df[, colSums(df)>0]
}

process_folder <- function(d, pattern, sep, idx, cosmic, method = "deconstructSigs", context_count = F, mask=NULL, sigcutoff = 0.05) {
  # d: directory
  # pattern: file name pattern
  # label: output label
  # alignment_check_d
  # sep, idx: split the file name by seperator and idx of the coverage in the splitted vector
  all_weights = data.frame()
  mutdist = data.frame()
  nmut= data.frame()
  oppo = NULL 
  rawdata = list()
  i = 1
  for (f in list.files(d, pattern)) {
    print(f)
    fields = unlist(str_split(f, sep))
    sample = fields[idx]
    dat = file.path(d, f)
    if (str_detect(f, ".maf$")) {
      tmpdf = read_maf(dat, nrows=2)
      if (nrow(tmpdf) == 0) {
        print(paste(f, "has no mutations"))
        next
      }
    }
    if (context_count) {
      if (idx > 1) {
        ctx_f = fields[1] 
        for (xx in 2:idx) {
          ctx_f = paste(ctx_f, fields[xx], sep=".")
        }
        ctx_f = paste(ctx_f, "context_count","txt", sep=".")
      } else {
        ctx_f = paste(sample, "context_count", "txt", sep=".")
      }
      print(file.path(d, ctx_f))
      res = get_signa(dat, sample, cosmic, file.path(d, ctx_f), method=method, mask=mask, sig_cutoff = sigcutoff)
    } else {
      res = get_signa(dat, sample, cosmic, method=method, mask=mask, sig_cutoff = sigcutoff)
    }
    oppo =  rbind(oppo, res$opportunities)
    w = weights(res$sigs)
    all_weights = rbind(all_weights, w) 
    mutdist =  rbind(mutdist, res$sigs$tumor)
    nmut = rbind(nmut, res$mut_cnt)
    rawdata[[i]] = res
    i = i+ 1
  }
  list(weights = all_weights, mutdist = mutdist, nmut=nmut, oppo=oppo, rawdata=rawdata)
}

##simple codec signature analysis
### 

get_sigs <- function(dat, cosmic, min_del = 1, min_del_frag = 30, sig_cutoff=0.05, method="deconstructSigs") {
  weights= data.frame()
  nmut = data.frame()
  sigs = data.frame()
  indels_all = data.frame()
  indels_count = data.frame()
  for (ss in unique(dat$sample)) {
    sig = get_signa(dat %>% filter(sample == ss), sample=ss, cosmic = cosmic, sig_cutoff = sig_cutoff, method=method)
    #indels = dat %>% filter(sample == ss) %>% filter((type == "DEL" & nchar(ref) > min_del & flen >= min_del_frag) | (type == "INS"))
    #indels = indels %>% select(chrom, ref_pos, ref, alt, type)
    #colnames(indels) = c("chr", "position", "REF", "ALT", "type")
    #indels = indels %>% mutate(type = ifelse(type == "DEL", "D", type )) %>% mutate(type = ifelse(type == "INS", "I", type ))
    #indels_table = tabToIndelsClassification(indels, ss)
    #indels_table$indels_classified$sample = ss
    #indels_all = rbind(indels_all, indels_table$indels_classified)
    weights = rbind(weights, sig$sigs$weights)
    nmut = rbind(nmut, sig$mut_cnt)
    sigs = rbind(sigs, sig$sigs$tumor)
    #indels_count = rbind(indels_count, indels_table$count_proportion)
  }
  #list(spec=sigs, nmut=nmut, weights=weights, indel=indels_all, indel_counts=indels_count) 
  list(spec=sigs, nmut=nmut, weights=weights) 
} 

#' Initial coercion to matrix for signatures/exposures/counts
to_matrix <- function(x, int = FALSE) {
  # If x is coming from retrieve_pars, get mean
  if (is.list(x) & "mean" %in% names(x))
    x <- x$mean
  # If x is a vector, transform to 1-row matrix
  if (is.vector(x))
    x <- matrix(x, nrow = 1, dimnames = list(NULL, names(x)))
  # Otherwise, try coercing to matrix
  if (!is.matrix(x))
    x <- as.matrix(x)
  # For counts matrix: if real-valued, round
  if (int) {
    x <- round(x)
  }
  x
}

mut_types <- function(strand = FALSE) {
  bases <- c("A", "C", "G", "T")
  muts <- paste0(rep(rep(bases, each = 4), 6),
                 rep(bases[c(2, 4)], each = 48),
                 rep(bases, 6 * 16 / 4),
                 ">",
                 rep(rep(bases, each = 4), 6),
                 c(rep(bases[-2], each = 16), rep(bases[-4], each = 16)),
                 rep(bases, 6 * 16 / 4))
  if (strand) {
    paste(c(rep("T", 96), rep("U", 96)), muts, sep = ":")
  }
  else {
    muts
  }
}


plot_signature <- function (spectra, pdf_path = NULL, pdf_width = 24, pdf_height = 8, 
                            name = "", max_y = NULL, colors = NULL, boxes = TRUE, xticklab=FALSE, ylab=FALSE, mutlab=FALSE) 
{
  if (is.list(spectra) & "mean" %in% names(spectra)) {
    spec <- to_matrix(spectra$mean)
    lwr <- to_matrix(spectra$lower)
    upr <- to_matrix(spectra$upper)
  }
  else {
    colnames(spectra) <-  trinuc_b2nonb(colnames(spectra))
    spectra <-  spectra[, mut_types()]
    spec <- to_matrix(spectra)
    lwr <- NULL
    upr <- NULL
  }
  NCAT <- ncol(spec)
  NSAMP <- nrow(spec)
  strand <- NCAT == 192
  counts <- any(spec > 1)
  TYPES <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
  COLORS <- c("deepskyblue", "black", "firebrick2", "gray76", 
              "darkolivegreen3", "rosybrown2")
  STRANDCOL <- c("deepskyblue3", "red3")
  BACKCOL <- c("#00BFFF33", "#00000033", "#EE2C2C33", "#C2C2C24D", 
               "#A2CD5A4D", "#EEB4B44D")
  LINECOL <- "gray60"
  XL <- c(0.2, 19.4, 38.6, 57.8, 77, 96.2)
  XR <- c(19.2, 38.4, 57.6, 76.8, 96, 115.2)
  BACKLIM <- c(0, 46.8, 93.2, 139.55, 186, 232.35, 279.2)
  if (!is.null(pdf_path)) {
    cairo_pdf(pdf_path, width = pdf_width, height = pdf_height, 
              onefile = TRUE)
    if (ncol(spec) %in% c(96, 192)) {
      par(mar = c(5.5, 7, 7.5, 2))
    }
    else {
      par(mar = c(10, 7, 7.5, 2))
    }
  }
  if (!(ncol(spec) %in% c(96, 192))) {
    if (is.null(colnames(spec))) {
      types <- paste("Mut. type", 1:ncol(spec))
    }
    else {
      types <- colnames(spec)
    }
    for (i in 1:NSAMP) {
      if (is.null(max_y)) {
        FACTOR <- 1.05
        samp_max_y <- ifelse(is.null(upr), max(spec[i, 
        ]) * FACTOR, max(upr[i, ]) * FACTOR)
      }
      else {
        samp_max_y <- max_y
      }
      if (is.null(colors)) {
        colors = "orangered3"
      }
      else {
        if ((length(colors) > 1) & (length(colors) != 
                                    ncol(spec))) {
          stop("'colors' must contain either a single value, or one value per mutation type.")
        }
      }
      bars <- barplot(spec[i, ], names.arg = types, col = colors, 
                      mgp = c(3, 0.8, 0), border = "white", las = 2, 
                      cex.names = 1, ylim = c(0, samp_max_y), yaxt = "n")
      if (counts) {
        axis(side = 2, cex.axis = 1.9, lwd = 2)
        label <- "Mutations"
        n_text <- paste0(" (", prettyNum(sum(spec[i, 
        ]), big.mark = ","), " mutations)")
      }
      else {
        axis(side = 2, cex.axis = 1.9, lwd = 2)
        label <- "Mutation probability"
        n_text <- ""
      }
      if (is.null(name)) {
        nme <- rownames(spec)[i]
      }
      else {
        nme <- name
      }
      mtext(label, side = 2, cex = 2.4, line = 3.5)
      title(paste0(nme, n_text), line = 4, cex.main = 2.5)
      if (!is.null(lwr)) {
        arrows(bars, upr[i, ], bars, lwr[i, ], length = 0, 
               lwd = 3, col = LINECOL)
      }
      if (boxes) {
        box(lwd = 2)
      }
    }
  }
  else {
    if (!strand) {
      for (i in 1:NSAMP) {
        if (is.null(max_y)) {
          FACTOR <- 1.095
          samp_max_y <- max(0.05, ifelse(is.null(upr), 
                                         max(spec[i, ]) * FACTOR, max(upr[i, ]) * 
                                           FACTOR))
        }
        else {
          samp_max_y <- max_y
        }
        if (boxes) {
          xlim <- c(-0.105, 115.5)
        }
        else {
          xlim <- c(-1, 116)
        }
        xtickshow = ifelse(xticklab, "s", "n")
        bars <- barplot(spec[i, ], names.arg = substr(mut_types(), 
                                                      1, 3), mgp = c(3, 0.8, 0), col = rep(COLORS, 
                                                                                           each = 16), border = "white", las = 2, ylim = c(0, 
                                                                                                                                           samp_max_y), xlim = xlim, yaxt = "n", cex.names = 1.6, 
                        xaxs = "i", family = "mono", xaxt=xtickshow) #ruolin
        if (xticklab) {
          for (j in 1:length(COLORS)) {
            idx <- ((j - 1) * 16 + 1):(j * 16)
            axis(side = 1, at = bars[idx], tick = FALSE,
                 cex.axis = 1.6, mgp = c(3, 0.8, 0), las = 2,
                 family = "mono", font = 2, col.axis = COLORS[j],
                 labels = paste0(" ", substr(mut_types()[idx],
                                             2, 2), " "))
          }
        }
        if (counts) {
          axis(side = 2, cex.axis = 1.9, lwd = 2)
          label <- "Mutations"
          n_text <- paste0(" (", prettyNum(sum(spec[i, 
          ]), big.mark = ","), " mutations)")
        }
        else {
          axis(side = 2, at = seq(0, samp_max_y, samp_max_y/5), #ruolin 
               cex.axis = 1.9, lwd = 2)
          label <- "Mutation probability"
          n_text <- ""
        }
        if (!ylab) label = ""
        if (is.null(name)) {
          nme <- rownames(spec)[i]
        }
        else {
          nme <- name
        }
        mtext(label, side = 2, cex = 2.4, line = 3.5)
        title(paste0(nme, n_text), line = 4, cex.main = 2.5)
        if (!is.null(lwr)) {
          arrows(bars, upr[i, ], bars, lwr[i, ], length = 0, 
                 lwd = 3, col = LINECOL)
        }
        if (mutlab) {
          text(x = (XL + XR)/2, y = 1.055 * samp_max_y, 
               labels = TYPES, cex = 2.4, xpd = TRUE)
          rect(xleft = XL, xright = XR, ybottom = 0.945 * 
                 samp_max_y, ytop = samp_max_y, col = COLORS, 
               border = "white")
        }
        if (boxes) {
          box(lwd = 2)
        }
      }
    }
    else {
      for (i in 1:NSAMP) {
        if (is.null(max_y)) {
          FACTOR <- 1.095
          samp_max_y <- max(0.05, ifelse(is.null(upr), 
                                         max(spec[i, ]) * FACTOR, max(upr[i, ]) * 
                                           FACTOR))
        }
        else {
          samp_max_y <- max_y
        }
        if (boxes) {
          xlim <- c(0, 279.2)
        }
        else {
          xlim <- c(-3, 280)
        }
        barplot(rbind(spec[i, 1:(NCAT/2)], spec[i, (NCAT/2 + 
                                                      1):NCAT]), beside = TRUE, col = NA, border = NA, 
                space = c(0.1, 0.8), xaxs = "i", yaxt = "n", 
                xaxt = "n", ylim = c(0, samp_max_y), xlim = xlim)
        for (j in 1:length(COLORS)) {
          rect(xleft = BACKLIM[j], xright = BACKLIM[j + 
                                                      1], ybottom = 0, ytop = samp_max_y, col = BACKCOL[j], 
               border = "white")
          text(x = (BACKLIM[j] + BACKLIM[j + 1])/2, y = 1.055 * 
                 samp_max_y, labels = TYPES[j], cex = 2.4, 
               xpd = TRUE)
          rect(xleft = BACKLIM[j], xright = BACKLIM[j + 
                                                      1], ybottom = 0.945 * samp_max_y, ytop = samp_max_y, 
               col = COLORS[j], border = "white")
        }
        legend("topright", bty = "n", inset = c(0.016, 
                                                0.03), legend = c("Transcribed", "Untranscribed"), 
               cex = 2.1, fill = NA, border = NA)
        legend("topright", bty = "n", inset = c(0.115, 
                                                0.03), pch = 15, pt.cex = 3.75, col = STRANDCOL, 
               legend = c("", ""), cex = 2.1)
        bars <- barplot(rbind(spec[i, 1:(NCAT/2)], spec[i, 
                                                        (NCAT/2 + 1):NCAT]), names.arg = substr(mut_types(), 
                                                                                                1, 3), beside = TRUE, space = c(0.1, 0.8), 
                        mgp = c(3, 0.8, 0), las = 2, col = STRANDCOL, 
                        border = "white", yaxt = "n", cex.names = 1.6, 
                        xaxs = "i", family = "mono", add = TRUE)
        for (j in 1:length(COLORS)) {
          idx <- ((j - 1) * 16 + 1):(j * 16)
          axis(side = 1, tick = FALSE, at = colMeans(bars[, 
                                                          idx]), cex.axis = 1.6, mgp = c(3, 0.8, 0), 
               las = 2, family = "mono", font = 2, col.axis = COLORS[j], 
               labels = paste0(" ", substr(mut_types()[idx], 
                                           2, 2), " "))
        }
        if (counts) {
          axis(side = 2, cex.axis = 1.9, lwd = 2)
          label <- "Mutations"
          n_text <- paste0(" (", prettyNum(sum(spec[i, 
          ]), big.mark = ","), " mutations)")
        }
        else {
          axis(side = 2, at = seq(0, samp_max_y, 0.05), 
               cex.axis = 1.9, lwd = 2)
          label <- "Mutation probability"
          n_text <- ""
        }
        if (is.null(name)) {
          nme <- rownames(spec)[i]
        }
        else {
          nme <- name
        }
        if (NSAMP > 1) {
          num <- paste0(" #", i)
        }
        else {
          num <- ""
        }
        mtext(label, side = 2, cex = 2.4, line = 3.5)
        title(paste0(nme, n_text), line = 4, cex.main = 2.5)
        if (!is.null(lwr)) {
          bars <- as.numeric(t(bars))
          arrows(bars, upr[i, ], bars, lwr[i, ], length = 0, 
                 lwd = 2.5, col = LINECOL)
        }
        if (boxes) {
          box(lwd = 2)
        }
      }
    }
  }
  if (!is.null(pdf_path)) {
    invisible(dev.off())
  }
}
