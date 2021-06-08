################
######### Figure 1B
####################

prepare_aggresult = function(dir, pattern, fsize_idx = 3) {
  dir = reformat_on_mac(dir)
  result = data.frame()
  family_sizes_files = list.files(dir, pattern = pattern )
  for (fsf in family_sizes_files) {
    fields = unlist(str_split(fsf, "_"))
    tb = as_tibble(fread(file.path(dir, fsf)))
    tb$fsize = as.integer(fields[fsize_idx])
    result = rbind(result, tb)
  }
  result$err = result$total_base_errors / result$total_bases_eval
  result
}

prepare_error <- function(folder, group, group_size =2, adjust = F) {
  dfall = data.frame()
  i = 0
  for (f in list.files(reformat_on_mac(folder), pattern="error_metrics.txt$")) {
    df = fread(reformat_on_mac(file.path(folder, f)))
    if (adjust) {
      adj_f = str_replace(f, "_error_metrics.txt", "_mutant_families.txt")
      adj_df = fread(reformat_on_mac(file.path(folder, adj_f)))
      df$total_base_errors = nrow(adj_df %>% filter(pass_filter == 1 & dist_from_end_of_frag != 245))
      df$errors_per_base_sequenced = df$total_base_errors / df$total_bases_eval
    }
    df$sample = unlist(str_split(f, "_"))[2]
    df$sid = i %% group_size
    dfall = rbind(dfall, df)
    i = i + 1
  }
  dfall$group = group 
  dfall
}


####Read Level Downsample Anlysis
dir = "/xchip/bloodbiopsy/ruolin/convention_duplex/accu_by_family_size/HD82/wrk/Miredas_new"
#dir = "/xchip/bloodbiopsy/ruolin/convention_duplex/accu_by_family_size/FFPE/wrk/Miredas"
#dir = "/xchip/bloodbiopsy/ruolin/convention_duplex/accu_by_family_size/early_BC_gdna/wrk/Miredas_new"

dsc = prepare_aggresult(dir, pattern = "DSC_fsize.*error_by_nuc_metrics.agg.txt")
dsc = dsc %>% filter(fsize != 1 & !is.na(fsize))
dsc_by_nuc = dsc %>% group_by(ref_allele, alt_allele, fsize) %>% summarise(sum_total_errors = sum(total_base_errors), sum_total_bases = sum(total_bases_eval)) 
dsc_by_nuc = dsc_by_nuc %>% mutate(monomer=paste(ref_allele, alt_allele, sep='>'))
dsc_by_nuc = dsc_by_nuc %>% mutate(err = ifelse(sum_total_errors == 0, 1/sum_total_bases, sum_total_errors/sum_total_bases), type = ifelse(sum_total_errors == 0, "DSC_theor_UL", "DSC_actual"))

ssc = prepare_aggresult(dir, pattern = "SSC_fsize.*error_metrics.agg.txt")
ssc = ssc %>% filter(!is.na(fsize))
ssc = ssc %>% filter(fsize %in% seq(2,20, 2))
ssc = ssc %>% mutate(fsize = fsize / 2)
ssc%>% filter (fsize == 1) %>% summarise(median = median(err))
ssc$monomer = "SSC"

ssc_by_nuc = ssc %>% group_by(ref_allele, alt_allele, fsize, sample_id ) %>% summarise(sum_total_errors = sum(total_base_errors), sum_total_bases = sum(total_bases_eval)) 
ssc_by_nuc = ssc_by_nuc %>% mutate(monomer=paste(ref_allele, alt_allele, sep='>'))
ssc_by_nuc = ssc_by_nuc %>% mutate(err = ifelse(sum_total_errors == 0, 1/sum_total_bases, sum_total_errors/sum_total_bases))
ssc_by_nuc$type = "SSC"
dsc_by_nuc$sample_id = "total"
all_by_nuc = rbind(dsc_by_nuc, ssc_by_nuc)

library(scales)
library(ggbeeswarm)
plot_theme <- theme_bw() + theme(text = element_text(size=10, color = "black"),
                                 axis.text = element_text(size = 10, color = "black"),
                                 axis.text.y = element_text(margin = margin(r=10)),
                                 axis.title = element_text(size = 10, color = "black"),
                                 title = element_text(size = 10, color = "black"),
                                 legend.text = element_text(size = 10, color = "black"),
                                 legend.title = element_text(size = 10, color = "black"),
                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank())
theme_set(plot_theme)

xlabs = as.character(1:20)
xlabs[seq(2,20,2)] = ""

ggplot(data=dsc_by_nuc, aes(x=factor(fsize, levels=1:20, labels=letters[1:20]), y=err)) + 
  geom_smooth(data = dsc_by_nuc, aes(group=monomer, color=monomer), span=0.4, show.legend = F, size=0.3, se=FALSE) + 
  geom_jitter(aes(color=monomer, shape=type), size = 0.3) + 
  geom_boxplot(data = ssc, aes(x=factor(fsize, levels=1:20, labels=letters[1:20]), y=err, fill=monomer), size=0.1,  outlier.size=0.1) +
  xlab("# reads per consensus sequence") + ylab("Error Rate") +
  scale_y_log10(limits = c(1e-8, 1e-3), breaks = c(1e-3,1e-4,1e-5,1e-6,1e-7, 1e-8), labels = c("1/1k", "1/10k", "1/100k", "1/1M", "1/10M", "1/100M")) + 
  #annotation_logticks(sides = "l", short = unit(0.2, "cm"), mid = unit(0.35, "cm"), long = unit(0.5, "cm")) +
  annotation_logticks(sides = "l", outside = T) +  coord_cartesian(clip = "off") + 
  scale_x_discrete(labels = xlabs) + 
  scale_shape_manual(values = c(19,1), guide = "none")+ labs(color = "", fill="") + 
  guides(colour = guide_legend(override.aes = list(size=3))) + 
  theme(legend.position="bottom")
ggsave("~/Documents/RO1.pdf", width = 4, height = 3)

###
#############Figure 3A
library(binom)
cds_dir = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_newsplit_highconf/miredas/fsize2"
cds_all = prepare_aggresult(cds_dir, pattern = "*error_metrics.txt")
cds_all$fsize = 2
dup_all_dir = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_dup/miredas/duplex"
dup_all = prepare_error(dupdsc_error, "DSC", adjust=T)
dup_all$fsize = 2
ssc_all_dir = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_newlinker_singleinsert/miredas/fsize2"
ssc_all = prepare_aggresult(ssc_all_dir, pattern = "*error_metrics.txt")
ssc_all$fsize = 2
singlepair_pat_dir = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_newlinker_singleinsert/fsize1_erate/miredas/CDS_05246_315_DC/"
singlepair_hd_dir = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_newlinker_singleinsert/fsize1_erate/miredas/CDS_HD82_DC/"
singlepair_pat = prepare_aggresult(singlepair_pat_dir, pattern = "*error_metrics.txt")
singlepair_hd_dir = prepare_aggresult(singlepair_hd_dir, pattern = "*error_metrics.txt")
singlepair_all = rbind(singlepair_hd_dir, singlepair_pat)

df = data.frame(group=c("NoOverlap", "CDS", "DupSeq", "PairMerge"), total_base_errors = rep(0,4), total_bases_eval = rep(0,4))
df[df$group == "CDS", "total_base_errors"] = sum(cds_all$total_base_errors)
df[df$group == "CDS", "total_bases_eval"] =  sum(cds_all$total_bases_eval)
df[df$group == "DupSeq", "total_bases_eval"] = sum(dup_all$total_bases_eval)
df[df$group == "DupSeq", "total_base_errors"] = sum(dup_all$total_base_errors)
df[df$group == "NoOverlap", "total_bases_eval"] = 152126578
df[df$group == "NoOverlap", "total_base_errors"] = 158343
df[df$group == "PairMerge", "total_bases_eval"] = sum(singlepair_all$total_bases_eval)
df[df$group == "PairMerge", "total_base_errors"] = sum(singlepair_all$total_base_errors)


dfconf = binom.confint(df$total_base_errors, df$total_bases_eval, method="wilson")
dfnew = cbind(df, dfconf[,4:6])
dfnew

ggplot(dfnew, aes(x = group, y = mean, color = group)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
  scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3), labels = c("1/10M", "1/1M", "1/100k", "1/10k", "1/1k")) +
  annotation_logticks(sides = "l", outside = T) +  coord_cartesian(clip = "off", ylim = c(1e-7, 1e-3)) + 
  labs(x = "", y = "Error Rate", color = "", shape = "") +
  theme(legend.position = "none")
ggsave("~/Documents/ErrorRate.pdf", width = 4, height = 3)


##
#######Figure 3C
highconf_detect = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_newsplit_highconf/detect_fingerprints"
duplex = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_dup/detect_fingerprints"
prepare_detect <- function(folder, group) {
  dfall = data.frame()
  i = 1
  for (f in list.files(reformat_on_mac(folder), pattern="detect_fingerprint.txt$")) {
    df = fread(reformat_on_mac(file.path(folder, f)))
    df = df %>% filter(n_alt_pass > 0)
    if (nrow(df) == 0) next
    df$sample = unlist(str_split(f, "_"))[2]
    df$sid = i %% 2 + 1
    dfall = rbind(dfall, df)
    i = i + 1
  }
  dfall$group = group 
  dfall
}
a = prepare_detect(highconf_detect, "CDS")
b = prepare_detect(duplex, "Duplex")
stat_box_data <- function(y) {
  return( 
    data.frame(
      y = median(y) *0.9 ,
      label = paste('count =', length(y), '\n',
                    'mean =', round(mean(y), 1), '\n')
    )
  )
}
detectdf = rbind(a,b)
plotdf = detectdf %>% group_by(sample, group, id) %>% summarise(saf = sum(n_alt_pass)/sum(depth_pass)) 

ggplot(data =  plotdf %>% filter(sample == "05246"), aes(x = group, y = saf, fill=sample)) + geom_boxplot() + 
  #stat_summary(fun.data = stat_box_data, geom = "text",hjust = 0.5, vjust = 0.9,  position = position_dodge(width = 0.75)) + 
  ylim(0,0.5) + 
   theme(axis.title.x=element_blank()) +
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  theme(legend.position = "none") +
  ylab("Somatic VAF") +   theme(legend.title=element_blank()) 
ggsave("~/Documents/recall.pdf", width = 4, height = 3)

ggplot(data = plotdf, aes(x=sample, fill=group)) + geom_bar(stat="count", position = position_dodge2())
