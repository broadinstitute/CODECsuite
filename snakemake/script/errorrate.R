
###Error position Evaluation

eof_error_table <- function(inputdir) {
  cdspan = prepare_error(inputdir, 'trim_eof')
  cdspan_noeof = prepare_error(paste0(inputdir, "/noeof"), 'notrim_eof')
  newd = rbind(cdspan, cdspan_noeof) %>% group_by(sample_id) %>% mutate(eof12_nerror = total_base_errors[group == 'notrim_eof'] - total_base_errors) %>%
    mutate(eof12_ntotal = total_bases_eval[group == 'notrim_eof'] - total_bases_eval) %>% filter(group == "trim_eof") %>% select("sample_id", "sample", "lane", "total_base_errors", "total_bases_eval", "eof12_nerror", "eof12_ntotal")
  colnames(newd) = c("sample_id", "sample", "lane", "trim_eof_nerror", "trim_eof12bp", "within_eof_nerror", "within_eof12bp")
  df = gather(newd %>% ungroup() %>% select(-c(trim_eof_nerror, within_eof_nerror)), group, total_bases_eval, trim_eof12bp, within_eof12bp) %>% mutate(total_base_errors =
        gather(newd %>% select(-c(trim_eof12bp, within_eof12bp)), group, total_base_errors, trim_eof_nerror, within_eof_nerror) %>% pull(total_base_errors))
  compute_ci(df)
}

dupseq_df = eof_error_table(dupseq)

withET_df = eof_error_table(cds)

total_df = rbind(dupseq_df, withET_df)

ggplot(total_df %>% filter(sample == "HD_82_ctDNA"), aes(x = group, y = mean, color = lane)) +
  geom_point(size = 1, position = position_dodge2(width=0.5)) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, position = position_dodge2(width=0.5)) +
  scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3), labels = c("1/10M","1/1M", "1/100k", "1/10k", "1/1k")) +
  annotation_logticks(sides = "l", outside = T) +  coord_cartesian(clip = "off", ylim = c(1e-7, 1e-3)) +
  labs(x = "", y = "Error Rate", color = "", shape = "")

library(ggplot2)

### Overall error rate
plot_theme <- theme_bw() + theme(text = element_text(size=10, color = "black"),
                                 axis.text = element_text(size = 10, color = "black"),
                                 axis.text.y = element_text(margin = margin(r=10)),
                                 axis.title = element_text(size = 10, color = "black"),
                                 title = element_text(size = 10, color = "black"),
                                 legend.text = element_text(size = 10, color = "black"),
                                 legend.title = element_text(size = 10, color = "black"),
                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank())
theme_set(plot_theme)


# ##### CDS PanCancerdata Error rate
header = c("sample_id", "total_base_errors", "total_bases_eval", "sample", "lane", "group")

total_df = rbind(dupseqdf, cdsdf, rawdf, r1r2df, ssc) 
total_df = total_df %>% mutate(group = paste(lane, group, sep="_"))
total_df = compute_ci(total_df)

ggplot(total_df, aes(x = group, y = mean, color = sample)) +
  geom_point(size = 1, position=position_dodge2(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, position=position_dodge2(width=1)) +
  scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3), labels = c("1/10M", "1/1M", "1/100k", "1/10k", "1/1k")) +
  annotation_logticks(sides = "l", outside = T) +  coord_cartesian(clip = "off", ylim = c(1e-7, 1e-3)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  labs(x = "", y = "Error Rate", color = "", shape = "") 
  #theme(legend.position = "none")
ggsave("~/Documents/PanCancerTotalErrorRate.pdf", width = 6, height = 4.5, dpi='retina')
write_csv(total_df %>% arrange(lane), "~/Documents/tmp.csv")

###
######### Downsample plot
#####################

total_ds = rbind(cdshd82_ds, cdspat_ds)
total_ds = compute_ci(total_ds)  
total_ds$family_size = substr(total_ds$sample_id, nchar(total_ds$sample_id), nchar(total_ds$sample_id))
total_ds$sample = 'cfDNA'
total_ds = total_ds %>% select(-c(lane, errors_per_base_sequenced))
colnames(total_ds)[5] = 'dna_type'
colnames(total_ds)[6] = 'sample'
ggplot(total_ds, aes(x = family_size, y = mean, color = sample)) +
  geom_point(size = 1, position=position_dodge2(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, position=position_dodge2(width=1)) +
  scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3), labels = c("1/10M", "1/1M", "1/100k", "1/10k", "1/1k")) +
  annotation_logticks(sides = "l", outside = T) +  coord_cartesian(clip = "off", ylim = c(1e-7, 1e-3)) + 
  labs(x = "Family Size", y = "Error Rate", color = "", shape = "") 
  #theme(legend.position = "none")
total_ds %>% write_csv(., "~/Documents/tmp.csv")

##########
############ monomer
#######################

r1r2_and_raw = rbind(r1r2_mono, raw_mono)
colnames(r1r2_and_raw)[1:5] = c("ref_allele", "alt_allele", "total_base_errors", "total_bases_eval", "sample_id")
r1r2_and_raw = r1r2_and_raw %>% mutate(sample_id = paste(lane, sample, sep="_"))
total_mono = rbind(dupseq_mono, cds_mono, ssc_mono, r1r2_and_raw)
total_mono = total_mono %>% mutate(mut_type = paste(ref_allele, alt_allele, sep="->"))
total_mono = compute_ci(total_mono)
total_mono = total_mono %>% mutate(lower = ifelse(lower < 1e-10, 0.0, lower))
total_mono %>% filter(group == "DupSeq_SSC") %>% group_by(sample)%>% summarise(sum(total_base_errors), sum(total_bases_eval)/3)

ggplot(total_mono, aes(x = mut_type, y = mean, color = group)) +
  geom_point(size = 1, position=position_dodge2(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, position=position_dodge2(width=1)) +
  scale_y_log10(breaks = c(1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3), labels = c("1/1B","1/100M", "1/10M", "1/1M", "1/100k", "1/10k", "1/k")) +
  annotation_logticks(sides = "l", outside = T) +  coord_cartesian(clip = "off", ylim = c(1e-9, 1e-3)) + 
  labs(x = "Mononucleotide Mutation Type", y = "Error Rate", color = "", shape = "") +
  facet_wrap(~sample)
  #theme(legend.position = "none")
write_csv(total_mono, "~/Documents/tmp.csv")


########FFPE WES errorrate
########FFPE WES errorrate
#CDS
total_df = compute_ci(wes)
write_csv(total_df, "~/Documents/WES_cds.csv")


colnames(ffpemono)[1:5] = c("ref_allele", "alt_allele", "total_base_errors", "total_bases_eval", "sample_id")
ffpemono = ffpemono %>% mutate(mut_type = paste(ref_allele, alt_allele, sep=">"))
ffpemono = compute_ci(ffpemono)

  
ffpe %>% select(n_bases_eval, n_errors, erate, lane)

wesreg_df = compute_ci(wesreg)

colnames(reg_ffpemono)[1:5] = c("ref_allele", "alt_allele", "total_base_errors", "total_bases_eval", "sample_id")
reg_ffpemono = reg_ffpemono %>% mutate(mut_type = paste(ref_allele, alt_allele, sep=">"))
reg_ffpemono = compute_ci(reg_ffpemono)
reg_ffpemono = reg_ffpemono %>% mutate(lane = paste(reg_ffpemono$lane, reg_ffpemono$sample, sep='_'))
totaldf = rbind(ffpemono, reg_ffpemono)
write_csv(totaldf, "~/Documents/WES_cds_ffpe_bycontext.csv")
ffpe = totaldf %>% filter(str_detect(lane, "ffpe") | str_detect(lane, "B000102"))
normal = totaldf %>% filter(str_detect(lane, "normal") | str_detect(lane, "3BC"))

ggplot(ffpe, aes(x = mut_type, y = mean, color = lane)) +
  geom_point(size = 1, position=position_dodge2(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, position=position_dodge2(width=1)) +
  scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3), labels = c("1/10M", "1/1M", "1/100k", "1/10k", "1/k")) +
  annotation_logticks(sides = "l", outside = T) +  coord_cartesian(clip = "off", ylim = c(1e-7, 1e-3)) + 
  labs(x = "Mononucleotide Mutation Type", y = "Error Rate", color = "", shape = "") 
ggsave("~/Documents/WESErrorRateByContext.FFPE.pdf", width = 6, height = 4.5, dpi='retina', device="pdf")

ggplot(normal, aes(x = mut_type, y = mean, color = lane)) +
  geom_point(size = 1, position=position_dodge2(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, position=position_dodge2(width=1)) +
  scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3), labels = c("1/10M", "1/1M", "1/100k", "1/10k", "1/k")) +
  annotation_logticks(sides = "l", outside = T) +  coord_cartesian(clip = "off", ylim = c(1e-7, 1e-3)) + 
  labs(x = "Mononucleotide Mutation Type", y = "Error Rate", color = "", shape = "") 
ggsave("~/Documents/WESErrorRateByContext.Normal.pdf", width = 6, height = 4.5, dpi='retina', device="pdf")
