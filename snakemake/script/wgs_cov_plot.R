library(ggplot2)
library(gtools)
plot_theme <- theme_bw() + theme(text = element_text(size=10, color = "black"),
                                 axis.text = element_text(size = 10, color = "black"),
                                 axis.text.y = element_text(margin = margin(r=10)),
                                 axis.title = element_text(size = 10, color = "black"),
                                 title = element_text(size = 10, color = "black"),
                                 legend.text = element_text(size = 10, color = "black"),
                                 legend.title = element_text(size = 10, color = "black"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())
theme_set(plot_theme)

#### Coverage plot
prepare_df = function(file, group) {
  df = fread(file, skip=9)%>% select(-unfiltered_baseq_count)
  df$frac_eql_larger = rev(cumsum(as.numeric(rev(df$high_quality_coverage_count)))) / sum(df$high_quality_coverage_count)
  df$frac = df$high_quality_coverage_count / sum(as.numeric(df$high_quality_coverage_count)) 
  df$group = group
  df
}

#dupseq_r1r2 = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/dupseq/wgs_NA12989_dec8/wrk/metrics/NA12878_merged.r1r2consensus.wgs_metrics.txt")
#dup_r1r2_df = prepare_df(dupseq_r1r2, "DupSeq_R1R2")

dup_dsc = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/dupseq/wgs_NA12989_dec8/wrk/metrics/NA12878_merged.duplex_consensus.wgs_metrics.txt")
dup_dsc_df = prepare_df(dup_dsc, "DupSeq_DSC")

cds = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/pcr_condition_v2_hiseq_250x_NA12878/wrk/metrics/CDS_V2_merged.cds_consensus.mol_consensus.wgs_metrics.txt")
cds_df = prepare_df(cds, "CDS")
df = rbind(dup_dsc_df, cds_df)

cds_filter = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/pcr_condition_v2_hiseq_250x_NA12878/vcwrk/output/CDS_V2.filtered.markdup.trimend.wgs_metrics.txt")
cds_filter_df = prepare_df(cds_filter, "CDS_filter")
df = rbind(cds_df, cds_filter_df)


ggplot(data=df) +
  geom_col(aes(x=coverage, y=frac, fill = group), position=position_dodge(), alpha=0.6) + 
  geom_line(aes(x=coverage, y=frac_eql_larger, color=group))+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 15), limits=c(-1,15)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) + 
  xlab("Coverage") + ylab("Fraction of unique depth >= or == ") +
  #geom_vline(xintercept = 3.96, color="red", linetype="dashed", size=1.5) +
  #geom_vline(xintercept = 0.025, color="blue", linetype="dashed", size=1.5) +
  theme(legend.position = c(0.8, 0.8))
ggsave("~/Documents/WGScov_cds_vs_filter.pdf", width = 6, height = 4.5, dpi='retina')


#### Error rate by context
get_monomer_erate <- function(input, group) {
  dup_raw = fread(reformat_on_mac(input))
  colnames(dup_raw)[1:3] = c("group", 'ref', 'alt')
  dup_raw$group = group
  dup_raw$total_bases_eval = as.numeric(dup_raw$total_bases_eval)
  dup_raw = dup_raw %>% mutate(total_bases_eval = floor(total_bases_eval/2), total_base_errors = floor(total_base_errors/2))
  dup_raw
}
cds_hiseq = monomer_fromc('/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/pcr_condition_v2_hiseq_250x_NA12878/wrk/accu_out/consensus', "CDS hiseq", 30)
cds_nova = monomer_fromc('/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/pcr_condition_v2_novaseq_250x_NA12878/wrk/accu_out/consensus', "CDS novaseq", 30)
cds = rbind(cds_hiseq, cds_nova)
cds = cds %>% group_by(group, ref, alt) %>% summarise(n = sum(n), den = sum(den))
cds$er = cds$n/cds$den
cds$sample="error_metrics"
cds$lane = "NA12878"
dup_dsc = monomer_fromc('/xchip/bloodbiopsy/ruolin/link_duplex/dupseq/wgs_NA12989_dec8/wrk/accu/duplex', 'DupSeq', 30)
dup_r1r2 = monomer_fromc('/xchip/bloodbiopsy/ruolin/link_duplex/dupseq/wgs_NA12989_dec8/wrk/accu/r1r2', 'R1R2', 30)
raw_q30 = monomer_fromc('/xchip/bloodbiopsy/ruolin/link_duplex/dupseq/wgs_NA12989_dec8/wrk/accu/raw', "DupSeq_raw_q30", 30)
#raw_q0 = get_monomer_erate('/xchip/bloodbiopsy/ruolin/link_duplex/dupseq/wgs_NA12989_dec8/wrk/accu/raw_q0/monomer_erate.txt', 'Raw_q0')
df_err = rbind(cds %>% ungroup(), dup_dsc, dup_r1r2, raw_q30) 
colnames(df_err)[4:5] = c("total_base_errors", "total_bases_eval")
df_err = compute_ci(df_err)
df_err = df_err %>% mutate(mut_type = paste(ref, alt, sep="->"))

#df_err %>% filter(group=='CDS') %>% select(mean) / (df_err %>% filter(group=='DupSeq') %>% select(mean))
df_err %>% write.csv(., file='~/Documents/tmp.csv', quote = FALSE)
df_err %>% group_by(group) %>% summarise(sum(total_base_errors), sum(total_bases_eval)/3)

ggplot(df_err, aes(x = mut_type, y = mean, color = group)) +
  geom_point(size = 1, position=position_dodge2(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, position=position_dodge2(width=1)) +
  scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3), labels = c("1/10M", "1/1M", "1/100k", "1/10k", "1/k")) +
  annotation_logticks(sides = "l", outside = T) +  coord_cartesian(clip = "off", ylim = c(1e-7, 1.5e-3)) + 
  labs(x = "Mononucleotide Mutation Type", y = "SNV Error Rate", color = "", shape = "") +
  ggtitle("NA12878 high confident regions, indel removed")
  #theme(legend.position = "none")
ggsave("~/Documents/WgsErrorRateByContext.pdf", width = 6, height = 4.5, dpi='retina')

#### Overall Error rate plot
#############
cds_hiseq = to_miredas_errordf(prepare_error("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/pcr_condition_v2_hiseq_250x_NA12878/wrk/accu_out/consensus", "CDS\nhiseq"))
cds_novaseq = to_miredas_errordf(prepare_error("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/pcr_condition_v2_novaseq_250x_NA12878/wrk/accu_out/consensus", "CDS\nnovaseq"))
cds = rbind(cds_hiseq, cds_novaseq)
cds = cds %>% group_by(group) %>% summarise(total_bases_eval = sum(total_bases_eval), total_base_errors = sum(total_base_errors)) %>% ungroup()
cds$lane = "NA12878"
dupseq = to_miredas_errordf(prepare_error("/xchip/bloodbiopsy/ruolin/link_duplex/dupseq/wgs_NA12989_dec8/wrk/accu/duplex", "DupSeq")) %>% select(-c(sample, sample_id))
dupseq_r1r2 = to_miredas_errordf(prepare_error("/xchip/bloodbiopsy/ruolin/link_duplex/dupseq/wgs_NA12989_dec8/wrk/accu/r1r2", "R1R2"))%>% select(-c(sample, sample_id))
dupseq_raw = to_miredas_errordf(prepare_error("/xchip/bloodbiopsy/ruolin/link_duplex/dupseq/wgs_NA12989_dec8/wrk/accu/raw", "Raw"))%>% select(-c(sample, sample_id))
total_df = rbind(cds, dupseq, dupseq_r1r2, dupseq_raw)
total_df = compute_ci(total_df)
total_df = total_df %>% separate(col = group, sep="_q", into = c("group", "qcut"))

ggplot(total_df %>% filter(qcut == 30 ), aes(x = group, y = mean, color=qcut)) +
  geom_point(size = 1, position=position_dodge2(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, position=position_dodge2(width=1)) +
  scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3), labels = c("1/10M", "1/1M", "1/100k", "1/10k", "1/1k")) +
  annotation_logticks(sides = "l", outside = T) +  coord_cartesian(clip = "off", ylim = c(1e-7, 1e-3)) + 
  labs(x = "", y = "Error Rate", color = "q-score cutoff", shape = "") 
  #theme(legend.title = element_text("q-score cutoff"))
ggsave("~/Documents/WGSTotalErrorRate.pdf", width = 6, height = 4.5, dpi='retina')
write_csv(total_df %>% arrange(lane), "~/Documents/tmp.csv")

