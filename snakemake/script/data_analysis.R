##########
############## CDS analysis
###########################

####
############## Byproducts
###########################
#inputdir = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_cds_highconf/hs_metrics_out"
#inputdir = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_newsplit_lowconf/hs_metrics_out"
blocker = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_capture_with_blocker/wrk_umi/hs_metrics_out"
noblocker = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_newsplit_highconf/hs_metrics_out"
wgsblocker = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/NA19238_with_blocker/wrk_nofastp2/metrics_out"
inputdir = wgsblocker
alldf = data.frame()
for (f in list.files(reformat_on_mac(inputdir), pattern = "*agg_metrics.txt")) {
  df = fread(reformat_on_mac(file.path(inputdir, f)))
  alldf = rbind(alldf, df)
}
barplotdf = alldf%>%filter %>% select(sample_id, pct_adp_dimer, pct_intermol, pct_cds_highconf, pct_cds_lowconf, pct_single_highconf, pct_single_lowconf, pct_double_ligation) %>% pivot_longer(cols = -c(sample_id))
ggplot(barplotdf, aes(x=sample_id, y= value, fill=name)) + geom_bar(stat="identity", position="stack") + 
  scale_y_continuous(breaks=seq(0,1, 0.1), labels=seq(0, 1, 0.1))+ scale_x_discrete(labels = c("median adpt", "long adpt")) +
  theme_classic(base_size = 32) 

inputdir = noblocker
alldf = data.frame()
for (f in list.files(reformat_on_mac(inputdir), pattern = "^CDS.*_1_.*agg_metrics.txt")) {
  df = fread(reformat_on_mac(file.path(inputdir, f)))
  alldf = rbind(alldf, df)
}
noblocker_alldf = alldf
noblocker_alldf

inputdir = wgsblocker
alldf = data.frame()
f20or (f in list.files(reformat_on_mac(inputdir), pattern = "*agg_metrics.txt")) {
  df = fread(reformat_on_mac(file.path(inputdir, f)))
  alldf = rbind(alldf, df)
}
wgsblocker_alldf = alldf

blocker_alldf$group="cpt_blocker"
wgsblocker_alldf$group="wgs_blocker"
noblocker_alldf$group="cpt"
df = rbind(blocker_alldf, wgsblocker_alldf, noblocker_alldf)
df
ggplot(df, aes(x=sample_id, y=pct_intermol, fill=group)) + geom_bar(stat="identity") + ylim(0,1) + 
  scale_x_discrete(labels=rep("", 6)) + theme_classic(base_size = 20)  

write.csv(alldf, "cds_with_blocker_lane1.summary.csv")
#write_csv(alldf, path = "cds_first_caputer_byprodcut.csv")
barplotdf = alldf%>%filter %>% select(sample_id, pct_adp_dimer, pct_intermol, pct_high_conf, pct_low_conf, pct_single_insert, pct_double_ligation) %>% pivot_longer(cols = -c(sample_id))
ggplot(barplotdf, aes(x=sample_id, y= value, fill=name)) + geom_bar(stat="identity", position="stack") + scale_y_continuous(breaks=seq(0,1, 0.1), labels=seq(0, 1, 0.1))+
  theme_classic(base_size = 22) + 
  #scale_fill_discrete(name = "", labels = c("Adpter Dimer", "Intermolecular", "Overalp Pair")) +
  theme(axis.ticks.x=element_blank(),axis.text.x=element_text(angle=90)) + ylab("Fraction") + xlab("Sample") + coord_cartesian(ylim = c(0,1)) + theme(legend.position = "top")

dbinom(0, prob=0.1, size=10)
inputdir = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/NA19238_with_blocker/wrk_nofastp/metrics_out"
alldf = data.frame()
for (f in list.files(reformat_on_mac(inputdir), pattern = "*insert_size_metrics.txt")) {
  sample = unlist(str_split(f, "_"))[1]
  df = fread(reformat_on_mac(file.path(inputdir, f)), skip=10)
  df$sample = sample
  alldf = rbind(alldf, df)
}
ggplot(data = alldf, aes(x=insert_size, y = All_Reads.fr_count, color = sample)) + geom_line() + xlim(0, 200) + theme_bw(base_size = 28)

pbinom(0, size=13, prob=0.01)
dhyper(1, m = 1, n = 1e6-1, k=1e5)

?phyper
###
#### insert_size hist
#########



####
#############pertarget
#####################
highconf_hsmetrics = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_cds_highconf/detect_fingerprints") 
dup_hsmetrics = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_dup/detect_fingerprints") 
comb_hsmetrics = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_cds_both_highlowconf/detect_fingerprints") 
prepare_efficiency <- function(folder, group) {
  dfall = data.frame()
  i = 1
  #for (f in list.files(folder, pattern="duplex.per_target_cov.txt$")) {
  for (f in list.files(folder, pattern="detect_fingerprint.txt$")) {
    df = fread(file.path(folder, f))
    df$sample = unlist(str_split(f, "_"))[2]
    df$sid = i %% 2
    dfall = rbind(dfall, df)
    i = i + 1
  }
  dfall$group = group 
  dfall
}
a = prepare_efficiency(highconf_hsmetrics, "HighConf")
b = prepare_efficiency(dup_hsmetrics, "Dup")
c = prepare_efficiency(comb_hsmetrics, "Combined")
df = rbind(a,b, c)
df$sample_id = paste(df$sample, df$sid, sep="_")

stat_box_data <- function(y) {
  return( 
    data.frame(
      y = mean(y),
      label = paste('mean =', round(mean(10^y), 1), '\n')
    )
  )
}
ggplot(data = df, aes(x= sample_id, y=depth, fill = group)) + geom_boxplot() + 
  theme_classic(base_size = 28) + #scale_y_continuous(breaks=(seq(0, 700,100))) + 
  ylab("prefilter depth") +
  stat_summary(fun.data = stat_box_data, geom = "text",hjust = 0.5, vjust = 0.9, size=6,  position = position_dodge(width = 0.75))+
  scale_y_log10(breaks = c(10, 100, 1e3, 1e4), labels = c("10", "100", "1k", "10k"))+ annotation_logticks()


#######
################# UMI families
############################
highconf_fsize_hist = fread(reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_cds_highconf/hs_metrics_out/CDS_HD82_DC_1_HCY3HBCX3.family_sizes.txt"))
lowconf_fsize_hist = fread(reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_cds_lowconf/hs_metrics_out/CDS_HD82_DC_1_HCY3HBCX3.family_sizes.txt"))
dup_fsize_hist = fread(reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_dup/hs_metrics_out/Dup_HD82_DC_1_HCY3HBCX3.family_sizes.txt"))

#highconf_fsize_hist = fread(reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_cds_highconf/hs_metrics_out/CDS_05246_315_DC_1_HCY3HBCX3.family_sizes.txt"))
#lowconf_fsize_hist = fread(reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_cds_lowconf/hs_metrics_out/CDS_05246_315_DC_1_HCY3HBCX3.family_sizes.txt"))
#dup_fsize_hist = fread(reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_dup/hs_metrics_out/Dup_05246_315_DC_1_HCY3HBCX3.family_sizes.txt"))
colnames(highconf_fsize_hist)[10] = "HighConf"
colnames(lowconf_fsize_hist)[10] = "LowConf"
colnames(dup_fsize_hist)[10] = "DUP"
#ggplot(highconf_fsize_hist, aes(x = family_size, y=HighConf)) + geom_point() + xlim(0, 1000) + ylim(0,0.8) + ylab("fraction_gt_or_eq_family_size") + theme_classic(base_size = 22)

library(ggplot2)
joined = highconf_fsize_hist %>% full_join(lowconf_fsize_hist, by = "family_size")
joined = joined %>% full_join(dup_fsize_hist, by= "family_size")
head(joined)
pivot_joined = joined %>% select(family_size, HighConf, LowConf, DUP) %>% pivot_longer(cols= - family_size, names_to = "Type")
ggplot(pivot_joined %>% filter(family_size <= 500), aes(x = family_size, y = value, color=Type)) + geom_point() + xlim(0, 500) + ylim(0,0.8) + ylab("fraction_gt_or_eq_family_size") + theme_classic(base_size = 22)


#######
##########Detection rate
######################
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
colnames(detectdf)
unique(detectdf$sample)
plotdf = detectdf %>% group_by(sample, group, id) %>% summarise(saf = sum(n_alt_pass)/sum(depth_pass)) 
#plotdf = detectdf %>% group_by(sample, group, id) %>% summarise(taf = sum(n_alt_pass)) 
hddf = plotdf %>% filter(sample == "HD82")
hddf[3,] = list("HD82", "Duplex", "1", 0)

ggplot(data =  plotdf %>% filter(sample == "05246"), aes(x = group, y = saf, fill=sample)) + geom_boxplot() + 
  stat_summary(fun.data = stat_box_data, size=8, geom = "text",hjust = 0.5, vjust = 0.9,  position = position_dodge(width = 0.75)) + 
  theme_classic(base_size = 28) + ylim(0,0.5) + 
   theme(axis.title.x=element_blank()) +
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank()) + 
  theme(legend.position = "none") +
  ylab("Somatic VAF") +   theme(legend.title=element_blank()) 

#ggplot(data = rbind(a,b) %>% filter(sample == "HD82"), aes(as.factor(sid) , y = n_alt_pass,  fill= group)) + geom_bar(stat = "identity", position = "dodge2") + 
#  xlab("replicate") + theme_classic(base_size = 22) + scale_y_continuous(breaks=c(1))

############
##################
###Error rates

highconf_error = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_newlinker_cds_highconf/miredas"
dupssc_error = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_dup/miredas/mol"
dupdsc_error = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_dup/miredas/duplex"
singleinsert_error = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_newlinker_singleinsert/miredas"
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

a = prepare_error(highconf_error, "CDS")
b = prepare_error(dupssc_error, "SSC")
d = prepare_error(dupdsc_error, "DSC", adjust=T)
e = prepare_error(singleinsert_error, "SingleInsert")

df = rbind(a, b, d, e)
df$sample_id = paste(df$sample, df$sid, sep="_")
plotdf = df %>% group_by(sample, group) %>% summarise(error_rate = sum(total_base_errors)/sum(total_bases_eval))
ggplot(data= plotdf, aes(y = error_rate, x = sample, color=group)) + geom_point(size=4) +
  scale_y_log10(limits = c(1e-7, 1e-3), breaks = c(1e-3,1e-4,1e-5,1e-6,1e-7), labels = c("1/k", "1/10k", "1/100k", "1/1M", "1/10M")) +
  #scale_x_discrete(labels = df$sample) + 
  xlab("")  + ylab("Error Rate") + 
  annotation_logticks(sides = "l", short = unit(0.2, "cm"), mid = unit(0.35, "cm"), long = unit(0.5, "cm")) + theme_classic(base_size = 28) 

###### By family size
pat_fsize_error = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_newlinker_cds_highconf/fsize_erate/miredas/CDS_05246_315_DC"
hd_fsize_error = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/wrk_newlinker_cds_highconf/fsize_erate/miredas/CDS_HD82_DC"
pat = prepare_error(pat_fsize_error, "05246_315", group_size = 3)
hd = prepare_error(hd_fsize_error, "HD_82", group_size = 3)
combined = rbind(pat, hd)
combined = combined %>% mutate(type = ifelse(str_detect(sample_id, "SSC"), "SSC", "Duplex"))
combined$sid = combined$sid + 1
#combined = combined %>% mutate(error_rate = ifelse(total_base_errors == 0, 1/total_bases_eval, errors_per_base_sequenced), 
#                    err_type = ifelse(total_base_errors == 0, "theory", "acutal"))
ggplot(data= combined %>% filter(total_base_errors > 0 ), aes(y = errors_per_base_sequenced, x = sid, shape=type, color=factor(group))) + geom_point(size=4, position = position_dodge(width=0.3)) + 
  scale_y_log10(limits = c(1e-7, 1e-3), breaks = c(1e-3,1e-4,1e-5,1e-6,1e-7), labels = c("1/k", "1/10k", "1/100k", "1/1M", "1/10M")) +
  scale_x_continuous(labels=c("1", "2", "3"), breaks=1:3) + xlab("Family Size") + ylab("Error Rate") + labs(color = "patient") + 
  annotation_logticks(sides = "l", short = unit(0.2, "cm"), mid = unit(0.35, "cm"), long = unit(0.5, "cm")) + theme_classic(base_size = 28) 


###########
####### 
a = fread(reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/data/allinput.tsv"))
head(a)
a = a %>% mutate(error_maf_file = ifelse(patient_id == "HD82", "/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/data/HD82.maf", "/xchip/bloodbiopsy/rhoades/analysis/05246/ERp_MBC_probe_design/design_probes/BRCA-05246_CCPM_0300315_probes.maf"))
write_tsv(a, path = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_first_caputre/data/allinput.tsv"))


a = read_tsv(reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/NA19238_with_blocker/inputall.tsv"))
a$patient_id = "NA19238"
a = a %>% select(colnames(a)[1:8])
colnames(a)[8] = "interval_list"
write_tsv(a, reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/NA19238_with_blocker/inputall.tsv"))
