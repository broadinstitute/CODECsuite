
library(gridExtra)
library(grid)

summ <- function(df) {
  df = df %>% mutate(group = ifelse(fraglen == 0 | fraglen > 100000, "intermol", "cds"))
  summa = df %>% group_by(group) %>% summarise(R1_q0_total_error = sum(R1_q0_nerror), R1_q0_total_bases = sum(R1_q0_efflen), 
                                      R2_q0_total_error = sum(R2_q0_nerror), R2_q0_total_bases = sum(R2_q0_efflen), 
                                      R1_q30_total_error = sum(R1_q30_nerror), R1_q30_total_bases = sum(R1_q30_efflen), 
                                      R2_q30_total_error = sum(R2_q30_nerror), R2_q30_total_bases = sum(R2_q30_efflen), 
                                      R1_q0_erate = R1_q0_total_error / R1_q0_total_bases,
                                      R2_q0_erate = R2_q0_total_error / R2_q0_total_bases,
                                      q0_erate = (R1_q0_total_error + R2_q0_total_error)/(R1_q0_total_bases + R2_q0_total_bases), 
                                      R1_q30_erate = R1_q30_total_error / R1_q30_total_bases,
                                      R2_q30_erate = R2_q30_total_error / R2_q30_total_bases, 
                                      q30_erate = (R1_q30_total_error + R2_q30_total_error)/(R1_q30_total_bases + R2_q30_total_bases))
  summa
}


novaseq_custom_primerV2 = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_capture/adapter_V2_first_capture/wrk/raw_error/adapterv2.index2.raw.readlevel.txt"
novaseq_custom_primerV2 = reformat_on_mac(novaseq_custom_primerV2)
df= fread(novaseq_custom_primerV2)
head(df)
df = df %>% mutate(R1_Error_rate = R1_q0_nerror / R1_q0_efflen, R2_Error_rate = R2_q0_nerror / R2_q0_efflen, 
                   R1_Q30_Error_rate = R1_q30_nerror / R1_q30_efflen, R2_Q30_Error_rate = R2_q30_nerror / R2_q30_efflen)
head(df)

idx = sample(1:nrow(df), 1000000)


novaseq_custom_primerV1 = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/novaseq_NA19238_custom_v1/wrk/raw_error/pre30_1.raw.readlevel.txt"
novaseq_custom_primerV1 = reformat_on_mac(novaseq_custom_primerV1)
df2= fread(novaseq_custom_primerV1)
df2 = df2 %>% mutate(R1_Error_rate = R1_q0_nerror / R1_q0_efflen, R2_Error_rate = R2_q0_nerror / R2_q0_efflen, 
                   R1_Q30_Error_rate = R1_q30_nerror / R1_q30_efflen, R2_Q30_Error_rate = R2_q30_nerror / R2_q30_efflen)
idx2 = sample(1:nrow(df2), 1000000)


plot_heatmap <- function(df, col, from = 150, to = 950) {
  df = df %>% filter(fraglen > from & fraglen < to & eval(parse(text=col)) < 0.2) 
  p = ggplot(data=df, aes(x=fraglen, y=eval(parse(text=col)))) + 
    #geom_density2d()
    geom_bin2d(binwidth= c(20, 0.01)) + 
    #scale_fill_gradient(low="lightblue1",high="darkblue",trans="log10") +
    scale_fill_distiller(palette=4, direction=1, trans="log10", limits=c(1, max=1e6)) +
    geom_smooth(color='red') +
    ylab(col) + 
    scale_y_continuous(breaks=seq(0, 0.2, 0.01), labels = seq(0, 0.2, 0.01)) 
  return(p) 
}

p1=plot_heatmap(df2[idx2, ], "R1_Error_rate")
p2=plot_heatmap(df2[idx2, ], "R2_Error_rate")
p3=plot_heatmap(df[idx, ], "R1_Error_rate")
p4=plot_heatmap(df[idx, ], "R2_Error_rate")
grid.arrange(p1, p3, p2, p4, nrow = 1)

hiseq_custom_primerV1 = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/NA19238_custom_v1/wrk/accu_out/cds/pre30C_1_HCVH2BCX3.raw.readlevel.txt"
hiseq_custom_primerV1 = reformat_on_mac(hiseq_custom_primerV1)
hiseq_V1.5_df = fread(hiseq_custom_primerV1, nrow=1000000)

hiseq_V1.5_df = hiseq_V1.5_df %>% mutate(R1_Error_rate = R1_q0_nerror / R1_q0_efflen, R2_Error_rate = R2_q0_nerror / R2_q0_efflen, 
                   R1_Q30_Error_rate = R1_q30_nerror / R1_q30_efflen, R2_Q30_Error_rate = R2_q30_nerror / R2_q30_efflen)
hiseq_V1.5_df = hiseq_V1.5_df %>% mutate(group = ifelse(fraglen == 0 | fraglen > 100000, "intermol", "cds"))
hiseq_V1.5_df %>% group_by(group) %>% summarise(R1_q0_total_error = sum(R1_q0_nerror), R1_q0_total_bases = sum(R1_q0_efflen), R1_q0_erate = R1_q0_total_error / R1_q0_total_bases,
                                      R2_q0_total_error = sum(R2_q0_nerror), R2_q0_total_bases = sum(R2_q0_efflen), R2_q0_erate = R2_q0_total_error / R2_q0_total_bases)
colnames(hiseq_V1.5_df)


###############
##########################################
hiseq_v1_cap = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_capture/cds_first_caputre/family_consensus_wrk/accu/cds/CDS_HD82_DC_1_HCY3HBCX3.raw.readlevel.txt"
hiseq_v1_cap = reformat_on_mac(hiseq_v1_cap)
hiseq_v1_cap_df = fread(hiseq_v1_cap) 
print(summ(hiseq_v1_cap_df), width=Inf)


p1 = plot_heatmap(df3[idx3, ], "R1_Error_rate")
p2 = plot_heatmap(df3[idx3, ], "R2_Error_rate")
p3 = plot_heatmap(df3[idx3, ], "R1_Q30_Error_rate")
p4 = plot_heatmap(df3[idx3, ], "R2_Q30_Error_rate")
grid.arrange(p1, p2, p3, p4, nrow=1)


###########

hiseq_dupseq = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_capture/cds_first_caputre/wrk_dup/accu/Dup_HD82_DC_1_HCY3HBCX3.fastp.readlevel.txt"
hiseq_dupseq = reformat_on_mac(hiseq_dupseq)
hiseq_dupseq = fread(hiseq_dupseq)
print(summ(hiseq_dupseq), width=Inf)




############
#inputdir = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/NA19238_custom_v1/wrk/accu_out/cds")
#inputdir = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/novaseq_NA19238_custom_v1/wrk/accu_out/cds")
inputdir = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/NA19238_with_blocker/wrk_new_pip_noumi/accu_out/cds")
dfall = data.frame()
for (f in list.files(inputdir, "error_metrics.txt")) {
  sampleid = unlist(strsplit(f, "\\."))[1]
  df = fread(file.path(inputdir, f)) 
  #colnames(df) = c(names(df)[2:27], "V1")
  df$sample = sampleid
  dfall = rbind(dfall, df)
}
dfall
dfall %>% select(sample,q0_erate, q0_R1_erate, q0_R2_erate, q30_erate, q30_R1_erate, q30_R2_erate) %>% write_csv(., path = "July27th.3lib_rawerrorrate.csv")
#library(binom)
#dfprint = dfall %>% rowwise() %>% mutate(lower = unlist(binom.confint(q0_n_errors, q0_n_bases_eval, methods = 'wilson'))[5], upper = unlist(binom.confint(q0_n_errors, q0_n_bases_eval, methods = 'wilson'))[6])
#print(dfprint, width = Inf)

#######
#inputdir = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/novaseq_NA19238_custom_v1/wrk/metrics_out")
inputdir = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/NA19238_with_blocker/wrk_new_pip_noumi/metrics_out")

cycle_df <-  function(inputdir) {
  dfall = data.frame()
  for (f in list.files(inputdir, "error_metrics.txt")) {
    sampleid = unlist(strsplit(f, "\\."))[1]
    df = fread(file.path(inputdir, f), skip=2) 
    #colnames(df) = c(names(df)[2:27], "V1")
    df$sample = sampleid
    df$type = substr(sampleid, 1, 3)
    dfall = rbind(dfall, df)
  }
  dfall
}
inputdir = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_capture/cds_first_caputre/family_consensus_wrk/accu/cds/")
input = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_capture/cds_first_caputre/wrk_dup/accu")
cds = cycle_df(inputdir)
dupseq = cycle_df(input)
dfall = rbind(cds, dupseq)
unique(cds$sample)
plotdf = dupseq %>% filter(sample %in% c("Dup_HD82_DC_1_HCY3HBCX3", "Dup_HD82_DC_2_HCY3HBCX3")) %>% filter(cycle < 244) %>% group_by(type, cycle) %>% summarise(R1_q0_erate = sum(R1_q0_error) / sum(R1_q0_cov), R2_q0_erate = sum(R1_q0_error) / sum(R2_q0_cov) )
ggplot(data=plotdf, aes(x=cycle + 7, y=R1_q0_erate, color=type)) + geom_point() 
View(dupseq)

v1.5dir = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/NA19238_custom_v1/newpip_wrk/accu_out/cds")
df = cycle_df(v1.5dir)
plotdf = df %>% filter(cycle < 229) %>% group_by(cycle) %>% summarise(R1_q0_erate = sum(R1_q0_error) / sum(R1_q0_cov), R2_q0_erate = sum(R1_q0_error) / sum(R2_q0_cov) )
View(plotdf)
ggplot(data=plotdf, aes(x=cycle + 7, y=R1_q0_erate, color=type)) + geom_point() + scale_y_log10()

############
########Capture with Blocker
inputdir = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_capture/cds_capture_with_blocker/new_pip_wrk/accu/cds")
cds = cycle_df(inputdir)
plotdf = cds %>% group_by(cycle) %>% summarise(R1_q0_erate = sum(R1_q0_error) / sum(R1_q0_cov), R2_q0_erate = sum(R1_q0_error) / sum(R2_q0_cov) )
View(plotdf)
ggplot(data=plotdf, aes(x=cycle + 7, y=R1_q0_erate, color=type)) + geom_point() 
+ scale_y_log10()


duplex = "/xchip/bloodbiopsy/ruolin/link_duplex/cds_capture/adapter_V2_first_capture/wrk/detect"

detect = prepare_detect(duplex, "HD82_Sep25")
detect %>% group_by(id, group) %>% summarise(depth_pass = sum(depth_pass)) %>% ggplot(data=., aes(y=depth_pass, fill=group)) + geom_boxplot() + theme_light(22) + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) 

cap1 = prepare_detect(highconf_detect, "HD82_Jun22")



readlevel = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_capture/cds_first_caputre/family_consensus_wrk/accu/cds/CDS_05246_315_DC_2_HCY3HBCX3.raw.readlevel.txt")
pat_rl_df = fread(readlevel)
head(pat_rl_df, n=30)
plotdf = pat_rl_df %>% mutate (intermol = ifelse(fraglen == 0 | fraglen > 5000, TRUE, FALSE)) %>% group_by(intermol) %>% summarise(R1_q0_er = sum(R1_q0_nerror) / sum(R1_q0_efflen), 
                                                                                                                                   R1_q30_er = sum(R1_q30_nerror)/sum(R1_q30_efflen), 
                                                                                                                                   R2_q0_er = sum(R2_q0_nerror)/sum(R2_q0_efflen), 
                                                                                                                                   R2_q30_er = sum(R2_q30_nerror)/sum(R2_q30_efflen))

v1dir = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_capture/cds_first_caputre/family_consensus_wrk/accu/cds/")
v2dir = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/pcr_condition_v2_hiseq_250x_NA12878/wrk/accu_out/consensus")
df = error_df(v2dir)
df$type = "CDSV2"
plotdf = rbind(df, dfv1.5, regular_df)
ggplot(data=plotdf, aes(x=cycle, y=R1_q30_erate, color=type)) + geom_point() + scale_y_log10(limits=c(1e-4,1)) + xlim(0,250)  + theme_light(base_size=20) 

regular = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/regular_wgs/hiseq2500RR_HG005/newpip_wrk/accu_out/HG005.readlevel.txt")
df2= fread(regular, nrow=1000000)
df2 = df2 %>% mutate(R1_Error_rate = R1_q0_nerror / R1_q0_efflen, R2_Error_rate = R2_q0_nerror / R2_q0_efflen, 
                   R1_Q30_Error_rate = R1_q30_nerror / R1_q30_efflen, R2_Q30_Error_rate = R2_q30_nerror / R2_q30_efflen)

plot_heatmap(df2, "R2_Q30_Error_rate")

error_df <- function(inputdir){
  dfall = data.frame()
  for (f in list.files(inputdir, ".error_metrics.txt")) {
    name = unlist(strsplit(f, "\\."))
    sample = name[2]
    lane = name[1]
    df = fread(file.path(inputdir, f), nrows =2) 
    #colnames(df) = c(names(df)[2:27], "V1")
    df$sample = sample
    df$lane = lane
    dfall = rbind(dfall, df)
  }
  dfall
}
v2dir_nov3 = reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/pcr_condition_v2_hiseq_250x_NA12878/wrk/accu_out/consensus")
df = error_df(v2dir_nov3)
df
df %>% select(sample, lane, q30_erate, num_pairs, q30_n_bases_eval) %>% write.csv(., "v2.consensus.Nov01.csv")

substr("abcd", 2, 3)
str_length("abcd")
