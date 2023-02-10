source("paper/utils.R")
### Overall error rate
plot_theme <- theme_bw() + theme(text = element_text(size=10, color = "black"),
                                 axis.text = element_text(size = 16, color = "black"),
                                 axis.text.y = element_text(margin = margin(r=10)),
                                 axis.title = element_text(size = 18, color = "black"),
                                 title = element_text(size = 10, color = "black"),
                                 legend.text = element_text(size = 14, color = "black"),
                                 legend.title = element_text(size = 10, color = "black"),

                                                                  legend.background = element_blank(),
                                 legend.box.background = element_blank(),
                                 panel.grid.major = element_blank(), panel.grid.minor = element_blank())
theme_set(plot_theme)

### Capture data Fig2a,b
###### CDS PanCancerdata Error rate, Fig2a
header = c("sample_id", "total_base_errors", "total_bases_eval", "sample", "lane", "group")
dupseqdf = prepare_error("figure2/a", pattern="dsc_mutation_metrics.txt", "DupSeq", sep="_sc_idt") %>% select(all_of(header))
cdsdf = prepare_error("figure2/a", pattern="cds_consensus.mol_consensus_mutation_metrics.txt", "CODEC", sep="_sc_idt") %>% select(header)
ssc = prepare_error("figure2/a", pattern="Dup_pan.*mutation_metrics.txt", "ssc", sep="_ctDNA.")
ssc = ssc %>% mutate(lane = str_replace(lane, "Dup_pan.", ""))
fig2adf = rbind(dupseqdf, cdsdf, ssc)
fig2adf = compute_ci(fig2adf)

ggplot(fig2adf , aes(x = group, y = mean, color = lane)) +
  geom_point(size = 1, position=position_dodge2(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, position=position_dodge2(width=1)) +
  scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3), labels= trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l", outside = T) +  coord_cartesian(clip = "off", ylim = c(1e-7, 1e-3)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_color_manual(values = c("black", "grey")) +
  labs(x = "", y = "Error Rate", color = "", shape = "") 

pt = prepare_error("figure2/b/05246_315_ctDNA", "patient")
hd = prepare_error("figure2/b/HD_82_ctDNA", "hd")
fig2bdf = rbind(pt, hd)
colnames(fig2bdf)[8] = "n_snv"
fig2bdf = fig2bdf %>% mutate(lane = str_replace(lane, "fsize_", ""))
fig2bdf = compute_ci(fig2bdf)

ggplot(fig2bdf , aes(x = lane, y = mean, color = group)) +
  geom_point(size = 1, position=position_dodge2(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, position=position_dodge2(width=1)) +
  scale_y_log10(breaks = c(1e-7, 1e-6, 1e-5, 1e-4, 1e-3), labels= trans_format("log10", math_format(10^.x))) +
  annotation_logticks(sides = "l", outside = T) +  coord_cartesian(clip = "off", ylim = c(5e-8, 1e-3)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_color_manual(values = c("grey", "black")) +
  labs(x = "", y = "Error Rate", color = "", shape = "") 
  #theme(legend.position = "none")

## Fig 2C, 
#check out duplex_yield.ipynb


## Fig 2D, Cost vs accuracy plot
cds_wgs = prepare_error("figure2/d/codec/", "CODEC",pattern = "NA12878.*mutation_metrics.txt")
cds_wgs = cds_wgs %>% select(sample, lane, group, n_bases_eval, n_snv) %>% mutate(sample_id = "NA12878.mutation_metrics")
colnames(cds_wgs)[4:5] = c("total_bases_eval", "total_base_errors")
cds_wgs = cds_wgs %>% group_by(group) %>% summarise(sample = "mutation_metrics", lane="NA12878", sample_id="NA12878.mutation_metrics", total_bases_eval = sum(total_bases_eval), total_base_errors=sum(total_base_errors))
dup_wgs = prepare_error("figure2/d/dupseq", "DupSeq")
dup_wgs = dup_wgs %>% select(sample, lane, group, n_bases_eval, n_snv) %>% mutate(sample_id = "NA12878.mutation_metrics")
colnames(dup_wgs)[4:5] = c("total_bases_eval", "total_base_errors")
std_wgs = to_miredas_errordf(prepare_error("figure2/d/standard", "Standard NGS"))

region_size = 2575064465
novaseq_factor = 160/17
q30_factor = 0.85
wgs_df = rbind(cds_wgs, dup_wgs, std_wgs)
wgs_df = compute_ci(wgs_df, 0.99)
wgs_df$cost_per_mb = 8600/novaseq_factor/(wgs_df$total_bases_eval/1e6 * 3.1e9/region_size)
wgs_df$group = str_replace(wgs_df$group, "_q30", "")
wgs_df[wgs_df$group == "Standard NGS",]$cost_per_mb = 1/(160*q30_factor)
wgs_df[wgs_df$group == "CODEC",]$cost_per_mb = wgs_df[wgs_df$group == "CODEC",]$cost_per_mb * 213657377/305459440


binom.confint(2188842, n=74362152460, conf.level = 0.95, method='wilson')
hifi = list(group="Pacbio HiFi", sample="mutation_metrics", lane="NA12878", sample_id = "NA12878.mutation_metrics", 
         total_bases_eval = 30 * 3.1e9, total_base_errors=NA, mean=2.94349e-05, upper=2.947391e-05, lower=2.939592e-05)
hifi$cost_per_mb = 13001 / (hifi$total_bases_eval / 1e6)
wgs_df = rbind(wgs_df, hifi)
ggplot(wgs_df, aes(x = cost_per_mb, y = mean, color = group)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
  geom_errorbar(aes(xmin = cost_per_mb/1.1, xmax = cost_per_mb/0.9), width = 0.5) +
  scale_x_log10(breaks = c(1e-3, 1e-2, 0.1, 1, 10, 100), labels = c("$0.001", "$0.01","$0.1", "$1", "$10", "$100"), limits=c(1e-3, 100)) +
  scale_y_log10(breaks = c(1e-6, 1e-5, 1e-4, 1e-3), labels = trans_format("log10", math_format(10^.x)), limits=c(1e-7, 1e-3)) +
  annotation_logticks(sides = "l", outside = F) +  coord_cartesian(clip = "off", xlim = c(0, 1e3)) + 
  annotation_logticks(sides = "b", outside = F) +  coord_cartesian(clip = "off", ylim = c(1e-7, 1e-3)) +
  scale_color_manual(values=c("deepskyblue3", "darkorange", "aquamarine3", "black")) +
  ylab("Error Rate per Base") + xlab("Cost per Million Base") + 
  theme(legend.title=element_blank()) + theme(legend.position = c(0.7, 0.74))

###Fig2e, Sperm
duprep_mono = MonomerErrorRate("figure2/e", group = "Duplex-Repair", patter=".*duprep.*mutation_metrics.txt", group_by_sample = TRUE)
codec_mono = MonomerErrorRate("figure2/e", group = "Commerical ER/AT", patter=".*R[1-2].mutation_metrics.txt", group_by_sample = TRUE)
codec_nanoseq = MonomerErrorRate("figure2/e", group="ddBTP-blocked", pattern="SpermNanoseq.*mutation_metrics.txt", group_by_sample = TRUE)

plotdf = rbind(codec_mono, duprep_mono, codec_nanoseq)
plotdf$muttype = factor(plotdf$muttype, levels =c("Overall","C>A", "C>G", "C>T", "T>A", "T>C", "T>G")) 
plotdf$group = factor(plotdf$group, levels =c("Commerical ER/AT", "Duplex-Repair", "ddBTP-blocked")) 
plotdf %>% filter(muttype == "Overall")

ggplot(plotdf, aes(x = muttype, y = mean, color = group)) +
  geom_point(size = 1, position=position_dodge2(width=0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5, position=position_dodge2(width=1)) +
  scale_y_log10(breaks = c(1e-9, 1e-8, 1e-7, 1e-6, 1e-5), labels = c("1/1B","1/100M", "1/10M", "1/1M", "1/100k")) +
  annotation_logticks(sides = "l", outside = T) +  coord_cartesian(clip = "off", ylim = c(1e-9, 1e-5)) + 
  labs(x = "Mononucleotide Mutation Type", y = "SNV Rate", color = "", shape = "") +
  #scale_x_discrete(labels=c(nanoseq_bb_mono$muttype)[1:6], "Overall")) +
  scale_color_manual(values=c("black", "firebrick", "forestgreen"), labels=c("Commerical ER/AT", "Duplex-Repair","ddBTP-blocked")) +
  theme_classic(20) +
  theme(legend.position = "bottom")

######Figure 3
#################
######################

###3b Germline
gatk = fread("figure3/b//wgs_low_pass_gatk.tsv")
g1 = ggplot(gatk %>% filter(depth < 6), aes(y=FNR, x=as.factor(depth), fill=Method)) + geom_bar(stat="identity", position=position_dodge2()) + xlab("")
g2 = ggplot(gatk %>% filter(depth < 6), aes(y=FPPM, x=as.factor(depth), fill=Method)) + geom_bar(stat="identity", position=position_dodge2()) + scale_y_reverse()
ggarrange(g1, g2, nrow=2, ncol=1)


##3c mutation rate vs age in buffy-coat samples
bc_nmut = fread("figure3/c/bc_mut_tables.tsv")
simp = lm(snv_rate ~ as.numeric(`Age at Diagnosis`), data=bc_nmut %>% filter(!str_detect(sample, "GBM")))
rsq = round(summary(simp)$r.square, digits=1)

slope = round(summary(simp)$coefficients[2,][1] * 6e9, digits=1)
intercept = round(summary(simp)$coefficients[1,][1] * 6e9,digits=1)

ggplot(bc_nmut, aes(y=snv_rate, x=`Age at Diagnosis`)) + geom_point(shape=23, size=3, fill="black") + 
  ylab("Redisual SNV / 10 Mb") +
  geom_abline(intercept=simp$coefficients[1],slope=simp$coefficients[2]) + 
  geom_abline(intercept=142.1/6e9, slope=19.8/6e9, color='red') +
  scale_y_continuous(breaks=c(1e-7, 1.5e-7, 2e-7, 2.5e-7, 3e-7, 3.5e-7), labels=c(1, 1.5, 2, 2.5, 3, 3.5), limits=c(1e-7, 4e-7)) +
  annotate(geom = "text", x=50, y=3.5e-7, label=paste0("R-square: ", rsq ,"\n slope: ", slope, "/ year / cell \n intercept: ", intercept," mut / year / cell"), size=4) +
  annotate(geom = "text", x=60, y=1.5e-7, label="Abascal et al. (granulocytes) \n slope : 19.8mut/year/cel\n intercept: 142 mut / year / cell", size=4, col='red') +
  theme_pubr(12)

###3d
### filter PMBC mutations
highinpt_PBMC_fingerprint = read_accu_files("figure3/d", "miredas_detect_fingerprint.txt")
highinpt_PBMC_fingerprint = highinpt_PBMC_fingerprint %>% mutate(sample = str_replace(sample, "HD_1", "HD1"))
highinpt_PBMC_fingerprint$sample = sapply(highinpt_PBMC_fingerprint$sample, function(x) {unlist(str_split(x, "_"))[1]})
highinpt_PBMC_fingerprint_merge = highinpt_PBMC_fingerprint %>% group_by(sample, contig, position, id, ref_allele, alt_allele) %>% summarise(total_alt_pass= sum(n_alt_pass), total_pass = sum(depth_pass))
highinpt_PBMC_fingerprint_vald = highinpt_PBMC_fingerprint_merge %>% filter(total_alt_pass > 0)  
ggplot(data=highinpt_PBMC_fingerprint_vald, aes(x=sample, y=total_alt_pass/total_pass)) + geom_violin() + geom_point(position=position_dodge2(width=0.2)) + 
  scale_y_log10(breaks = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1), labels = trans_format("log10", math_format(10^.x))) 

###3e
fig3e_sigs = get_sigs(highinpt_PBMC_fingerprint_vald, cosmic=cosmic_signatures_v3.2)
par(mfrow=c(4,2))
plot_signature(fig3e_sigs$spec, max_y=0.4, mutlab = TRUE)

###Figure 4a
#first download 8 breast cancer 80x WGS Mutect2 MAF files to figure4/mutect2"
###CODEC
bc8_codec_recall_precision = process_recall_prec("figure4/a/codec", 
                                     "figure4/mutect2", 
                                     codec_pat = ".variants_called.txt",
                                     mutect2_pat = "-filtered.annotated.maf")

####Stdnard NGS first 8 tumor
standdard_recall_precision = process_recall_prec("figure4/a/standard", 
                                     "figure4/mutect2", 
                                     codec_pat = ".variants_called.txt",
                                     sep="\\.",
                                     idx=1,
                                     mutect2_pat = "-filtered.annotated.maf")

totdf = rbind(standdard_recall_precision$prec %>% filter(mut == "SNV"), bc8_codec_recall_precision$prec %>% filter(mut == "SNV"))
totdf = totdf %>% mutate(method=c(rep("Standard", 16), rep("CODEC", 16)))
ggplot(totdf %>% filter(concordant == "high"), aes(x=s, y=perc, fill=method)) + geom_bar(stat="identity", position = position_dodge2()) + 
  scale_x_discrete(labels=c(1,2,3,4,6,7,9,10)) + ylab("Precision") + theme_pubr(20) + xlab("Patient") + scale_y_continuous(breaks=c(0, 0.2, 0.4, 0.6, 0.8), limits=c(0, 0.8))

#####Figure 4h,g
### Mutation signatures
brac_cosmicv2 = cosmic_signatures_v2[c(1,2,3,5,6,8,13,17,18,20,26,30),]
bc_mutect2_cosmic2 = process_folder("figure4/mutect2", "-filtered.annotated.maf$", sep = "-f", idx = 1, context_count = F,
                                      cosmic = cosmic_signatures_v2)
codec_bc8 = process_folder("figure4/a/codec", ".*variants_called.txt$", sep = "\\.", idx = 1, context_count = T,
                               cosmic = brac_cosmicv2)
standard_bc8 = process_folder("figure4/a/standard", ".*variants_called.txt$", sep = "\\.", idx = 1, context_count = F,
                               cosmic = brac_cosmicv2)
codec_cossim = pairwise_cossim(codec_bc8$mutdist, bc_mutect2_cosmic2$mutdist)
codec_cossim$method = "CODEC"
standard_cossim = pairwise_cossim(standard_bc8$mutdist, bc_mutect2_cosmic2$mutdist)
standard_cossim$method = "Standard"
ggplot(rbind(codec_cossim, standard_cossim), aes(x=left, y=cossim, fill=method)) + geom_bar(stat="identity", position = position_dodge2()) + scale_y_continuous(breaks=1:10/10) +
  scale_x_discrete(labels=c(1,2,3,4,6,7,9,10)) + ylab("Cosine Similarity") + theme_pubr(20) + xlab("Patient")

codec_bc8_weight = codec_bc8$weights
codec_bc8_weight$sample = rownames(codec_bc8$weights)
codec_bc8_weight$method = "CODEC"
codec_bc8_weight$mutect2weight = bc_mutect2_cosmic2$weights[,"Signature_3"]
standard_bc8_weight = standard_bc8$weights
standard_bc8_weight$sample = rownames(standard_bc8$weights)
standard_bc8_weight$method = "standard"
standard_bc8_weight$mutect2weight = bc_mutect2_cosmic2$weights[,"Signature_3"]

codeclm = lm(codec_bc8$weights[,"Signature_3"]~ bc_mutect2_cosmic2$weights[,"Signature_3"])
standardlm = lm(standard_bc8$weights[,"Signature_3"]~ bc_mutect2_cosmic2$weights[,"Signature_3"])
pc = cor(codec_bc8$weights[,"Signature_3"], bc_mutect2_cosmic2$weights[,"Signature_3"])
ps = cor(standard_bc8$weights[,"Signature_3"], bc_mutect2_cosmic2$weights[,"Signature_3"])
summary(codeclm)
summary(standardlm)

ggplot(rbind(codec_bc8_weight, standard_bc8_weight), aes(y=Signature_3, x=mutect2weight, color=method)) + geom_point() + 
  geom_abline(intercept=codeclm$coefficients[1],slope=codeclm$coefficients[2], color="pink") + 
  geom_abline(intercept=standardlm$coefficients[1],slope=standardlm$coefficients[2], color="skyblue") + 
  ylim(0, 0.4) + xlim(0, 0.4) +
  ylab("CODEC or Standrad 1-2x") + xlab("Mutect2") + theme_classic(12) + 
  theme(legend.position = "top") +
  annotate("text", x= 0.1, y=0.2, label= paste0("CODEC Pearson's r = ", pc, "\n p = ", summary(codeclm )$coefficients[2,4])) + 
  annotate("text", x= 0.3, y=0.05, label= paste0("Standard Pearson's r =", ps ,"\n p = ", summary(standardlm )$coefficients[2,4]))


####4b
#first download 175142_ABS_MAF.short.vcf.highcmplexregion.maf and 175142.filtered.annotated.maf DUOS 

##sensitivity
###alt>=1
coverage = fread("figure4/b/sensitivity_vs_depth.tsv")
high_coverage2 = process_recall_prec("figure4/b", 
                                     "figure4/mutect2", 
                                     sep = "_T",
                                     codec_pat = "*.variants_called.txt",
                                     mutect2_pat = "_ABS_MAF.short.vcf.highcmplexregion.maf")
###alt>=2
high_coverage_alt2 = process_recall_prec("figure4/b/maf_alt2", 
                                     "figure4/mutect2", 
                                     sep = "_T",
                                     codec_pat = "175142_T.*ds.variants_called.txt.filtered.maf",
                                     mutect2_pat = "_ABS_MAF.short.vcf.highcmplexregion.maf")

#####Precision
high_coverage_allmut = process_recall_prec("figure4/b", 
                                     "figure4/mutect2", 
                                     sep = "_T",
                                     codec_pat = "175142_T.*.variants_called.txt",
                                     mutect2_pat = "-filtered.annotated.maf")

high_coverage_allmut_alt2 = process_recall_prec("figure4/b/maf_alt2", 
                                     "figure4/mutect2", 
                                     sep = "_T",
                                     codec_pat = "175142_T.*.filtered.maf",
                                     mutect2_pat = "-filtered.annotated.maf")

alt1recall = high_coverage2$recall %>% filter(mut == "SNV" & concordant)
alt1recall$cutoff="alt>=1"
alt2recall = high_coverage_alt2$recall %>% filter(mut == "SNV" & concordant)
alt2recall$cutoff="alt>=2"

prec_alt1 = high_coverage_allmut$prec %>% filter(mut == "SNV" & concordant=="high")
prec_alt1$cutoff="alt>=1"

prec_alt2 = high_coverage_allmut_alt2$prec %>% filter(mut == "SNV" & concordant=="high")
prec_alt2$cutoff="alt>=2"

codec_recall = rbind(alt1recall, alt2recall)
codec_recall = codec_recall %>% mutate(data="CODEC", type="Recall") %>% select(perc, cutoff, type, data)
codec_prec = rbind(prec_alt1, prec_alt2)
codec_prec = codec_prec %>% mutate(data="CODEC", type="PPV") %>% select(perc, cutoff, type, data)
codec_all = rbind(codec_recall, codec_prec)
codec_all$PF_cov = coverage$PF_cov 
codec_all$cov = coverage$cov 


theory = data.frame(perc = 1-dbinom(0, size=1:10, prob=0.86*0.5), cutoff = "alt>=1", PF_cov=1:10, type="Recall", data="Theory", cov=NA)
suppalt1 = data.frame(perc = 1-pbinom(0, size=4:10, prob=0.86*0.5*0.75), cutoff = "alt>=1", PF_cov=4:10, type="Recall", data="predict", cov=NA)
suppalt2 = data.frame(perc = 1-pbinom(1, size=4:10, prob=0.86*0.5*0.9), cutoff = "alt>=2", PF_cov=4:10, type="Recall", data="predict", cov=NA)
allall = rbind(codec_all, suppalt1, suppalt2, theory)
xlabels = round(coverage$PF_cov, digits=1)
ggplot(allall %>% filter(type=="Recall"), aes(x=as.numeric(PF_cov), y=perc, color=factor(data,c("CODEC", "Theory", "predict")))) + geom_point() + 
  geom_line(aes(group=interaction(factor(data,c("CODEC", "Theory", "predict")), cutoff), )) +
  #geom_line(data=allall%>%filter(type=="PPV"), aes(x=PF_cov, y=perc, group=cutoff), color="deepskyblue3", linetype='dotted') +
  geom_point(data=allall%>%filter(type=="PPV"), aes(x=PF_cov, y=perc, group=cutoff), , shape=24) +
  geom_point(aes(x=3.7, y=0.986), color="#b64848") +
  geom_point(aes(x=3.7, y=0.928), color="#b64848", shape=24) +
  #geom_hline(yintercept = 0.95, linetype = "dashed") +
  scale_x_continuous(breaks=c(1.1, 1.9, 2.6, 3.7, 5:10), labels=c("1\n5", "1.9\n9", "2.6\n12", "3.7\n18", "5\n25", "6\n30", "7\n35", "8\n40", "9\n45", "10\n50")) +
  scale_y_continuous(
    breaks=1:10/10, limits=c(0,1.0), name="PPV",
    sec.axis = sec_axis(~., breaks=1:10/10, name = "Sensitivity")
  )+
  scale_color_manual(values=c("deepskyblue3", "azure2", "darkorange"), labels=c("CODEC", "Theoretical", "Projection"), name="")+
  xlab("coverage") + ylab("Sensitivity") +
  theme_pubr(7) 

############
##MSI
#### Fig 4d,e
#From the terra workspace download gs://fc-8585a400-07f7-43c5-a267-bb13476c2b09/ct19/figure4/

good_sigs = c("SBS15", "SBS21", "SBS26", "SBS1")

codec_tumor = process_folder("figure4/d/codec",
               ".variants_called.txt", 
               sep=".bam.",
               idx=2, 
               cosmic = cosmic_signatures_v3.2,
               context_count = F)
plot_signature(codec_tumor$mutdist["1.0000",], mutlab = TRUE, xticklab = TRUE)

mutect2 = process_folder("figure4/d/mutect2",
               "x.*.maf$", 
               sep="\\.", 
               idx=6,
               context_count = F,
               cosmic = cosmic_signatures_v3.2)

plot_signature(mutect2$mutdist["12x",], mutlab = TRUE, xticklab = TRUE)

std_allmut = process_folder("figure4/d/standard",
               ".variants_called.txt$", 
               sep="\\.", 
               idx=2,
               context_count = F,
               cosmic = cosmic_signatures_v3.2)

plot_signature(std_allmut$mutdist["12x",], mutlab = TRUE, xticklab = TRUE)

prepare_ds_plot <- function(test, true, depthidx, good_sigs=c("SBS15", "SBS21", "SBS26"), group="") {
  cossim = apply(test$mutdist, 1, function(x) {cosine_sim(as.numeric(x), as.numeric(true))})
  print(names(cossim))
  a = sapply(names(cossim), function(x){unlist(strsplit(x, "_(?![^_]*_)", perl=TRUE))})
   
  df = data.frame("group" = group, "depth"=as.numeric(str_replace(a, "x", "")), "cossim"=cossim, "sum_relevant_weight" = rowSums(test$weights[, good_sigs]))
  df = arrange(df, depth) 
  df$"x" = 1:length(cossim)
  df
}

dfdsall = data.frame()
dfdsall = rbind(dfdsall, prepare_ds_plot(codec_tumor, mutect2$mutdist[12,], group="CODEC"))
dfdsall = rbind(dfdsall, prepare_ds_plot(mutect2,  mutect2$mutdist[12,], group="MuTect2"))
dfdsall = rbind(dfdsall, prepare_ds_plot(std_allmut,  mutect2$mutdist[12,], group="Standard"))

dfdsall %>% ggplot(., aes(x=depth, y=cossim, color=group)) + geom_point() + geom_line() +
  #scale_x_continuous(breaks=seq(0, 2, 0.1), limits=c(0, 2)) +
  scale_x_log10(limits=c(1e-4,13)) +
  scale_y_continuous(breaks=seq(0, 1, 0.1)) +
  theme_classic(base_size = 20) + 
  xlab("Depth") + ylab("Cosine Similarity") +
  scale_color_manual(values=c("deepskyblue3", "black", "grey")) +
  theme(legend.position="top", legend.title = element_blank())

##Fig 4f
###Mutect2
msimutect2_maf = read_maf("figure4/d/mutect2/MSI_tumor.stdngs.aligned.dupliecates_marked.bam.12x.filtered.maf")
msimutect2 = process_folder("figure4/d/mutect2/", "MSI_tumor.stdngs.aligned.dupliecates_marked.bam.12x.filtered.maf", sep="\\.", idx=1, cosmic=cosmic_signatures_v3.2)
plot_weights(msimutect2$weights, good_sigs = good_sigs, "MSI")

##CODEC tumor all mutations
msitumor = read_accu_files(reformat_on_mac("figure4/d/codec", "ct19_tumor1.merged.cds_consensus.mol_consensus.aligned.bam.2.5000.bam.variants_called.tx"), "tumor.variants_called.txt")
msitumor_sigs = get_sigs(msitumor, cosmic = cosmic_signatures_v3.2)
plot_weights(msitumor_sigs$weights, good_sigs = good_sigs, "MSI")

### remove mutect2 variants
msitumor = add_varid(msitumor)
msimutect2_maf = add_varid(msimutect2_maf)
msitumor_notinm2 = msitumor %>% filter(!id %in% msimutect2_maf$id)
msitumor_notinm2_sigs = get_sigs(msitumor_notinm2, cosmic = cosmic_signatures_v3.2)
plot_weights(msitumor_notinm2_sigs$weights, good_sigs = good_sigs, "MSI")

###CODEC adj-normal
msinormal = process_folder("figure4/f", "normal.variants_called.txt", sep="\\.", idx=1, cosmic=cosmic_signatures_v3.2)
plot_weights(msinormal$weights, good_sigs = good_sigs, "MSI")

###standard NGS all mutations
msistdngs = process_folder("figure4/d/standard", "2x.variants_called.txt", sep="\\.", idx=1, cosmic=cosmic_signatures_v3.2)
plot_weights(msistdngs$weights, good_sigs = good_sigs, "MSI")


#####Figure 5
#check CODECsuite/msi for Figure 5