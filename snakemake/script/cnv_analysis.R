###FILE PATH: /xchip/bloodbiopsy/ruolin/link_duplex/snakemake/script/
### PreProcess:
## /xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/merged_pcr_condition_v2/cnv -r $HG19 -b $BAM -L $BED -m minmapq -V phased.vcf -q minbaseq
######
bins = readRDS(reformat_on_mac("/xchip/bloodbiopsy/ruolin/nucleosome/frag_pattern/bins.rds"))
newStyle <- mapSeqlevels(seqlevels(bins), "NCBI")
bins <- renameSeqlevels(bins, newStyle)
bins
library(rtracklayer)
library(ggplot2)
#export.bed(bins, con=reformat_on_mac("/xchip/bloodbiopsy/ruolin/nucleosome/frag_pattern/100kbbins.bed"))
binsdf <- data.frame(seqnames=seqnames(bins),
                 starts=start(bins)-1,
                 ends=end(bins),
                 domain = bins$domain, 
                 arm = bins$arm,
                 gc = bins$gc)
binsdf$seqnames = (as.character(binsdf$seqnames))
#write_tsv(binsdf, path = reformat_on_mac("/xchip/bloodbiopsy/ruolin/nucleosome/frag_pattern/100kbbins.bed"), col_names = FALSE)

#cnvdf = fread(reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/pcr_condition_v2_hiseq_250x_NA12878/wrk/regular.cnv.lowprec.out"))
cnvdf = fread(reformat_on_mac("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/merged_pcr_condition_v2/cnv/na12878_regular_4x.tsv"))
#cdscnvdf
cnvdf$gr = paste(cnvdf$chrom, cnvdf$begin, cnvdf$end, sep=":")
cnvdf = cnvdf %>% mutate(BAF = nB / (nA + nB))
cnvdf$chrom = as.character(cnvdf$chrom)
cnvdffull = cnvdf %>% left_join(binsdf, by = c("chrom"="seqnames", "begin"="starts", "end"="ends" ) )

binsize=200
cnvdffull = cnvdffull %>% group_by(arm) %>% 
    mutate(combine = ifelse(grepl("p", arm), ceiling(((1:length(arm))/binsize)), 
                            ceiling(rev(1:length(arm))/binsize)))
plotdf = cnvdffull %>%  filter(nloci >= 10) %>% filter(abs(BAF-0.5) < 0.1) %>%
  group_by(chrom, arm, combine) %>% summarise(n=n(), BAF_eact = sum(nA)/ (sum(nA) + sum(nB)), BAF_median = median(BAF, na.rm=TRUE))
plotdf$x = paste(plotdf$arm, plotdf$combine, sep='.')
plotdf = plotdf %>% mutate(mBAF=abs(BAF_eact - 0.5) + 0.5)
#png(filename="~/Documents/regular_4x_cnv.1.png", width=10, heigh=8, unit='in', res=300)
ggplot(data=plotdf, aes(x=x, y=mBAF)) + geom_point() + theme_classic(base_size= 20) +  
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=12)) + xlab("Chr1, 5Mb bins") +
  #scale_y_continuous(breaks=seq(0.5, 0.60, 0.01), limits=c(0.49, 0.6)) + 
  scale_y_continuous(breaks=seq(0.5, 0.75, 0.05), limits=c(0.49, 0.75)) + 
  geom_hline(yintercept = hemloss, color="red") +   geom_text(aes(10, hemloss,label = "hemi loss, TF = 0.12", vjust=-0.5)) +
  geom_hline(yintercept = scgain, color="blue") + geom_text(aes(25, scgain,label = "single copy gain, TF=0.13", vjust=-0.4)) +
  geom_hline(yintercept = cnloh, color="green") + geom_text(aes(40, cnloh, label = "CN_LOH, TF=0.06", vjust=1)) 
dev.off()
  

nrow(plotdf)

cn_loh <- function(x) {
  (2 - x)/2
}
cn_loh(0.95)

hem_loss <- function(x) {
  1 / (1 + x)
}
sc_gain <- function(x) {
  (2 - x)/(3 - x)
}

hemloss = hem_loss(0.88)
scgain = sc_gain(0.87)
cnloh = cn_loh(0.99)


log2(1)
