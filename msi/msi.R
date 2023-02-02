library("rjson")
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

empdis <- function(l, nuc, el, profile, sd=0.5) {
  #el: expect length
  r = dnorm(l, mean=el, sd=sd, log=T) - log(sum(dnorm(0:(2*el), mean=el, sd=sd)))
  exp(r)
}

delep <- function(nuc, el, profile, ndel=1) {
  rep = paste(nuc, el, sep= ',') 
  rcrept = paste(rev.comp(nuc), el, sep=',')
  dat = addList(profile[[rep]], profile[[rcrept]])
  
  corr = dat[[as.character(el)]]
  tot = corr
  for (ii in names(profile[[rep]])) {
    if (as.numeric(ii) < el) {
      tot = tot + dat[[ii]]
    }
  }
  minus = dat[[as.character(el-ndel)]]
  minus/tot
}

# for (i in 10:25) {
#   print(paste(i, delep('A', i, profile, 2)))
# }

empdis2 <- function(profile, l, nuc, el) {
  # fit a normal distribution. But currently the fitting is bad 
  opt <- function(sigma, dat, rl) {
    #dist: dist for a single repeat type
    #rl: repeat length 
    s = sum(unlist(dat))
    idx = 1
    cost = 0
    scale = 1
    for (cnt in dat) {
      obl = as.numeric(names(dat)[idx])
      idx = idx+ 1
      #print(names(profile$"AT,9"[ii]))
      cost = cost + scale* abs(empdis(obl, nuc, rl, profile, sd = sigma) - cnt/s)
    }
    cost
  }
  
  rept = paste(nuc, el, sep=',')
  dat = profile[[rept]]
  sigma = optimize(opt, c(0,1), dat=dat, rl=el)
  empdis(l,nuc, el, profile, sd=sigma$minimum)
}

addList <- function(lista, listb)  {
  if (is.null(lista)) {
    lista = list()
  }
  if (is.null(listb)) {
    listb = list()
  }
  for (i in names(lista)) {
    if (i %in% names(listb)) {
      listb[[i]] = listb[[i]] + lista[[i]]
    } else {
      listb[[i]] = lista[[i]] 
    }
  }
  listb
}

rev.comp<-function(x,rev=TRUE)
{
  #https://www.r-bloggers.com/2008/11/r-function-to-reverse-and-complement-a-dna-sequence/
  x<-toupper(x)
  y<-rep("N",nchar(x))
  xx<-unlist(strsplit(x,NULL))
  for (bbb in 1:nchar(x))
  {
    if(xx[bbb]=="A") y[bbb]<-"T"        
    if(xx[bbb]=="C") y[bbb]<-"G"        
    if(xx[bbb]=="G") y[bbb]<-"C"        
    if(xx[bbb]=="T") y[bbb]<-"A"
  }
  if(rev==FALSE) 
  {
    for(ccc in (1:nchar(x)))
    {
      if(ccc==1) yy<-y[ccc] else yy<-paste(yy,y[ccc],sep="")
    }
  }
  if(rev==T)
  {
    zz<-rep(NA,nchar(x))
    for(ccc in (1:nchar(x)))
    {
      zz[ccc]<-y[nchar(x)+1-ccc]
      if(ccc==1) yy<-zz[ccc] else yy<-paste(yy,zz[ccc],sep="")
    }
  }
  return(yy)    
}

empdis3 <- function(l, nuc, el, profile, correct = FALSE, merge_rc = FALSE) {
  #l: observe len 
  #nuc: repeat nucleotides
  #el: expected len
  #profile: msi profile
  ### just use the emprical distribution. Not fitting a model.
  rept = paste(nuc, el, sep=',')
  rcrept = paste(rev.comp(nuc), el, sep=',')
  if (merge_rc) {
    if (!(rept %in% names(profile)) & !(rcrept %in% names(profile))) {
      return (empdis(l,nuc, el, profile, sd=0.5))
    }
  } else {
    if (!(rept %in% names(profile))) {
      return (empdis(l,nuc, el, profile, sd=0.5))
    }
  }
  if (merge_rc) {
    dat1 = copy(profile[[rept]])
    dat2 = copy(profile[[rcrept]]) 
    dat = addList(dat1, dat2)
  } else {
    dat = copy(profile[[rept]])
  }
  if (sum(unlist(dat)) < 500) {
    return (empdis(l,nuc, el, profile, sd=0.5))
  }
  lens_char = names(dat)
  if (correct) {
    for (cur in lens_char){
      obl = as.numeric(cur)
      if (obl < el) {
        last = as.character(obl+1)
        if (!(last %in% lens_char) || dat[[cur]] > dat[[last]]) {
          dat[[cur]] = 0 
        }
      }
      if (obl > el) {
        last = as.character(obl-1)
        if (!(last %in% lens_char) || dat[[cur]] > dat[[last]]) {
          dat[[cur]] = 0 
        }
      }
    }
  }
  if (as.character(l) %in% lens_char) {
    as.numeric(dat[[as.character(l)]] / sum(unlist(dat)))
  } else {
    empdis(l,nuc, el, profile, sd=0.5)
    #0.0
  } 
}



NormGtPrior <- function(gt, b = 0.5) {

    #gt:  0 homozygous reference, 1 het
  #b: fraction of germline indel in a population
    return (ifelse(gt, b, 1-b))
}

TumGtPrior <- function(gt, c=0.5) {

    #gt:  0 homozygous reference, 1 het
  #c: prob of MSI-H
    return (ifelse(gt, c, 1-c))
}


NormLik <- function(y, gt, log=FALSE) {
  #gt:  0 homozygous reference, 1 het
  #y num. wildtype in normal sample  
  if (y == 0) {logp = 0}
  else {logp = ifelse(gt, y*log(0.5), 0)}
  ifelse(log, logp, exp(logp))
} 

HaploLik <- function(h, gt) {
  #gt:  0 homozygous reference, 1 het
  #h: 0 wildtype haplotype, 1 mutant haplotype
  if (gt) {
    0.5  
  } else {
    ifelse(h, 1, 0) 
  }
} 
SingleReadLik <- function(rept, l, gt, H1, H2, profile) { #Given genotype
  #rept: repeat seq
  #l: observed hp len
  #gt:  0 homozygous reference, 1 het
  #H1: reference hyplotype in homopolyer length
  #H2: alternative hyplotype in homopolyer length
  pH1 = ifelse(gt, 0.5, 1)
  pH2 = 1- pH1
  empdis3(l, rept, H1, profile) * pH1 + empdis3(l, rept, H2, profile) * pH2
} 


SomaLik <- function(rept, L, a, G_n, G_t, H1, H2, profile, log=FALSE) {
  #rept: repeat seq
  #L: observed hp length in tumor
  #a: fraction of normal cells in tumor sample
  #G: Genotypes (G_n, G_t) genotype for normal and tumor
  #H1: reference hyplotype in homopolyer length #H2: alternative hyplotype in homopolyer length
  sum = 0
  for (l in L) {
    sum = sum + log(SingleReadLik(rept, l, G_n, H1, H2, profile) * a + SingleReadLik(rept, l, G_t, H1, H2, profile) * (1-a))
  }
  ifelse(log, sum, exp(sum))
}


DataLik <- function(rept, L, y, a, G_n, G_t, H1, H2, profile, log=FALSE) {
  #rept: repeat seq
  #L: observed hp length in tumor
  #y: num. wildtype reads in normal
  #a: fraction of normal cells in tumor sample
  #G: Genotypes (G_n, G_t) genotype for normal and tumor
  #H1: reference hyplotype in homopolyer length
  #H2: alternative hyplotype in homopolyer length
  p = SomaLik(rept, L, a, G_n, G_t, H1, H2, profile, log=TRUE) + NormLik(y, G_n, log=TRUE) 
  ifelse(log, p, exp(p))
}

Poster_s <- function(rept, L, y, a, g_n, g_t, H1, H2, b, c, profile, log=FALSE) {
  logp = DataLik(rept, L, y, a, g_n, g_t, H1, H2, profile, log=T) +  log(NormGtPrior(g_n, b)) + log(TumGtPrior(g_t, c)) 
  ifelse(log, logp , exp(logp))
}

Poster <- function(rept, L, y, a, H1, H2, b, c, profile, verbose=F) {
  #rept: repeat seq
  #L: msi-informative repeat lengths in tumor
  #y: num. wildtype reads in normal
  #H1: haplotype 1, normal repeat expected length
  #H2: haplotype 2, tumor repeat expected length
  #a: fraction of normal cells in tumor sample
  #b: fraction of germline indel in a population
  #c: prob. of MSI-H
  
  #G=AA, S=AB 
  p1 = Poster_s(rept, L, y, a, 0, 1, H1, H2, b, c, profile, log=T) 
  den = 0
  for(g_n in 0:1) {
    for(g_t in 0:1) {
      pp = Poster_s(rept, L, y, a, g_n, g_t, H1, H2, b, c, profile)
      if (verbose) {
        print(paste(g_n, g_t, pp))
      }
      den = den + pp
    }
  }
  exp(p1 - log(den))
}


#Poster("T", L=c(6,6,6), y=2, a=0.9, H1=10, H2=6, b=0.0, c=0.5, profile=profile)
#Poster("A", L=c(13), y=0, a=0.9999, H1=17, H2=13, b=1e-5, c=0.5, profile=profile)
#Poster("A", L=c(10), y=0, a=0.9999, H1=13, H2=10, b=1e-5, c=0.5, profile=profile)



Main <- function(sample_file, str_profile, nf, verbose=0) {
  nf = min(nf, 1-1e-4)
  total_P = 0  
  profile = fromJSON(file=reformat_on_mac(str_profile))
  df = fread(reformat_on_mac(sample_file), header=F, col.names = c("chrom", "pos", "rept", "length" ,"tumor", "normal", "prefix", "suffix", "popAF", "germline"))
  if (nrow(df) == 0) return (total_P)
  for (i in 1:nrow(df)) {
    if (str_detect(df[i,]$tumor, ",")) {
      L = as.numeric(unlist(strsplit(df[i, ]$tumor ,",")))
    } else {
      L = as.numeric(df[i, ]$tumor)
    }
    if (is.na(df[i,]$normal)) {
      y = 0
      H1 = df[i,]$length
    } else {
      if (str_detect(df[i,]$normal, ",")) {
        normLen = as.numeric(unlist(strsplit(df[i, ]$normal ,",")))
      } else {
        normLen = as.numeric(df[i, ]$normal)
      }
      y = length(normLen)
      H1 = getmode(normLen)    
    }
    H2 = min(L)    
    rept = df[i,]$rept
    p = Poster(rept, L=L, y=y, a=nf, H1=H1, H2=H2, b=df[i,]$popAF, c=0.5, profile=profile, verbose > 1)
    if (is.nan(p)) {
      next
    }
    total_P = total_P + p
    if (verbose) print(paste(df[i,]$rept, df[i,]$tumor, df[i,]$normal, p, sep='          '))
  }
  total_P
}

MsiDetect <- function(dirin, str_profile, verbose=F) {
  sid = c()
  for (f in list.files(reformat_on_mac(dirin))) {
    fields = unlist(strsplit(f, ".", fixed = TRUE))
    tf = as.numeric(paste0("0.", fields[2])) * 100
    sample = fields[3]
    sid = c(sid, paste(sample, tf))
  }
  result = data.frame("sid" = unique(sid), "all_filtered" = 0, "msi_filtered" = 0, "all" =0 , "msi" = 0, score=0)
  for (f in list.files(reformat_on_mac(dirin))) {
    germfilter = 0
    if (str_detect(f, "filtered")) {
      germfilter = 1
    }
    # if (str_detect(f, "all")) {
    #   allsite = 1
    # }
    fields = unlist(strsplit(f, ".", fixed = TRUE))
    tf = as.numeric(paste0("0.", fields[2])) * 100
    sample = fields[3]
    sid = paste(sample, tf)
    result[result$sid == sid, "sample"] = sample
    if (length(fields) == 5)  cat = paste(fields[4], fields[5], sep="_")
    else cat = fields[4]
    print(cat)
    if (germfilter) {
      df = fread(file.path(reformat_on_mac(dirin), f), header=F, col.names = c("chrom", "pos", "rept", "length" ,"tumor", "normal", "prefix", "suffix", "popAF", "germline"))
    } else {
      df = fread(file.path(reformat_on_mac(dirin), f), header=F, col.names = c("chrom", "pos", "rept", "length" ,"tumor", "normal", "prefix", "suffix", "popAF"))
    }
    if (nrow(df) == 0) {
      result[result$sid == sid, cat] = n
      if (cat == "msi_filtered") {
        result[result$sid == sid, "score"] = 0
      }
      next
    }
    df = df %>% filter(str_length(rept) == 1 & prefix < 0.7 & suffix < 0.7)
    n = nrow(df)
    result[result$sid == sid, cat] = n
    if (cat == "msi_filtered") {
      stopifnot(result[result$sid == sid, "all_filtered"] >0)
      est_tf = result[result$sid == sid, cat] * 8 / result[result$sid == sid, "all_filtered"] 
      score = Main(file.path(dirin, f), str_profile, 1- est_tf)
      result[result$sid == sid, "score"] = score
    }
  }
  result = result %>% separate(sid, c("sample", "tf"), sep=" ") %>% pivot_longer(c("msi", "msi_filtered"), names_to="Germ. Filtering")
  result$tf = as.numeric(result$tf)
  result
}
str_profile =""


###Compare standard vs CODEC
reg_result = MsiDetect("", str_profile)
reg_result$Library = "Standard NGS"
orig_result$Library = "CODEC"
df = rbind(reg_result, orig_result)
ggplot(data=df %>% filter(`Germ. Filtering`=="msi_filtered") , aes(x=sample, y=value, fill=Library)) + geom_bar(stat="identity", position=position_dodge2()) +
  theme_classic(24) + scale_y_log10() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  xlab("") + ylab("# Informative Fragment") + scale_fill_manual(values = c("deepskyblue3", "darkorange")) +
  theme(legend.position="bottom", axis.text = element_text(size = 22, color = "black"))



###NA12878
na12878 = MsiDetect("", str_profile)
Main("", str_profile, 1-1e-6, verbose=2)

#################END TEMP RESULTS


post_process <- function(dtable, class, lb= 7e-3) {
  tmp = dtable
  if ("Germ. Filtering" %in% colnames(dtable)) {
    tmp = dtable %>% filter(`Germ. Filtering` == "msi_filtered")
  }
  tmp$class = class
  tmp$observed = tmp$score / max(tmp$score) * max(tmp$tf) 
  tmp$observed = ifelse(tmp$observed < lb, lb, tmp$observed)
  tmp
}

t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=TRUE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}

get_pvalues <- function(plotdf, adjacent=F) {
  plotdf$pvalue = NaN
  if (adjacent) {
    for(i in 1:(nrow(plotdf)-1)) {
      res = t.test2( plotdf[i, ]$MSIscore, plotdf[i+1,]$MSIscore, plotdf[i,]$sd, plotdf[i+1, ]$sd, 3 ,3 )
      plotdf[i+1, ]$pvalue = res["p-value"]
    }
  } else {
    for (i in 1:nrow(plotdf)) {
      res = t.test2( plotdf[i, ]$MSIscore, plotdf[plotdf$tf == 0,]$MSIscore, plotdf[i,]$sd, plotdf[plotdf$tf == 0, ]$sd, 3 ,3 )
      plotdf[i, ]$pvalue = res["p-value"]
    }
  }
  plotdf
}

###LOAD ALL DATA!!!!!!!!!!!!

##Standard 
#######Dilution series below
#### Compare Model based vs count based


std2xr1 = MsiDetect("reg_2x_tumor1_normal1/result", str_profile)
std2xr2 = MsiDetect("reg_2x_tumor1_normal2/result", str_profile)
std2xr3 = MsiDetect("reg_2x_tumor2_normal1/result", str_profile)
#std2xr4 = MsiDetect("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/ct19_msih_tumor_normal/codec/msi/reg_2x_tumor2_normal2/result", str_profile)

std2x = MsiDetect("stdngs_dilu2x_diffnorm/result", str_profile)
codec = MsiDetect("codec_dilu1x_diffbam/result", str_profile)
codec2x = MsiDetect("codec_dilu2x_diffbam/result", str_profile)

#std2x_neg = MsiDetect("/xchip/bloodbiopsy/ruolin/link_duplex/cds_wgs/ct19_msih_tumor_normal/codec/msi/stdngs_addition_nc/result/2x", str_profile)

stdplot2xr1 = post_process(std2xr1, "Standard")
stdplot2xr2 = post_process(std2xr2, "Standard")
stdplot2xr3 = post_process(std2xr3, "Standard")
codecplot2x = post_process(codec2x, "Codec")

###ADD MsMutect
msmutectr1 = "msmutect/r1"
msmutectr2 = "msmutect/r2"
msmutectr3 = "msmutect/r3"
msmutect_parse <- function(msdir) {
  tfs = c()
  score = c()
  for (f in list.files(reformat_on_mac(msdir), "mut.tsv")) {
    fields = unlist(str_split(f, "\\."))
    tf = as.numeric(fields[2])
    tfs = c(tfs, tf)
    msdf = fread(file.path(reformat_on_mac(msdir), f))
    msdf = msdf %>% filter(CALL == 1)
    score = c(score, nrow(msdf))
  }
  msmutect_result = data.frame("score" = score, tf=tfs/100)
  msmutectdf = post_process(msmutect_result, "Standard_MSMutect")
  msmutectdf
}
msmutectdf1 = msmutect_parse(msmutectr1)
msmutectdf2 = msmutect_parse(msmutectr2)
msmutectdf3 = msmutect_parse(msmutectr3)
msmutect_totaldf = rbind(msmutectdf1, msmutectdf2, msmutectdf3)
msmutect_totaldf = msmutect_totaldf %>% group_by(class, tf) %>% summarise(MSIscore = mean(score), sd = sd(score), observed_tf = mean(observed), tfsd = sd(observed)) %>% ungroup()
msmutect_totaldf = get_pvalues(msmutect_totaldf, adjacent = F)

##For main plot. agg for all date
####compare stdngs vs codec
stdplot2xall = rbind(stdplot2xr1, stdplot2xr2, stdplot2xr3)
stdplot = stdplot2xall %>% group_by(class, tf) %>% summarise(MSIscore = mean(score), sd = sd(score), observed_tf = mean(observed), tfsd = sd(observed)) %>% ungroup()
stdplot = get_pvalues(stdplot, adjacent = F)

codecplot =  codecplot2x %>% group_by(class, tf) %>% summarise(MSIscore = mean(score), sd = sd(score), observed_tf = mean(observed), tfsd = sd(observed)) %>% ungroup()
codecplot$pvalue = NaN

plotdf = rbind(stdplot, codecplot, msmutect_totaldf)
plotdf$x = paste0(plotdf$tf,"%")

###################

### Dilution dot plot 
#V1:  expected vs observed tf

ggplot(data=plotdf %>% mutate(tf = ifelse(tf == 0, 5e-3, tf)), aes(x=tf, y=observed_tf, color = class)) + geom_point(size=3) + 
  geom_line(aes(group=class)) +
  theme_classic(24) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  geom_errorbar( aes(x=tf, ymin=observed_tf+tfsd, ymax=observed_tf-tfsd), width=0.3, colour="orange", alpha=0.9, size=0.5) +
  scale_color_manual(name="", values = c("deepskyblue3", "darkorange", "orchid3"), labels = c("CODEC", "Standard", "Standard_MSMutect")) +
  #scale_x_discrete(labels=unique(plotdf$x)) +
  scale_x_log10(breaks = plotdf$tf, labels = comma_format(accuracy = 0.01), limits=c(-1, 1e2)) +
  scale_y_log10(breaks = plotdf$tf, labels = comma_format(accuracy = 0.01), limits=c(1e-2, 1e2)) +
  annotation_logticks(sides = "l", outside = F) +  coord_cartesian(clip = "off", ylim = c(1e-2, 1e2)) + 
  xlab("Expected TF") + ylab("Observed TF") +
  geom_segment(x=-2, y=-2, xend=2, yend = 2, color='grey') +
  theme(legend.position="bottom", axis.text = element_text(size = 22, color = "black"))

#V2: barplot 
ggplot(data=plotdf, aes(x=as.factor(tf), y=MSIscore, fill = class)) + geom_bar(stat="identity", position=position_dodge2()) + 
  geom_errorbar( aes(x=as.factor(tf), ymin=MSIscore-sd, ymax=MSIscore+sd), width=0.9, colour="orange", alpha=0.9, size=1, position=position_dodge2()) +
  theme_classic(24) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(values = c("deepskyblue3", "darkorange"), labels = c("CODEC", "Standard")) +
  scale_x_discrete(labels=unique(plotdf$x)) +
  scale_y_log10(breaks = c(1e-1, 1, 10, 100, 1000, 1e4), labels = trans_format("log10", math_format(10^.x)), limits=c(1e-1, 1e4)) +
  annotation_logticks(sides = "l", outside = F) +  coord_cartesian(clip = "off", ylim = c(1e-1, 1e4)) + 
  xlab("") + ylab("Signal") +
  theme(legend.position="bottom", axis.text = element_text(size = 22, color = "black"))

#V3:
ggplot(data=plotdf, aes(x=as.factor(tf), y=MSIscore, color = class)) + 
  geom_line(aes(group=class)) +
  geom_point(size=3) + 
  geom_errorbar( aes(x=as.factor(tf), ymin=MSIscore-sd, ymax=MSIscore+sd), width=0.8, alpha=0.9, size=0.5) +
  theme_classic(24) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(name="", values = c("deepskyblue3", "darkorange", "orchid3"), labels = c("CODEC", "Standard", "Standard_MSMutect")) +
  scale_x_discrete(labels=unique(plotdf$x)) +
  scale_y_log10(breaks = c(1e-1, 1, 10, 100, 1000, 1e4), labels = trans_format("log10", math_format(10^.x)), limits=c(1e-1, 1e4)) +
  annotation_logticks(sides = "l", outside = F) +  coord_cartesian(clip = "off", ylim = c(1e-1, 1e4)) + 
  annotate("text", x = 2:10, y=300*(2:10), label= format(round(plotdf$pvalue[2:10], 4), nsmall=4)) +
  xlab("") + ylab("MSIscore") +
  theme(legend.position="bottom", axis.text = element_text(size = 22, color = "black"))

plotdf$pvalue = NaN
for(i in 1:(nrow(plotdf)-1)) {
  if (plotdf[i+1, ]$class == "Standard") {
    res = t.test2( plotdf[i, ]$MSIscore, plotdf[i+1,]$MSIscore, plotdf[i,]$sd, plotdf[i+1, ]$sd, 3 ,3 )
    plotdf[i+1, ]$pvalue = res["p-value"]
  }
}
View(plotdf)
###MISC


#Downsample standard NGS dilution series
dsbytf <- function(wanttf, dsdep, otdep, ondep, origtf, tumor, normal, output, dupaware) {
  ## wanttf: desired tumor fraction in output mixture
  ## otdep: original tumor depth
  ## ondep: original normal depth
  ## dsdep: desired depth in output mixture
  ## tumor: tumor bam
  ### normal: normal bam
  ### output: output prefix
  ### dupaware: duplicate aware 
  ds1 = wanttf * dsdep / origtf / otdep
  stopifnot(ds1 <= 1)
  ds2 = (dsdep - ds1*otdep)/ondep
  stopifnot(ds2 <= 1)
  ndigit = 5
  wanttf = format(round(wanttf, 4), nsmall=4, scientific = F)
  ds1 = format(round(ds1, ndigit), nsmall=ndigit, scientific = F)
  ds2 = format(round(ds2, ndigit), nsmall=ndigit, scientific = F)
  c(ds1, ds2)
  if (dupaware) {
    cm1 = paste0("picard_downsample ", tumor, " ", ds1, " ", tumor, ".", ds1, "ds true" )
    cm2 = paste0("picard_downsample ", normal, " ", ds2, " ", normal, ".", ds2, "ds true" )
  } else {
    cm1 = paste0("picard_downsample ", tumor, " ", ds1, " ", tumor, ".", ds1, "ds" )
    cm2 = paste0("picard_downsample ", normal, " ", ds2, " ", normal, ".", ds2, "ds" )
  }
  cm3 = paste0("samtools merge ", output, ".", wanttf, ".bam ", tumor, ".", ds1, "ds.bam ",  normal, ".", ds2, "ds.bam && samtools index ", output, ".", wanttf, ".bam")
  cm4 = paste0("rm ", tumor, ".", ds1, "ds.ba? ",  normal, ".", ds2, "ds.ba?")
  print(paste(cm1, cm2, cm3, cm4, sep=" && "))
  #print(paste(cm1))
}

#downsample to certain coverage
dsbycov <- function(sample, otdep, dsdep, dupaware) {
  ds1 = dsdep / otdep
  stopifnot(ds1 <= 1)
  ndigit = 5
  dsdep = format(round(dsdep, 4), nsmall=4, scientific = F)
  ds1 = format(round(ds1, ndigit), nsmall=ndigit, scientific = F)
  if (dupaware)
    cm1 = paste0("picard_downsample ", sample, " ", ds1, " ", sample, ".", dsdep, " true" )
  else
    cm1 = paste0("picard_downsample ", sample, " ", ds1, " ", sample, ".", dsdep)
  print(cm1)
}

out = "ct19_MSI_tumor.stdngs.insilico.dilu"
norm7=5.779632
norm8=6.220596
tumor7=6.275061
tumor8=6.713597
normal_bam = "8_MSI_normal_HGCTGCCX2.8.aligned.duplicates_marked.bam"
tumor_bam = "8_MSI_tumor_HGCTGCCX2.8.aligned.duplicates_marked.bam"

norm2=2.02
tumor1=2.21
normal_bam="ct19_normal2.merged.cds_consensus.mol_consensus.aligned.bam"
tumor_bam="ct19_tumor1.merged.cds_consensus.mol_consensus.aligned.bam"

out = "ct19_tumor2.merged.cds_consensus.mol_consensus.aligned.bam"
sink("ds.run")
#series = c(0.66, 0.22, 0.11, 0.06, 0.03, 0.01, 0.006, 0.003, 0.001, 0.0005, 0.0003, 0.0001, 0)
series = c(1, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001, 0.0005, 0.00025, 0.0001)
for (x in series) {
  dsbycov(out, 1.57, x, F)
}
sink()

###Mutation signature simulation
tumor_cov= 0.64 
normal_cov = 5.1
target_cov = 2
tumor="13-383_ONC147579.cds_consensus.mol_consensus.aligned.bam"
normal="CDS_V2.merged.cds_consensus.mol_consensus.aligned.bam"
series = c(0.2, 0.1, 0.05, 0.01, 0.005, 0.003, 0.001, 0.0005, 0.0003, 0.0001, 0)
sink("ds.run")
for (x in series) {
  dsbytf(x, target_cov, tumor_cov, normal_cov, 0.66, tumor, normal, "ct19_tumor2_into_normal_2x", F)
}
sink()
