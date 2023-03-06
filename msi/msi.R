library("rjson")
library(data.table)
library(tidyverse)

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
  profile = fromJSON(file=str_profile)
  df = fread(sample_file, header=F, col.names = c("chrom", "pos", "rept", "length" ,"tumor", "normal", "prefix", "suffix", "popAF", "germline"))
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
  for (f in list.files(dirin)) {
    fields = unlist(strsplit(f, ".", fixed = TRUE))
    tf = as.numeric(paste0("0.", fields[2])) * 100
    sample = fields[3]
    sid = c(sid, paste(sample, tf))
  }
  result = data.frame("sid" = unique(sid), "all_filtered" = 0, "msi_filtered" = 0, score=0)
  for (f in list.files(dirin)) {
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
      df = fread(file.path(dirin, f), header=F, col.names = c("chrom", "pos", "rept", "length" ,"tumor", "normal", "prefix", "suffix", "popAF", "germline"))
    } else {
      df = fread(file.path(dirin, f), header=F, col.names = c("chrom", "pos", "rept", "length" ,"tumor", "normal", "prefix", "suffix", "popAF"))
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
  result = result %>% separate(sid, c("sample", "tf"), sep=" ") %>% pivot_longer(c("msi_filtered"), names_to="Germ. Filtering")
  result$tf = as.numeric(result$tf)
  result
}

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

msmutect_parse <- function(msdir) {
  tfs = c()
  score = c()
  for (f in list.files(msdir, "mut.tsv")) {
    fields = unlist(str_split(f, "\\."))
    tf = as.numeric(fields[2])
    tfs = c(tfs, tf)
    msdf = fread(file.path(msdir, f))
    msdf = msdf %>% filter(CALL == 1)
    score = c(score, nrow(msdf))
  }
  msmutect_result = data.frame("score" = score, tf=tfs/100)
  msmutectdf = post_process(msmutect_result, "Standard_MSMutect")
  msmutectdf
}

###LOAD ALL DATA!!!!!!!!!!!!
########################
setwd("paper/data/msi")
str_profile ="msi_profile.json"
codec2x = MsiDetect("codec", str_profile)
codecplot2x = post_process(codec2x, "CODEC")
msmutectdf1 = msmutect_parse("msmutect")

plotdf= rbind(codecplot2x %>% select(score, tf, class, observed), msmutectdf1)
ggplot(data=plotdf, aes(x=as.factor(tf), y=score, color = class)) + 
  geom_line(aes(group=class)) +
  geom_point(size=3) + 
  theme_classic(24) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_color_manual(name="", values = c("deepskyblue3", "grey"), labels = c("CODEC", "Standard_MSMutect")) +
  scale_x_discrete(labels=unique(plotdf$tf)) +
  scale_y_log10(breaks = c(1e-1, 1, 10, 100, 1000, 1e4), labels = trans_format("log10", math_format(10^.x)), limits=c(1e-1, 1e4)) +
  annotation_logticks(sides = "l", outside = F) +  coord_cartesian(clip = "off", ylim = c(1e-1, 1e4)) + 
  xlab("") + ylab("MSIscore") +
  theme(legend.position="bottom", axis.text = element_text(size = 22, color = "black"))
