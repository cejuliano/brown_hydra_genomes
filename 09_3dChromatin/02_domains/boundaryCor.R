library(rstudioapi)
library(ggplot2)

setwd(dirname(getActiveDocumentContext()$path))

#import information on the nearest HIC boundary for each gene
tadLink <- read.delim('genesCloseTads.bed', header = F)

#drop any genes that fall within a boundary
tadLink <- tadLink[tadLink$V13 != 0,]

#save original df to use later
tadLink.orig <- tadLink

#split genes by their nearest TAD
tadLink <- split(tadLink,tadLink$V10)


crossPairs <- lapply(tadLink, function(x){
  #get all genes that lie downstream of the domain boundary
  rightG <- x[x$V13 < 0,]               #weirdly, bedtools gave downstream genes negative distance values
  
  #pick the downstream gene that is closest to the domain boundary
  rightG <- rightG[which.max(rightG$V13),'V4']
  
  #get numeric equivalent of gene ID (to check if left and right genes are consecutive)
  rightG.n <- as.numeric(gsub('HVAEP1_G','',rightG))
  
  #get all genes that lie upstream of the domain boundary 
  leftG <- x[x$V13 > 0,]
  
  #pick the upstream gene that is closest to the domain boundary
  leftG <- leftG[which.min(leftG$V13),'V4']
  
  #get numeric equivalent of gene ID
  leftG.n <- as.numeric(gsub('HVAEP1_G','',leftG))
  
  #if the two genes flanking a boundary are consecutive, return the gene pair
  #otherwise do nothing
  if(length(rightG.n) > 0 & length(leftG.n) > 0){
    if((rightG.n - leftG.n) == 1){
      return(c(leftG,rightG))
    }
  }
})

#drop empty results
crossPairs <- crossPairs[sapply(crossPairs,length) > 0]

#collapse into table
crossPairs <- do.call(rbind,crossPairs)

#import NMF gene scores
gScore <- read.delim('../../ds/nmf/final/whole_unfilt_fine_broad.gene_spectra_score.k_28.dt_0_2.txt',row.names = 1)

gScore <- t(gScore)

#fix gene name formatting
rownames(gScore) <- gsub('[.]','_',rownames(gScore))

#drop gene pairs that don't have gene scores
crossPairs <- as.data.frame(crossPairs[crossPairs[,1] %in% rownames(gScore) &
                           crossPairs[,2] %in% rownames(gScore),])

#generate gene ID for genes that are downstream of the righthand genes in the crosspairs list
#these genes will be in the same domain as the righthand crosspairs gene
cisPairs <- as.numeric(substr(crossPairs[,2],9,14)) + 1

cisPairs <- formatC(cisPairs, width = 6, format = "d", flag = "0")

cisPairs <- paste0('HVAEP1_G',cisPairs)

cisPairs <- data.frame(g1 = crossPairs[,2],g2=cisPairs)

#make sure cispairs both have gscores
cisPairs <- as.data.frame(cisPairs[cisPairs[,1] %in% rownames(gScore) &
                                     cisPairs[,2] %in% rownames(gScore),])

#check and make sure the two cispair genes are indeed in the same domain
cisPairs <- cisPairs[cisPairs$g2 %in% tadLink.orig$V4,]

cisPairs$g1Bid <- tadLink.orig[match(cisPairs$g1, tadLink.orig$V4),'V10']
cisPairs$g2Bid <- tadLink.orig[match(cisPairs$g2, tadLink.orig$V4),'V10']

cisPairs <- cisPairs[cisPairs$g1Bid == cisPairs$g2Bid,]

#limit the crosspairs to only those genes that also had a valid cis pair
crossPairs <- crossPairs[crossPairs$V2 %in% cisPairs$g1,]

#compare gene scores across metagenes for crosspairs
xPairCor <- apply(crossPairs,1,function(x){
  cor(gScore[x[1],],gScore[x[2],])
})

#compare gene scores across metagenes for cispairs
cisPairCor <- apply(cisPairs,1,function(x){
  cor(gScore[x[1],],gScore[x[2],])
})

#generate dataframe of correlation scores for plotting
plotDF <- data.frame(corVal = c(cisPairCor,xPairCor), lab = as.factor(rep(c('cis','cross'),c(nrow(cisPairs),nrow(cisPairs)))))

#define cross as first level of factors (specifying plotting order)
plotDF$lab <- relevel(plotDF$lab, "cross")

ggplot(plotDF,aes(x=lab,y=corVal,fill=lab)) + geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.15) + theme_bw()
ggsave('tadExpCor.pdf',width = 4, height = 6)

#use t-test to see if cispairs have significantly higher correlation scores than crosspairs
t.test(x=cisPairCor,y=xPairCor,alternative = 't',var.equal = F)

#look at domain sizes real quick
boundaries <- read.delim('../aep16k_domains.bed',header = F)

boundaries$size <- boundaries$V3 - boundaries$V2

hist(boundaries$size)

median(boundaries$size)   #median size of 176 kb
