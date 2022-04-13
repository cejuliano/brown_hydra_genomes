library(rstudioapi)
library(plyr)

setwd(dirname(getActiveDocumentContext()$path))

#import NMF results to get coexpressed genes

clNMF <- read.delim('../nmf/cl_fine.gene_spectra_score.k_37.dt_0_13.txt',row.names = 1)

#load genome annotation 

clGtf <- read.delim('../../clFinal.merge.longestIso.rename.gtf',header = F)

clGtf <- clGtf[clGtf$V3 == 'gene',]

#reset coordinates to correspond to putative promoter regions (1000 bp upstream to TSS)

newCoords <- clGtf[,c(4,5,7)]

newCoords.p <- newCoords[newCoords$V7 == '+',]

newCoords.p$V5 <- newCoords.p$V4 + 100

newCoords.p$V4 <- newCoords.p$V4 - 1000

newCoords.m <- newCoords[newCoords$V7 == '-',]

newCoords.m$V4 <- newCoords.m$V5 - 100

newCoords.m$V5 <- newCoords.m$V5 + 1000

clGtf[clGtf$V7 == '+',c(4,5)] <- newCoords.p[,1:2]

clGtf[clGtf$V7 == '-',c(4,5)] <- newCoords.m[,1:2]

#convert from 1 based to 0 based coordinates

clGtf$V4 <- clGtf$V4 - 1

#convert to bed-like table

clBed <- clGtf[,c(1,4,5,9,6,7)]

clBed$V9 <- gsub(' ;.*','',clBed$V9)

clBed$V9 <- gsub('gene_id ','',clBed$V9)

clBed$V4[clBed$V4 < 0] <- 0

#import chrom sizes
cSizes <- read.delim('../../../cl/clytiaG.fa.fai', header = F)

#truncate any coordinates that go beyond the chrom boundaries
clBed$max <- as.numeric(mapvalues(clBed$V1, from = cSizes$V1, to = cSizes$V2, warn_missing = F))

clBed$V5[clBed$V5 > clBed$max] <- clBed[clBed$V5 > clBed$max,'max']

clBed$max <- NULL

#generate versions of bed file for each metagene with a score for each gene 

clNMF <- as.data.frame(t(clNMF))

clNMF$ID <- rownames(clNMF)

clNMF$ID <- gsub('[.]','_',clNMF$ID)

unlink('mgBeds',recursive = T)
dir.create('mgBeds',showWarnings = F)

options(scipen = 999)

lapply(1:37,function(x){
  clBed.mgs <- clBed
  
  #drop genes without a score
  clBed.mgs <- clBed.mgs[clBed.mgs$V9 %in% clNMF$ID,]
  
  clBed.mgs$score <- as.numeric(mapvalues(clBed.mgs$V9,from=clNMF$ID,to=clNMF[,x],warn_missing = F))
  
  #need to invert score so that genes strongly associated with the mg have small values
  
  clBed.mgs$score <- -clBed.mgs$score
  
  #order by score
  clBed.mgs <- clBed.mgs[order(clBed.mgs$score),]
  
  clBed.mgs$score <- (clBed.mgs$score - min(clBed.mgs$score)) + 1e-6
  
  #combine score with name
  clBed.mgs$V9 <- paste(clBed.mgs$V9,clBed.mgs$score,sep = ' ')
  
  clBed.mgs$score <- NULL
  write.table(clBed.mgs,paste0('mgBeds/','mg',x,'_scoredProms.bed'),row.names = F, col.names = F, sep = '\t', quote = F)
})
