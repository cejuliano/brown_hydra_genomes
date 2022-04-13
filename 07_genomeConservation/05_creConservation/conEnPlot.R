library(rstudioapi)
library(ggplot2)
setwd(dirname(getActiveDocumentContext()$path))

#import list of conserved H3K4me1 peaks
#this file includes how far each peak is from the nearest TSS
enAnnot.h41 <- read.delim('uropaOut/conPeaksH41_finalhits.txt')

#import list of conserved H3K4me1 peaks that did not intersect with a H3K4me3 peak
enAnnot.h41.nonh43 <- read.delim('h41.conPeaks.en.bed',header=F)

#drop H3K4me1 peaks that overlapped H3K4me3 peaks (likely core promoter regions)
enAnnot.h41 <- enAnnot.h41[enAnnot.h41$peak_id %in% enAnnot.h41.nonh43$V4,]

#only focus on peaks that did not fall within a gene
enAnnot.h41 <- enAnnot.h41[enAnnot.h41$relative_location %in% c('Downstream','Upstream'),]

#for upstream genes, set distance to negative
enAnnot.h41[enAnnot.h41$relative_location == 'Upstream','distance'] <- -enAnnot.h41[enAnnot.h41$relative_location == 'Upstream','distance']

#ggplot(enAnnot.h41,aes(distance)) + geom_density() + geom_vline(xintercept = -2000) + xlim(-5e4,5e4)


#repeat the same basic process for ATAC-seq peaks
enAnnot.ATAC <- read.delim('uropaOut/conPeaksATAC_finalhits.txt')

enAnnot.ATAC.nonh43 <- read.delim('ATAC.conPeaks.en.bed',header=F)

enAnnot.ATAC <- enAnnot.ATAC[enAnnot.ATAC$peak_id %in% enAnnot.ATAC.nonh43$V4,]

enAnnot.ATAC <- enAnnot.ATAC[enAnnot.ATAC$relative_location %in% c('Downstream','Upstream'),]

enAnnot.ATAC[enAnnot.ATAC$relative_location == 'Upstream','distance'] <- -enAnnot.ATAC[enAnnot.ATAC$relative_location == 'Upstream','distance']

#ggplot(enAnnot.ATAC,aes(distance)) + geom_density() + geom_vline(xintercept = -2000) + xlim(-5e4,5e4)

#combine H41 and ATAC data into a single table for plotting
enAnnot.h41$mark <- 'H41'

enAnnot.ATAC$mark <- 'ATAC'

enAnnot.both <- rbind(enAnnot.h41,enAnnot.ATAC)

ggplot(enAnnot.both,aes(distance,color=mark,fill=mark)) + 
  geom_density(alpha=0.2) + 
  geom_vline(xintercept = -10000,linetype='longdash') + 
  geom_vline(xintercept = 10000,linetype='longdash') + 
  xlim(-5e4,5e4) + 
  theme_bw() + facet_wrap(.~mark)
ggsave('enDistribution.pdf',width = 9,height=4)

#basic stats on conserved distal peaks
summary(abs(enAnnot.ATAC$distance))

length(which(abs(enAnnot.ATAC$distance) > 1e4))
length(which(abs(enAnnot.ATAC$distance) > 1e4))/nrow(enAnnot.ATAC)

length(which(abs(enAnnot.h41$distance) > 1e4))
length(which(abs(enAnnot.h41$distance) > 1e4))/nrow(enAnnot.h41)

