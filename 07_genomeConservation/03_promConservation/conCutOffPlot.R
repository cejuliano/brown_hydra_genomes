library(rstudioapi)
library(ggplot2)
library(mixtools)
library(dplyr)

setwd(dirname(getActiveDocumentContext()$path))

#import conservation matrix (rows are genes, columns are bins/position relative to gene)
#need to import header separately because of parsing problems
conMat <- read.delim('broadConGeneMatrix.names.txt',skip = 3, sep = '\t',header = F)

conMatCol <- t(read.delim('broadConGeneMatrix.names.txt',skip = 2, nrows = 1, sep = '\t',header = F))[-1]

#fix colnames issue
colnames(conMat) <- conMatCol

#calcualte the average conservation for each bin (average conservation profile for all genes)
conMat.ave <- colMeans(conMat,na.rm = T)

conMat.ave <- data.frame(pos=seq(from=0,to=20740,by=10),val=conMat.ave)

#determine the baseline conservation level (average at +/- 10 Kb from genes)
end5 <- mean(conMat.ave[conMat.ave$pos < 1000,'val'])
end3 <- mean(conMat.ave[conMat.ave$pos > 19075,'val'])

#determine the most up and downstream regions from genes that have a conservation values above the baseline
con5 <- max(conMat.ave[conMat.ave$val < end5 & conMat.ave$pos < 10000,'pos'])
con3 <- min(conMat.ave[conMat.ave$val < end3 & conMat.ave$pos > 10000,'pos'])

library(MESS)

#within the region of elevated conservation, calculate the area under the curve
#relative to distance from the TSS or TTS
#that is, how far upstream (for instance) do you have to go from the TSS to encompass
#50% of the total upstream conservation signal

#first look upstream
upAuc <- conMat.ave[conMat.ave$pos >= con5 & conMat.ave$pos <= 10000,]
upAuc <- upAuc[rev(1:nrow(upAuc)),]

upAuc$pos <- abs(upAuc$pos - 10000)

upAuc$auc <- vapply(upAuc$pos, function(x){
  
  auc(x=upAuc$pos,
      y=upAuc$val,
      from=0,
      to = x,
      type = 'spline',
      absolutearea = T)
  
},numeric(1))

#convert auc to percent of total auc
upAuc$auc.p <- upAuc$auc/max(upAuc$auc)

#get position of half auc of upstream chunk
upHalf <- max(upAuc[upAuc$auc.p <= 0.5,'pos'])

#get position of 90% auc of upstream chunk
upNinety <- max(upAuc[upAuc$auc.p <= 0.9,'pos'])

#now look downstream
downAuc <- conMat.ave[conMat.ave$pos <= con3 & conMat.ave$pos >= 10750,]
#downAuc <- downAuc[rev(1:nrow(downAuc)),]

downAuc$pos <- abs(downAuc$pos - 10750)

downAuc$auc <- vapply(downAuc$pos, function(x){
  
  auc(x=downAuc$pos,
      y=downAuc$val,
      from=0,
      to = x,
      type = 'spline',
      absolutearea = T)
  
},numeric(1))

downAuc$auc.p <- downAuc$auc/max(downAuc$auc)

downHalf <- max(downAuc[downAuc$auc.p <= 0.5,'pos'])
downNinety <- max(downAuc[downAuc$auc.p <= 0.9,'pos'])

#generate plot with AUC colored according to chunks calculated above

conMat.ave$fill <- 'null'
#conMat.ave[conMat.ave$pos > con5 & conMat.ave$pos <= 10000,'fill'] <- 'upAll'
conMat.ave[conMat.ave$pos > 10000-upNinety & conMat.ave$pos <= 10000,'fill'] <- 'upNinety'
conMat.ave[conMat.ave$pos > 10000-upHalf & conMat.ave$pos <= 10000,'fill'] <- 'upFifty'
#conMat.ave[conMat.ave$pos < con3 & conMat.ave$pos >= 10750,'fill'] <- 'downAll'
conMat.ave[conMat.ave$pos < 10750+downNinety & conMat.ave$pos >= 10750,'fill'] <- 'downNinety'
conMat.ave[conMat.ave$pos < 10750+downHalf & conMat.ave$pos >= 10750,'fill'] <- 'downFifty'
conMat.ave[conMat.ave$pos > 10000 & conMat.ave$pos < 10750,'fill'] <- 'gene'
conMat.ave[conMat.ave$pos >= con3,'fill'] <- 'null2'

ggplot(conMat.ave,aes(x=pos,y=val)) + 
  theme_bw() + 
  geom_line(color='red',size=2) + 
  geom_vline(xintercept = con5,size=1) + 
  geom_vline(xintercept = con3,size=1) + 
  scale_x_continuous(breaks = c(0,con5,10000,10750,con3,20750),minor_breaks = NULL,labels = c('-10 Kb',paste0('-',10000-con5),'TSS','TTS',paste0('+',10750-con3),'+10 Kb')) + 
  theme(axis.text.x = element_text(face="bold",size=11, angle=315,color='black'),axis.text.y = element_text(face="bold", size=11,color='black')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave('conservationCutoff.pdf',width = 5,height=5)  



ggplot(conMat.ave,aes(x=pos,y=val)) + 
  theme_bw() + 
  geom_area(aes(fill = fill), position = 'identity',alpha=0.8) + 
  geom_line(color='red',size=2) + 
  geom_vline(xintercept = con5,size=1) + 
  geom_vline(xintercept = con3,size=1) + 
  scale_x_continuous(breaks = c(0,con5,10000-upNinety,10000-upHalf,10000,10750,10750 + downHalf,10750 + downNinety,con3,20750),minor_breaks = NULL,
                     labels = c('-10 Kb',paste0('-',10000-con5),-upNinety,-upHalf,'TSS',
                                'TTS',paste0('+',downHalf),paste0('+',downNinety),paste0('+',con3-10750),'+10 Kb')) + 
  theme(axis.text.x = element_text(face="bold",size=11, angle=300,color='black'),axis.text.y = element_text(face="bold", size=11,color='black')) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "none")

ggsave('conservationCutoffFilled.pdf',width = 5,height=5)  


