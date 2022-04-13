setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ggplot2)

#import STAR mapping logs for AEP data mapped to the 105 reference
a2mRNA <- list.files(path='rnaMapStats/a2m',recursive = T,full.names = T)

a2mRNA <- lapply(a2mRNA,read.delim)

#pull total read number
a2mRNA.total <- sum(sapply(a2mRNA,function(x){
  as.numeric(x[4,2])
}))

#calculate portion of reads that were unmapped
a2mRNA.unmapped <- sum(sapply(a2mRNA,function(x){
  sum(as.numeric(x[c(24,27,29,31),2]))
}))/a2mRNA.total

#calculate portion of reads that were uniquely mapped
a2mRNA.unique <- sum(sapply(a2mRNA,function(x){
  as.numeric(x[7,2])
}))/a2mRNA.total

#calculate portion of reads that were multi-mappers
a2mRNA.mm <- sum(sapply(a2mRNA,function(x){
  as.numeric(x[22,2])
}))/a2mRNA.total

#import STAR mapping logs for AEP data mapped to the AEP reference
a2aRNA <- list.files(path='rnaMapStats/a2a',recursive = T,full.names = T)

a2aRNA <- lapply(a2aRNA,read.delim)

#pull total read number
a2aRNA.total <- sum(sapply(a2aRNA,function(x){
  as.numeric(x[4,2])
}))

#calculate portion of reads that were unmapped
a2aRNA.unmapped <- sum(sapply(a2aRNA,function(x){
  sum(as.numeric(x[c(24,27,29,31),2]))
}))/a2aRNA.total

#calculate portion of reads that were uniquely mapped
a2aRNA.unique <- sum(sapply(a2aRNA,function(x){
  as.numeric(x[7,2])
}))/a2aRNA.total

#calculate portion of reads that were multi-mappers
a2aRNA.mm <- sum(sapply(a2aRNA,function(x){
  as.numeric(x[22,2])
}))/a2aRNA.total

#fold increase in unmapped reads with 105 reference
a2mRNA.unmapped/a2aRNA.unmapped


#import error logs for AEP data mapped to the 105 reference
a2mATAC <- list.files(path='atacMapStats/a2m',recursive = T,full.names = T)

a2mATAC <- lapply(a2mATAC,read.delim)

#isolate bowtie output from error logs
a2mATAC <- lapply(a2mATAC,function(x) {
  retVec <- x[c(97:99),]
  retVec <- as.numeric(gsub(' \\(.*','',retVec))
  return(retVec)
})

#get total read number
a2mATAC.total <- sum(sapply(a2mATAC,function(x) sum(x)))

#calculate portion of reads that were unmapped
a2mATAC.unmapped <- sum(sapply(a2mATAC,function(x) x[1]))/a2mATAC.total

#calculate portion of reads that were uniquely mapped
a2mATAC.unique <- sum(sapply(a2mATAC,function(x) x[2]))/a2mATAC.total

#calculate portion of reads that were multi-mappers
a2mATAC.mm <- sum(sapply(a2mATAC,function(x) x[3]))/a2mATAC.total

#import error logs for AEP data mapped to the AEP reference
a2aATAC <- list.files(path='atacMapStats/a2a',recursive = T,full.names = T)

a2aATAC <- lapply(a2aATAC,read.delim)

#isolate bowtie output from error logs
a2aATAC <- lapply(a2aATAC,function(x) {
  retVec <- x[c(97:99),]
  retVec <- as.numeric(gsub(' \\(.*','',retVec))
  return(retVec)
})

#get total read number
a2aATAC.total <- sum(sapply(a2aATAC,function(x) sum(x)))

#calculate portion of reads that were unmapped
a2aATAC.unmapped <- sum(sapply(a2aATAC,function(x) x[1]))/a2aATAC.total

#calculate portion of reads that were uniquely mapped
a2aATAC.unique <- sum(sapply(a2aATAC,function(x) x[2]))/a2aATAC.total

#calculate portion of reads that were multi-mappers
a2aATAC.mm <- sum(sapply(a2aATAC,function(x) x[3]))/a2aATAC.total

#plot to compare ratios for the AEP and 105 reference
plot.df <- data.frame(percent=c(a2aRNA.unique,a2aRNA.mm,a2aRNA.unmapped,
                                a2mRNA.unique,a2mRNA.mm,a2mRNA.unmapped,
                                a2aATAC.unique,a2aATAC.mm,a2aATAC.unmapped,
                                a2mATAC.unique,a2mATAC.mm,a2mATAC.unmapped),
                      label=rep(c("unique","mm","unmapped"),4),
                      mapping=c(rep('a2aRNA',3),rep('a2mRNA',3),rep('a2aATAC',3),rep('a2mATAC',3)))

plot.df$label <- factor(plot.df$label, levels = c("unmapped", "mm",'unique'))

plot.df$mapping <- factor(plot.df$mapping, levels = c('a2aRNA','a2mRNA','a2aATAC','a2mATAC'))

ggplot(plot.df,aes(x=mapping,y=percent,fill=label)) + geom_col() + theme_bw()

ggsave('mappingStats.pdf',width=6,height=8)

#fold increase in unmapped reads with 105 reference
a2mATAC.unmapped/a2aATAC.unmapped


