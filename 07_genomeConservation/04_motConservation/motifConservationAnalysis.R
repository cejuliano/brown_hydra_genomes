library(rstudioapi)
library(ggplot2)
library(gridExtra)
library(plyr)
library(plotly)

setwd(dirname(getActiveDocumentContext()$path))

#aep motifs that were in an ATAC peak and didn't intersect coding sequence
allA <- read.delim('motifDB/aepMotifsATACnCDS.sort.bed',header = F)

#suffled equivalent
allA.S <- read.delim('motifDB/aepShufMotifsATACnCDS.sort.bed',header = F)

#get file list of conserved motif bed files
mFiles <- list.files(path='motifDB/', pattern = 'OlapCon.bed',recursive = T, full.names = T)

#shuffled subset
smFiles <- mFiles[grepl('Shuf',mFiles)]

#non-shuffled subset
mFiles <- mFiles[!grepl('Shuf',mFiles)]

#import motif hits
realMots <- lapply(mFiles,read.delim,header=F)

shufMots <- lapply(smFiles,read.delim,header=F)


#generate list of unique IDs for each motif hit
realMotID <- lapply(realMots,function(x){
  apply(x[,1:4],1,paste,collapse='_')
})

shufMotID <- lapply(shufMots,function(x){
  apply(x[,1:4],1,paste,collapse='_')
})

#get motifs that are conserved from AEP to 105 and as well as in either viridissima or oligactis
realMotIDCon <- which(realMotID[[1]] %in% unlist(realMotID[2:3]))

shufMotIDCon <- which(shufMotID[[1]] %in% unlist(shufMotID[2:3]))

#make table of conserved motifs
realMotCon <- realMots[[1]][realMotIDCon,]

shufMotCon <- shufMots[[1]][shufMotIDCon,]

#get conservation frequencies for each motif 
realMotConTab <- as.data.frame(table(realMotCon$V4))

shufMotConTab <- as.data.frame(table(shufMotCon$V4))

#get frequencies for all motif instances in the AEP genome (not just conserved)
realMotTab <- as.data.frame(table(allA$V4))

shufMotTab <- as.data.frame(table(allA.S$V4))

#drop motifs that aren't in the other motifs
#some motifs had no hits and some shuffled motif had no hits
#so the two lists aren't perfectly identical
motUse <- union(shufMotTab$Var1,realMotTab$Var1)

#combine all tables
#or each motif it summarizes:
#the number of conserved sites in the AEP genome
#the number of conserved sites in the AEP genome for the shuffled version
#the total number of sites in the AEP genome
#the total number of in the AEP genome for the shuffled version 
allTab <- lapply(list(realMotConTab,shufMotConTab,realMotTab,shufMotTab), function(x){
  x[match(motUse,x[,1]),2]
})

allTab <- as.data.frame(do.call(cbind,allTab))

colnames(allTab) <- c('conReal','conShuf','allReal','allShuf')

rownames(allTab) <- motUse

#calculate the odds that an instance of the genuine motif is conserved
allTab$realOds <- allTab$conReal/(allTab$allReal - allTab$conReal)

#calculate the odds that an instance of the shuffled motif is conserved
allTab$shufOdds <- allTab$conShuf/(allTab$allShuf - allTab$conShuf)

#compare the odds and calculate the log odds ratio
allTab$lor <- log(allTab$realOds/allTab$shufOdds)

allTab <- allTab[complete.cases(allTab),]

#load motif names/information
motInfo <- read.csv('MotifDB/motifInfo.csv', row.names = 1)

allTab$name <- mapvalues(rownames(allTab),from=motInfo$ID,to=motInfo$name,warn_missing = F)

allTab$fam <- mapvalues(rownames(allTab),from=motInfo$ID,to=motInfo$family,warn_missing = F)

#order by log odds ratio, with the highest values being ranked first
allTab <- allTab[order(-allTab$lor),]

#index column to encode rank of each motif
allTab$idx <- 1:nrow(allTab)

#function to perform chi-square test to determine if conservation frequency 
#is different between the shuffled and real motif
doChiTest <- function(x){
  conTab <- matrix(as.numeric(x[1:4]),nrow = 2, byrow = T)
  chiRes <- chisq.test(conTab)
  return(chiRes$p.value)
}

#test for enrichment depletion
allTab$pval <- apply(allTab,1,doChiTest)

#adjust for multiple testing
allTab$fdr <- p.adjust(allTab$pval,method = 'fdr')

#classify motifs as enriched, neutral, or depleted
allTab$res <- 'neutral'

allTab[allTab$lor > 0 & allTab$fdr <= 0.01,'res'] <- 'enriched'

allTab[allTab$lor < 0 & allTab$fdr <= 0.01,'res'] <- 'depleted'

#generate plot of log odds ratio values for motifs
ggplot(allTab,aes(x=idx,y=lor,color=res,name=name)) + geom_point() + theme_bw()
ggsave('conMotPlot.pdf',width=6,height=5)

#highlight particular motifs of interest in the plot
highlightList <- c('NEUROD1','FOXA2','PAX1','MYC','POU4F3','TCF7L1','Smad2::Smad3','Su(H)','SOX4','TCF7')
highlightList <- highlightList[match(allTab$name,highlightList)]
highlightList <- highlightList[!is.na(highlightList)]

allTab$highlight <- allTab$name %in% highlightList

allTab[allTab$name %in% highlightList,'res'] <- 'highlight'

allTab <- allTab[order(allTab$highlight),]

allTab$highlight[allTab$highlight] <- 0.6
allTab$highlight[allTab$highlight ==0] <- 0.5

ggplot() + geom_point(data = allTab,aes(x=idx,y=lor,color=res,name=name,size=highlight)) + theme_bw() + scale_size(range = c(2,3))
ggsave('conMotPlotHighlight.pdf',width=6,height=5)


#export conservation results
write.csv(allTab,'motifConservationStats.csv',row.names = T)

#export list of conserved motif instances in the AEP genome
write.table(realMotCon[,1:6],'conMotsATAC.bed',col.names = F,row.names = F,sep = '\t',quote = F)
write.table(shufMotCon[,1:6],'conShufMotsATAC.bed',col.names = F,row.names = F,sep = '\t',quote = F)

