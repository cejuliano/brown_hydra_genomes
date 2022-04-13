library(rstudioapi)
library(ggplot2)
library(dplyr)

setwd(dirname(getActiveDocumentContext()$path))

getConPeaks <- function(markName){
  
  #import conservation score matrix from deeptools
  valueMat <- paste0('geneMatrix.',markName,'.specCon.names.txt')
  
  #(need to import header separately because of formatting problems)
  conMat <- read.delim(valueMat,skip = 3, sep = '\t',header = F)
  
  conMatCol <- t(read.delim(valueMat,skip = 2, nrows = 1, sep = '\t',header = F))[-1]
  
  colnames(conMat) <- conMatCol
  
  #this second matrix includes peak names and coordinates
  #which we'll need for exporting the results
  #again nead to import the header separately
  nameMat <- paste0('geneMatrix.',markName,'.specCon.regions.txt')
  
  conMatRow <- read.delim(nameMat, sep = '\t')
  
  rownames(conMat) <- conMatRow$name
  
  #split the conservation matrix by species (1100 columns per species)
  conMat <- lapply(list(c(1,1100),c(1101,2200),c(2201,3300),c(3301,4400)),function(x){
    return(conMat[,x[1]:x[2]])
  })
  
  #each matrix includes flanking regions we don't need,
  #so we just drop those
  conMat.peak <- lapply(conMat,function(x) x[,501:600])
  
  #calculate the average conservation across each peak
  conMat.peak.cnsvd <- lapply(conMat.peak,function(x) apply(x,1,mean))
  
  #combine the average scores from each species into a single table
  conMat.peak.cnsvd <- do.call(cbind,conMat.peak.cnsvd)
  
  #make a plot of conservation score distribution
  #general basis for classification approach
  conMat.peak.cnsvd.plot <- data.frame(conVal = as.numeric(conMat.peak.cnsvd),
                                       spec = rep(c('105','olig','virid','clytia'),rep(nrow(conMat.peak.cnsvd),4)))
  
  conMat.peak.cnsvd.plot$spec <- factor(conMat.peak.cnsvd.plot$spec, levels = c("105", "olig", "virid","clytia"))
  
  conMat.peak.cnsvd.plot <- conMat.peak.cnsvd.plot[conMat.peak.cnsvd.plot$conVal > 0,]
  
  print(ggplot(conMat.peak.cnsvd.plot,aes(x=conVal,col=spec)) + 
    geom_density() + 
    #scale_y_continuous(trans='log2') + 
    facet_wrap(.~spec))
  ggsave(paste0(markName,'.conDist.png'),width=6,height=5,dpi=300)
  
  #perform K-means clustering to partition scores into two populations
  #then return only the members of the higher score population
  conMat.clust <- apply(conMat.peak.cnsvd,2,function(x){
    kRes <- kmeans(x,2)
    cUse <- which.max(kRes$centers)
    return(kRes$cluster == cUse)
  })
  
  #subset our peak list to include only those peaks that were in the 'conserved'
  #population in at least two pairwise comparisons
  conMat.clust <- cbind(as.data.frame(conMat.clust),con=(rowSums(conMat.clust) > 1))
  
  conMat.clust$mark <- markName
  
  conMat.clust.pos <- conMat.clust[conMat.clust[,5],]
  
  conMat.peak.cnsvd.exp <- conMatRow[conMatRow$name %in% rownames(conMat.clust.pos),]
  
  outName <- paste0(markName,'.conPeaks.bed')
  
  write.table(conMat.peak.cnsvd.exp[,1:6],outName,sep = '\t',quote = F, col.names = F, row.names = F)
  
  return(conMat.clust)
}

conRes <- lapply(c('CDS','H41','H43','H273','ATAC'),function(x) getConPeaks(x))
