library(rstudioapi)
library(plyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(viridis)

setwd(dirname(getActiveDocumentContext()$path))

#import clytia enrichment
enScores <- readRDS('../cl/remap/enrichment/clEnScores.rds')

clEn.DF <- read.csv('../cl/remap/enrichment/clEnrichmentRes.csv')

#load hydra motif enrichment scores
hvEn <- readRDS('../nmf/motScores.rds')

#disregard non-conserved motifs
motifCon <- read.csv('../../alignment_conservation/motifConservationStats.csv',row.names = 1)

motifCon.keep <- rownames(motifCon[motifCon$res == 'enriched',])

hvEn <- hvEn[,colnames(hvEn) %in% motifCon.keep]

#look for similar enrichment patterns using psuedocells
hvCl <- readRDS('aepClInt.rds')

hvCl <- FindClusters(hvCl, resolution = 10, graph.name = 'integrated_snn')

DefaultAssay(hvCl) <- 'SCT'

DimPlot(hvCl, group.by = 'integrated_snn_res.10') + NoLegend() + NoAxes()

#export version of hydra enrichment table with pseudocell info for cells
hvEn.exp <- hvEn

hvEn.exp$pseudoCell <- mapvalues(rownames(hvEn), from = colnames(hvCl), to = hvCl$integrated_snn_res.10,warn_missing = F)

#drop any cells not grouped into a psuedocell in the enrichment results matrix
hvEn.exp <- hvEn.exp[grepl('^\\d',hvEn.exp$pseudoCell),]

#drop non-conserved/non-enriched motifs
hvEn <- hvEn[,colnames(hvEn) %in% clEn.DF$motif_ID]

#use motif name instead of motif ID
colnames(hvEn) <- mapvalues(colnames(hvEn),from = clEn.DF$motif_ID,to = clEn.DF$motif_alt_ID,warn_missing = F)

#export version of clytia enrichment table with pseudocell info for cells
enScores.pc <- as.data.frame(enScores[rownames(enScores) %in% colnames(hvCl),])
enScores.pc.id <- mapvalues(rownames(enScores.pc),from=colnames(hvCl),to=hvCl$integrated_snn_res.10,warn_missing = F)

enScores.pc.exp <- as.data.frame(enScores.pc)
enScores.pc.exp$pseudoCell <- enScores.pc.id

#calculate average enrichment score for each motif for each pseudo-cell
#first for clytia
enScores.pc <- split(enScores.pc,enScores.pc.id)
enScores.pc <- lapply(enScores.pc,function(x) apply(x,2,mean))

enScores.pc <- do.call(rbind,enScores.pc)

#then for hydra
hvEn.pc <- as.data.frame(hvEn[rownames(hvEn) %in% colnames(hvCl),])
hvEn.pc.id <- mapvalues(rownames(hvEn.pc),from=colnames(hvCl),to=hvCl$integrated_snn_res.10,warn_missing = F)

hvEn.pc <- split(hvEn.pc,hvEn.pc.id)
hvEn.pc <- lapply(hvEn.pc,function(x) apply(x,2,mean))

hvEn.pc <- do.call(rbind,hvEn.pc)

#restrict psuedo-cells and motifs to only those that are present in both species
hvEn.pc <- hvEn.pc[,colnames(hvEn.pc) %in% colnames(enScores.pc)]
enScores.pc <- enScores.pc[,colnames(enScores.pc) %in% colnames(hvEn.pc)]

hvEn.pc <- hvEn.pc[rownames(hvEn.pc) %in% rownames(enScores.pc),]
enScores.pc <- enScores.pc[rownames(enScores.pc) %in% rownames(hvEn.pc),]

hvEn.pc <- hvEn.pc[,colnames(enScores.pc)]

write.csv(hvEn.pc,'pcAveHvEn.csv')
write.csv(enScores.pc,'pcAveClEn.csv')

#calculate correlation scores across pseudo-cells for each motif
crossMotCor <- vapply(1:ncol(hvEn.pc),function(x){
  cor(hvEn.pc[,x], enScores.pc[,x], method = 'pearson')
},numeric(1))

crossMotCor <- data.frame(motID = colnames(enScores.pc),crossMotCor)

#reduce redundancy of cor scores
lrMots <- read.csv('../nmf/hydraEnrichmentMatrixLR.csv')[,1]

lrMots <- gsub('.* ','',lrMots)

crossMotCor <- crossMotCor[crossMotCor$motID %in% lrMots,]

#subset to look only at high scoring enrichment results
crossMotCor.match <- crossMotCor[crossMotCor$crossMotCor > 0.5,]

#export correlation results
write.csv(crossMotCor,'crossMotCor.csv',row.names = F)

#function for plotting motif enrichment scores for each species using the aligned seurat object
motScorePlot <- function(x){
  newMD <- hvCl@meta.data
  
  newMD$cID <- rownames(newMD)
  
  motCheck <- x
  
  motScoreDF <- data.frame(cID = c(rownames(enScores),rownames(hvEn)), motScore = c(enScores[,motCheck],hvEn[,motCheck]))
  
  newMD <- merge(newMD,motScoreDF,by='cID',all.x=T)
  
  newMD[is.na(newMD$motScore),'motScore'] <- 0
  
  rownames(newMD) <- newMD$cID
  
  newMD$cID <- NULL
  
  newMD <- newMD[colnames(hvCl),]
  
  hvCl@meta.data <- newMD
  
  gg <- FeaturePlot(hvCl,'motScore',split.by = 'species',pt.size = 0.6,order = T,combine=F)
  gg[[1]] + theme_void() + theme(plot.title = element_blank())
  ggsave(paste0('enPlot_',x,'_aep.png'),width = 5.5,height=5,dpi=450)
  
  gg[[2]] + theme_void() + theme(plot.title = element_blank())
  ggsave(paste0('enPlot_',x,'_cl.png'),width = 5.5,height=5,dpi=450)
}


motScorePlot('EBF3')

gg <- FeaturePlot(hvCl,'HVAEP1-G005780',split.by = 'species',order = T,pt.size=0.6,combine=F)
gg[[1]] + theme_void() + theme(plot.title = element_blank())
ggsave('ebfTran_aep.png',width = 5.5,height=5,dpi=450)

gg[[2]] + theme_void() + theme(plot.title = element_blank())
ggsave('ebfTran_cl.png',width = 5.5,height=5,dpi=450)


motScorePlot('PAX6')

gg <- FeaturePlot(hvCl,'HVAEP1-G008423',split.by = 'species',order = T,pt.size=0.6,combine=F)
gg[[1]] + theme_void() + theme(plot.title = element_blank())
ggsave('paxTran_aep.png',width = 5.5,height=5,dpi=450)

gg[[2]] + theme_void() + theme(plot.title = element_blank())
ggsave('paxTran_cl.png',width = 5.5,height=5,dpi=450)



motScorePlot('POU4F1')

gg <- FeaturePlot(hvCl,'HVAEP1-G023106',split.by = 'species',order = T,pt.size=0.6,combine=F)
gg[[1]] + theme_void() + theme(plot.title = element_blank())
ggsave('pouTran_aep.png',width = 5.5,height=5,dpi=450)

gg[[2]] + theme_void() + theme(plot.title = element_blank())
ggsave('pouTran_cl.png',width = 5.5,height=5,dpi=450)


motScorePlot('FOXN3')

gg <- FeaturePlot(hvCl,'HVAEP1-G005557',split.by = 'species',order = T,pt.size=0.6,combine=F)
gg[[1]] + theme_void() + theme(plot.title = element_blank())
ggsave('atohTran_aep.png',width = 5.5,height=5,dpi=450)

gg[[2]] + theme_void() + theme(plot.title = element_blank())
ggsave('atohTran_cl.png',width = 5.5,height=5,dpi=450)


motScorePlot('ATOH7')

gg <- FeaturePlot(hvCl,'HVAEP1-G021588',split.by = 'species',order = T,pt.size=0.6,combine=F)
gg[[1]] + theme_void() + theme(plot.title = element_blank())
ggsave('atohTran_aep.png',width = 5.5,height=5,dpi=450)

gg[[2]] + theme_void() + theme(plot.title = element_blank())
ggsave('atohTran_cl.png',width = 5.5,height=5,dpi=450)

