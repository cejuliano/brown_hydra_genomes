library(rstudioapi)
library(plyr)
library(Seurat)
library(ggplot2)
library(gplots)
library(viridis)

setwd(dirname(getActiveDocumentContext()$path))

#import clytia enrichment
enScores <- readRDS('../cl/remap/enrichment/clEnScoresShuf.rds')

clEn.DF <- read.csv('../cl/remap/enrichment/clEnrichmentResShuf.csv')

#load hydra motif enrichment scores
hvEn <- readRDS('../nmf/motScores_shuf.rds')

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

write.csv(hvEn.pc,'pcAveHvEnShuf.csv')
write.csv(enScores.pc,'pcAveClEnShuf.csv')

#calculate correlation scores across pseudo-cells for each motif
crossMotCor <- vapply(1:ncol(hvEn.pc),function(x){
  cor(hvEn.pc[,x], enScores.pc[,x], method = 'pearson')
},numeric(1))

crossMotCor <- data.frame(motID = colnames(enScores.pc),crossMotCor)

#export correlation results
write.csv(crossMotCor,'crossMotCorShuf.csv',row.names = F)

#subset to look only at high scoring enrichment results
crossMotCor.match <- crossMotCor[crossMotCor$crossMotCor > 0.5,]
print(nrow(crossMotCor.match))
