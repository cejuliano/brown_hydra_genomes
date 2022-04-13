library(Seurat)
library(tidyverse)
library(rstudioapi)
library(glmGamPoi)
library(plotly)
library(RColorBrewer)
library(patchwork)
library(plyr)
library(RColorBrewer)
library(gplots)
library(viridis)
library(fgsea)

#utility to convert transcript ID to gene ID
t2g <- function(x){
  vapply(x, function(y) gsub('HVAEP1_T(\\d+)[.]\\d','HVAEP1-G\\1',y),"")
}

setwd(dirname(getActiveDocumentContext()$path))

####prep gene sets and gene scores for gsea####

#load NMF gene Z scores
gscores <- read.delim('final/whole_unfilt_fine_narrow.gene_spectra_score.k_56.dt_0_13.txt',row.names = 1)

#reformat gene names because cnmf got rid of the underscores
colnames(gscores) <- gsub('[.]','_',colnames(gscores))

#load annotated conserved transcription factor binding sites
mots <- read.delim('../../alignment_conservation/motifAnnot/conMotATAC_finalhits.txt')

#drop motifs that weren't near a gene
mots <- mots[complete.cases(mots),]

#drop all genes that have no conserved motifs
gscores <- gscores[,colnames(gscores) %in% mots$gene_id]

#create a gene set for each binding motif
mots <- split(mots$gene_id,mots$peak_id)

#drop duplicated gene IDs in each gene set
mots <- lapply(mots,unique)

####ID enriched motifs####

#go through the gene scores of each metagene and use them as weights
#to order genes for a gsea analysis
#genes are linked to motifs though the motif gene set we made above
set.seed(12345)
enrichList <- apply(gscores,1,function(x){
  gscores.use <- x
  names(gscores.use) <- colnames(gscores)
  fgseaMultilevel(mots, gscores.use, minSize = 25, nproc = 6,scoreType = 'pos')
})

#Set all normalized enrichment scores (NES) to zero if
#they fail to pass the significance cutoff
#then subset the results table to include just the 
#motif name and enrichment score
enrichList <- lapply(enrichList,function(x){
  x[x$padj > 0.01,'NES'] <- 0
  x[,c(1,6)]
})

#Fix issue where enrichment returned NAs (just set to zero)
enrichList <- lapply(enrichList,function(x){
  x[is.na(x$NES),]$NES <- 0
  return(x)
})

#save motif names column from output to use for rownames later
rowN <- enrichList[[1]][,1]

#extract just the enrichment scores for each metagene
enrichList <- lapply(enrichList,function(x){
  x[,2]
})

#combine enrichment scores into an enrichment matrix spanning all metagenes
enrichDF <- as.data.frame(do.call(cbind,enrichList))

rownames(enrichDF) <- rowN$pathway

colnames(enrichDF) <- paste0('mg',1:ncol(enrichDF))

#save results to save time later
saveRDS(enrichDF,'enrichDFnCDSATAC.rds')
enrichDF <- readRDS('enrichDFnCDSATAC.rds')

####Propagate enrichment scores to seurat object####

#import NMF metagene cell scores
k56.usage <- read.delim('final/whole_unfilt_fine_narrow.usages.k_56.dt_0_13.consensus.txt',row.names = 1)

#normalize cell scores so that that all the scores for a single cell sum to 1
k56.usage <- t(apply(k56.usage,1,function(x) x/sum(x)))

motScores <- apply(k56.usage, 1, function(x) apply(enrichDF,1, function(y) sum(x*y), simplify = T))

motScores <- as.data.frame(t(motScores))

motScores <- motScores[,colSums(motScores) != 0]

#load doublet-free, annotated atlas Seurat object
ds <- readRDS('../dropSeqMapping/nonDubLabeledSeurat.rds')

motScores <- motScores[rownames(ds@meta.data),]

saveRDS(motScores,'motScores.rds')
motScores <- readRDS('motScores.rds')

####motif enrichment heatmap####

#load metagene descriptions
mgD <- read.csv('k56_mg_annot.csv',header = F)

mgD$V1 <- paste0('mg',mgD$V1)

#initialize plotting data object
nesPlot <- enrichDF

#rename the columns to have the more descriptive metagene name
colnames(nesPlot) <- mgD$V2

#import statistics on motif conservation
#we only want to look at motifs that showed
#evidence of conservation
motifCon <- read.csv('../../alignment_conservation/motifConservationStats.csv',row.names = 1)

motifCon.keep <- rownames(motifCon[motifCon$res == 'enriched',])

nesPlot <- nesPlot[rownames(nesPlot) %in% motifCon.keep,]

#drop an motifs not enriched in any clusters
#and drop any clusters without any enriched motifs
nesPlot <- nesPlot[rowSums(nesPlot) != 0,]
nesPlot <- nesPlot[,colSums(nesPlot) != 0]

#normalize scores by row
#(each motif can have a max enrichment score of 1)
nesPlot <- t(apply(nesPlot,1,function(x) x/max(x)))

#bring in motif name/family information
#used to make the rownames more readable/interpretable
motInfo <- read.csv('../../alignment_conservation/motifDB/motifInfo.csv',row.names = 1)

nesPlot.rnames <- mapvalues(rownames(nesPlot),from=motInfo$ID,to=motInfo$name,warn_missing = F)

nesPlot.rnames <- paste0(rownames(nesPlot),' ',nesPlot.rnames)

rownames(nesPlot) <- nesPlot.rnames

#order metagenes to be grouped more logically
mgOrder <- read.csv('mg_order.csv',header = F,row.names = NULL)
mgOrder$V2 <- mgD[mgOrder$V1,'V2']

mgOrder <- mgOrder[mgOrder$V2 %in% colnames(nesPlot),]

nesPlot <- nesPlot[,mgOrder$V2]

#because one of the metagenes got dropped (batch effect mg)
#we have to go back and make sure we have no all 0 rows or columns
nesPlot <- nesPlot[rowSums(nesPlot) != 0,]
nesPlot <- nesPlot[,colSums(nesPlot) != 0]

#plot heatmap
pdf('motifHeatmap.pdf',width = 15,height = 60)
heatmap.2(as.matrix(nesPlot),
          Colv = F,
          scale = 'none',
          dendrogram = 'none',
          col = viridis(30),
          trace='none',
          key = F,
          keysize = 0.1,
          margins = c(10,10),
          colsep = c(3,19,29,38,47),
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="average"))
dev.off()
write.csv(nesPlot,'hydraEnrichmentMatrix.csv')

####reduced redundancy heatmap####

#look at correlation in enrichment patterns across motifs
motEnCor <- apply(enrichDF,1,function(x) apply(enrichDF,1, function(y) cor(x,y)))

motEnCor <- motEnCor[!is.na(motEnCor[,1]),!is.na(motEnCor[,1])]

#generate distance matrix from correlation coefficient
d <- as.dist(1 - motEnCor)

#use hierarchical clustering to group motifs based on correlation
hc1 <- hclust(d, method = "average")

pdf('motDistClustEn.pdf',width=150,height=20)
plot(hc1)
abline(h = 0.4, col = 'red')
dev.off()

#cut the branches at 0.4 and group accordingly
clusts <- cutree(hc1, h = 0.4)

repMotifs <- data.frame(ID = names(clusts), clust = clusts)

#load motif clustering (based on sequence composition alone)
motClust <- read.csv('../../alignment_conservation/motifDB/motif_clusters.csv',row.names = 1)

#extract just the JASPAR part of the motif ID
motClust$jID <- gsub('.*_','',motClust$ID)

#subset to just focus on motifs still being considered for this analysis
motClust <- motClust[motClust$jID %in% repMotifs$ID,]

#make row order identical between the two clustering results
motClust <- motClust[match(repMotifs$ID,motClust$jID),]

#merge clustering results
motClust$enClust <- repMotifs$clust

#combine the two clustering results
#only consider motifs redundant if they
#both clustering analyses grouped them together
motClust$bClust <- paste(motClust$clust,motClust$enClust,sep='_')
write.csv(motClust,file='finalMotClust.csv')

#generating a reduced redundancy version of the heatmap
nesPlotLr <- as.data.frame(nesPlot)

#extract just jaspar ID
nesPlotLr$mot <- gsub(' .*','',rownames(nesPlotLr))

#bring in clustering results
nesPlotLr$motClust <- motClust[match(nesPlotLr$mot,motClust$jID),'bClust']

#sum together enrichment scores across all metagenes for each motif
nesPlotLr$cumScore <- apply(nesPlotLr[,1:53],1,sum)

#select a representative motif from each cluster
#base it on which motif had the strongest
#overall enrichment signal (highest cumulative score)

#for the actuall heatmap, make the values an average
#of all the motifs that were grouped into that cluster

#first get the name of the representative motif
nesPlotLr <- nesPlotLr[order(nesPlotLr$motClust,-nesPlotLr$cumScore),]

nesPlotLr.motNames <- rownames(nesPlotLr[!duplicated(nesPlotLr$motClust),])

#group motifs by cluster
nesPlotLr <- split(nesPlotLr[1:53],nesPlotLr$motClust)

#get average enrichment profile across all motifs with a cluster
nesPlotLr <- lapply(nesPlotLr,function(x) apply(x,2,mean))

#collapse into new low redundancy heatmap matrix
nesPlotLr <- do.call(rbind,nesPlotLr)

rownames(nesPlotLr) <- nesPlotLr.motNames

pdf('motifHeatmapLR.pdf',width = 15,height = 40)
heatmap.2(nesPlotLr,
          Colv = F,
          scale = 'none',
          dendrogram = 'none',
          col = viridis(30),
          trace='none',
          key = F,
          keysize = 0.1,
          margins = c(10,10),
          colsep = c(3,19,29,38,47),
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="average"))
dev.off()
write.csv(nesPlotLr,'hydraEnrichmentMatrixLR.csv')

####Link TFs to Motifs####

#look for putatively linked TFs
#based on pfam domain prediction and
#expression profile

#generate imputed gene expression based on NMF
#makes correlation easier

#reload gene scores
gscores <- t(read.delim('final/whole_unfilt_fine_narrow.gene_spectra_tpm.k_56.dt_0_13.txt',row.names = 1))
rownames(gscores) <- gsub('[.]','_',rownames(gscores))

#multiply the gene and cell scores to create the NMF-derrived approximation for the gene
#expression matrix
gscores <- apply(k56.usage, 1, function(x) apply(gscores,1, function(y) sum(x*y), simplify = T))

gscores <- t(gscores)

gscores <- gscores[rownames(ds@meta.data),]

ds@meta.data$nmfG <- gscores[,which(colnames(gscores) == 'HVAEP1_G023106')]

FeaturePlot(ds,'nmfG',order = T) + NoLegend() + NoAxes()

FeaturePlot(ds,'HVAEP1-G023106',order = T) + NoLegend() + NoAxes()

saveRDS(gscores,'nmfNormExp.rds')
gscores <- readRDS('nmfNormExp.rds')

#get list of genes linked to motifs by pfam domain
pfamMot <- readRDS('../../alignment_conservation/motifDB/motifBindPfam.rds')

pfamMot <- pfamMot[names(pfamMot) %in% gsub(' .*','',rownames(nesPlot))]

#convert to a dataframe
pfamMot <- data.frame(motID = unlist(sapply(1:length(pfamMot), function(x) rep(names(pfamMot)[x],length(pfamMot[[x]])))),
                      geneID = unlist(pfamMot))

pfamMot$geneID <- t2g(pfamMot$geneID)

#subset to genes present in drop-seq object
pfamMot <- pfamMot[gsub('_','-',pfamMot$geneID) %in% rownames(ds),]

#list of TFs
tfList <-read.delim('../../Genome_annotation/functionalAnnotation/tfIDs.txt')[,1,drop=T]

gscores.tf <- gscores[,colnames(gscores) %in% tfList]

#check similarity between motif enrichment and gene expression
corRes <- apply(gscores.tf,2,function(x) apply(motScores, 2, function(y) cor(x,y)))

saveRDS(corRes,'tfMotifCor.rds')

corRes <- readRDS('tfMotifCor.rds')

####Generate Summary Table####
pfamMot.res <- pfamMot

pfamMot.res <- pfamMot.res[pfamMot.res$motID %in% rownames(corRes),]
pfamMot.res <- pfamMot.res[pfamMot.res$geneID %in% colnames(corRes),]

#add in the correlation score calculated above
pfamMot.res$cor <- apply(pfamMot.res,1,function(x){
  corRes[x[1],x[2]]
})

write.csv(pfamMot.res,'pfamMotCor.csv')

#set a correlation score cutoff of 0.5
pfamMot.res <- pfamMot.res[pfamMot.res$cor >= 0.5,]

#bring in motif name
pfamMot.res$motName <- mapvalues(pfamMot.res$motID,from=motInfo$ID,to=motInfo$name,warn_missing = F)

#collapse by gene to simplify
pfamMot.res <- split(pfamMot.res,pfamMot.res$geneID)

#concatenate correlated motifs into a single string, ordered to have most correlated motifs first
pfamMot.res <- lapply(pfamMot.res, function(x) {
  newDf <- x[order(-x$cor),]
  newDf$comboMot <- paste0(newDf$motID,' (',newDf$motName,'):',round(newDf$cor,digits = 3))
  newDf <- data.frame(geneID = newDf[1,'geneID'], corMots = paste(newDf$comboMot,collapse = '; '))
  return(newDf)
})

pfamMot.res <- do.call(rbind,pfamMot.res)

#bring in gene annotations (orthology, domains, blast hits)
gInfo <- read.csv('../../Orthofinder/HVAEP1_annotation.csv')

pfamMot.res$ortho <- mapvalues(pfamMot.res$geneID,from=gsub('HVAEP1_T(\\d+)[.].*','HVAEP1_G\\1',gInfo$H_vulgarisAEP),to = gInfo$EnsemblLongName, warn_missing = F)

pfamMot.res$pfam <- mapvalues(pfamMot.res$geneID,from=gsub('HVAEP1_T(\\d+)[.].*','HVAEP1_G\\1',gInfo$H_vulgarisAEP),to = gInfo$PFAM_NAME, warn_missing = F)

pfamMot.res$gb <- mapvalues(pfamMot.res$geneID,from=gsub('HVAEP1_T(\\d+)[.].*','HVAEP1_G\\1',gInfo$H_vulgarisAEP),to = gInfo$genBankAnnotation, warn_missing = F)

write.csv(pfamMot.res,'regulatorResults.csv')

####Plot correlated motifs/TFs####
unlink('motOlapPlots',recursive=T)
dir.create('motOlapPlots',showWarnings = F)

#ebf motif
ds@meta.data$testMot <- motScores$`MA1637.1`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/ebfM.png',width=9,height=8,dpi=300)
#ebf
FeaturePlot(ds,t2g('HVAEP1_T005780.2'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/ebf.png',width=9,height=8,dpi=300)

#gata motif
ds@meta.data$testMot <- motScores$`MA0766.2`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/gataM.png',width=9,height=8,dpi=300)
#GATA1/2/3
FeaturePlot(ds,t2g('HVAEP1_T022640.1'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/gata.png',width=9,height=8,dpi=300)

#tcf in head
ds@meta.data$testMot <- motScores$`MA1421.1`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/tcfM.png',width=9,height=8,dpi=300)

FeaturePlot(ds,t2g('HVAEP1_T010730.1'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/wnt3.png',width=9,height=8,dpi=300)



#pou
ds@meta.data$testMot <- motScores$`MA0791.1`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/pouM.png',width=9,height=8,dpi=300)

FeaturePlot(ds,t2g('HVAEP1_T023106.1'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/pou.png',width=9,height=8,dpi=300)

#ets extremities
ds@meta.data$testMot <- motScores$`MA0474.2`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/ergM.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1_T001083.1'),order = T,min.cutoff = 1.5,max.cutoff = 2.5,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/erg_cut.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1_T001083.1'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/erg.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1-G001385'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/edl.png',width=9,height=8,dpi=300)

#rfx in gland
ds@meta.data$testMot <- motScores$`MA0509.2`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes() + NoAxes()
ggsave('motOlapPlots/rfx1M.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1-G002581'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/rfx4.png',width=9,height=8,dpi=300)

#zic
ds@meta.data$testMot <- motScores$`MA0697.1`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/zic3M.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1_T004456.1'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/zic4.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1_T004305.1'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/zic1.png',width=9,height=8,dpi=300)

#e2f
ds@meta.data$testMot <- motScores$`MA0471.2`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/e2fM.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1_T001247.1'),order = T,pt.size = 0.8, min.cutoff = 0.7) + NoAxes()
ggsave('motOlapPlots/tfdp_Cut.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1_T001247.1'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/tfdp.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1-G023703'),order = T,pt.size = 0.8, min.cutoff = 0.7) + NoAxes()
ggsave('motOlapPlots/e2f7-8_cut.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1-G023703'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/e2f7-8.png',width=9,height=8,dpi=300)

#fos ecto
ds@meta.data$testMot <- motScores$`MA0478.1`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/fos2lM.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1-G011284'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/fos.png',width=9,height=8,dpi=300)


#otx ecto
ds@meta.data$testMot <- motScores$`MA0712.2`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/otxM.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1-G011623'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/otx.png',width=9,height=8,dpi=300)

#myc
ds@meta.data$testMot <- motScores$`MA0147.3`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/mycM.png',width=9,height=8,dpi=300)

#myc3
FeaturePlot(ds, t2g('HVAEP1-G003730'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/myc3.png',width=9,height=8,dpi=300)

#myc2
FeaturePlot(ds, t2g('HVAEP1_T023705.1'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/myc2.png',width=9,height=8,dpi=300)

#myc1
FeaturePlot(ds, t2g('HVAEP1_T004403.1'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/myc1.png',width=9,height=8,dpi=300)

#foxn
ds@meta.data$testMot <- motScores$`MA0480.1`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/foxoM.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1-G005557'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/foxn.png',width=9,height=8,dpi=300)

#atoh
ds@meta.data$testMot <- motScores$`MA1468.1`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/atoh7M.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1_T021588.1'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/atoh8.png',width=9,height=8,dpi=300)

#nrf nemato
ds@meta.data$testMot <- motScores$`MA0728.1`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/Nr2f6M.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1-G022496'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/Nr2f-like.png',width=9,height=8,dpi=300)


#brachyury 1 head organizer
ds@meta.data$testMot <- motScores$`MA0690.1`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/bra1M.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1-G006413'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/bra1.png',width=9,height=8,dpi=300)

#nk-2/prdl-b nematocyte/endofoot
ds@meta.data$testMot <- motScores$`MA0253.1`

FeaturePlot(ds,'testMot',order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/vndM.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1-G009188'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/nk2.png',width=9,height=8,dpi=300)

FeaturePlot(ds, t2g('HVAEP1-G001197'),order = T,pt.size = 0.8) + NoAxes()
ggsave('motOlapPlots/prdlb.png',width=9,height=8,dpi=300)
