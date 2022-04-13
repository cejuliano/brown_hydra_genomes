library(Seurat)
library(tidyverse)
library(rstudioapi)
library(glmGamPoi)
library(plotly)
library(RColorBrewer)
library(patchwork)
library(plyr)
library(SeuratDisk)
library(cowplot)
library(gplots)
library(viridis)

#utility to convert transcript ID to gene ID
t2g <- function(x){
  vapply(x, function(y) gsub('HVAEP1_T(\\d+)[.]\\d','HVAEP1-G\\1',y),"")
}

setwd(dirname(getActiveDocumentContext()$path))

#import integrated hydra and clytia seurat object
aepCl.int <- readRDS('aepClInt.rds')

#perform a high resolution louvain clustering analysis
#this generates ad-hoc 'pseudo-cells' that group 
#together small groups of cells with similar expression
#(total of 132 clusters at this resolution)
aepCl.int <- FindClusters(aepCl.int, resolution = 10, graph.name = 'integrated_snn')

#pseudo-cell UMAP
DimPlot(aepCl.int) + NoAxes() + NoLegend()

#extract normalized expression data
expDat <- t(as.matrix(aepCl.int@assays$SCT@data))

#get the average expression of each gene, aggregated by both pseudo-cell ID
#and species
expDat.cl <- aggregate(expDat,list(aepCl.int$seurat_clusters,aepCl.int$species),mean)

#only keep the clusters that had cells from both species
expDat.cl <- expDat.cl[expDat.cl$Group.1 %in% expDat.cl[expDat.cl$Group.2 == 'Cl','Group.1'],]
expDat.cl <- expDat.cl[expDat.cl$Group.1 %in% expDat.cl[expDat.cl$Group.2 == 'AEP','Group.1'],]

write.csv(expDat.cl,file='pcExpression.csv',row.names = F)

#split data by species (make an AEP df and a Clytia DF)
expDat.cl.list <- split(expDat.cl,expDat.cl$Group.2)

#droup group.1 and group.2 columns
expDat.cl.list <- lapply(expDat.cl.list,function(x){
  x <- x[,c(-1,-2)]
  return(x)
})

#calculate the correlation in pseudocell expression patterns across
#the two species
corRes <- lapply(1:ncol(expDat.cl.list[[1]]),function(x){
  cor(expDat.cl.list[[1]][,x],expDat.cl.list[[2]][,x],method = 'pearson')
})

#generate DF of correlation scores for all genes
corRes <- do.call(c,corRes)

corRes <- data.frame(id = colnames(expDat.cl.list[[1]]),cor=corRes)

#set NAs to 0
corRes[is.na(corRes$cor),'cor'] <- 0

#order by correlation score
corRes <- corRes[order(-corRes$cor),]

write.csv(corRes,file='geneCoreRes.csv',row.names = F)

#after identifying genes with similar expression patterns
#convey the data using a heatmap based on a lower resolution clustering
#otherwise it's just an insane number of columns

#import more descriptive names for lower resolution clusters
corPlotLabels <- read.csv('corPlotClustsLables.csv',header=F)

#make umap with descriptive names for clusters used for heatmap plots
aepCl.int$renameClust <- mapvalues(aepCl.int$integrated_snn_res.0.4,from = corPlotLabels$V1, corPlotLabels$V2, warn_missing = F)

DimPlot(aepCl.int,group.by = 'renameClust',label = T, repel = T) + NoLegend() + NoAxes()
ggsave('corPlotClusts.pdf',width=8,height=8)

#generate averaged expression data for lower resolution clustering
heatObj <- aggregate(expDat,list(aepCl.int$integrated_snn_res.0.4,aepCl.int$species),mean)

heatObj <- heatObj[heatObj$Group.1 %in% heatObj[heatObj$Group.2 == 'Cl','Group.1'],]
heatObj <- heatObj[heatObj$Group.1 %in% heatObj[heatObj$Group.2 == 'AEP','Group.1'],]

#bring in the descriptive cluster names
heatObj$Group.1 <- mapvalues(heatObj$Group.1,
                             from=corPlotLabels$V1,
                             to=corPlotLabels[,2],
                             warn_missing = F)

#split by species
heatObj.list <- split(heatObj,heatObj$Group.2)

#drop group.1 and group.2 columns
heatObj.list <- lapply(heatObj.list,function(x){
  rownames(x) <- x[,1]
  x <- x[,c(-1,-2)]
  return(x)
})

#get data for AEP heatmap
heatObj.1 <- heatObj.list[[1]]

#only plot genes with high correlation scores
heatObj.1 <- heatObj.1[,colnames(heatObj.1) %in% corRes[corRes$cor > 0.65,'id']]

#import table that specifies a reasonably logical way of ordering
#the heatmap columns (cell clusters)
corPloOrder <- read.csv('corPlotClustsOrder.csv',header=F)

heatObj.1 <- heatObj.1[corPloOrder$V2,]

pdf('heat1.pdf',width = 8, height = 40)
#plot heatmap, save the output, which specifies how rows are ordered
#which is needed to keep the clytia heatmap consistent
dgm <- heatmap.2(t(heatObj.1),
                 scale = 'row',
                 dendrogram = 'none',
                 Colv = F,
                 col = viridis(30),
                 trace='none',
                 key = F,
                 keysize = 0.1,
                 margins = c(10,6),
                 distfun = function(x) as.dist(1-cor(t(x))),
                 hclustfun = function(x) hclust(x, method="average"))
dev.off()

write.csv(t(heatObj.1),file='crossSpecMatHydra.csv',row.names = F)

#extract clytia heatmap data
heatObj.2 <- heatObj.list[[2]]

#keep only genes with good correlation scores
heatObj.2 <- as.data.frame(t(heatObj.2[,colnames(heatObj.2) %in% corRes[corRes$cor > 0.65,'id']]))

#reorder data to have the same order as the AEP matrix
heatObj.2 <- heatObj.2[rev(dgm$rowInd),corPloOrder$V2]

#change gene names back to original clytia ortho name
orthos <- read.delim('../../Orthofinder/Results_Sep15_1/Phylogenetic_Hierarchical_Orthogroups/N14.tsv')

orthos <- orthos[,c('C_hemisphaerica','H_vulgarisAEP')]

orthos <- orthos[orthos$C_hemisphaerica != '' & orthos$H_vulgarisAEP != '',]

orthos <- orthos[!(grepl(',',orthos$C_hemisphaerica) | grepl(',',orthos$H_vulgarisAEP)),]

orthos$H_vulgarisAEP <- gsub('HVAEP1_T(\\d+)[.].*','HVAEP1-G\\1',orthos$H_vulgarisAEP)

rownames(heatObj.2) <- mapvalues(rownames(heatObj.2), from = orthos$H_vulgarisAEP, to = orthos$C_hemisphaerica,warn_missing = F)

#plot clytia heatmap
pdf('heat2.pdf',width = 8, height = 40)
heatmap.2(as.matrix(heatObj.2),
          Rowv = F,
          Colv = F,
          scale = 'row',
          dendrogram = 'none',
          col = magma(30),
          trace='none',
          key = F,
          keysize = 0.1,
          margins = c(10,6),
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="average"))
dev.off()

write.csv(t(heatObj.2),file='crossSpecMatClytia.csv',row.names = F)

#generate table of correlation scores for export (only high scoring genes)

#bring in functional annotation information to make browsing results easier to interpret
annots <- read.csv('../../Orthofinder/HVAEP1_annotation.csv')

annots$H_vulgarisAEP <- t2g(annots$H_vulgarisAEP)

corRes.annot <- merge(corRes,annots,by.x = 'id', by.y = 'H_vulgarisAEP', all.x = T)

corRes.annot <- corRes.annot[order(-corRes.annot$cor),]

corRes.annot <- corRes.annot[corRes.annot$cor > 0.65,]

corRes.annot$clOrtho <- mapvalues(corRes.annot$id, from = orthos$H_vulgarisAEP, to = orthos$C_hemisphaerica, warn_missing = F)
write.csv(corRes.annot,'crossSpecExpCor.csv',row.names = F)

#subset correlation heatmap to just look at trancription factors
tfList <- read.delim('../../Genome_annotation/functionalAnnotation/tfIDs.txt',header=F)[,1,drop=T]

tfList <- gsub('_','-',tfList)

heatObj.1 <- heatObj.list[[1]]

heatObj.1 <- heatObj.1[,colnames(heatObj.1) %in% corRes.annot[corRes.annot$cor > 0.65,'id']]

heatObj.1 <- heatObj.1[,colnames(heatObj.1) %in% tfList]

heatObj.1 <- heatObj.1[corPloOrder$V2,]

pdf('tfHeat1.pdf',width = 8, height = 10)
dgm <- heatmap.2(t(heatObj.1),
                 scale = 'row',
                 Colv = F,
                 dendrogram = 'none',
                 col = viridis(30),
                 trace='none',
                 key = F,
                 keysize = 0.1,
                 margins = c(5,10),
                 distfun = function(x) as.dist(1-cor(t(x))),
                 hclustfun = function(x) hclust(x, method="average"))
dev.off()

write.csv(t(heatObj.1),file='crossSpecTfMatHydra.csv',row.names = F)

heatObj.2 <- heatObj.list[[2]]

heatObj.2 <- as.data.frame(t(heatObj.2[,colnames(heatObj.2) %in% corRes.annot[corRes.annot$cor > 0.65,'id']]))

heatObj.2 <- heatObj.2[rownames(heatObj.2) %in% tfList,]

heatObj.2 <- heatObj.2[rev(dgm$rowInd),corPloOrder$V2]

pdf('tfHeat2.pdf',width = 8, height = 10)
heatmap.2(as.matrix(heatObj.2),
          Rowv = F,
          Colv = F,
          scale = 'row',
          dendrogram = 'none',
          col = magma(30),
          trace='none',
          key = F,
          keysize = 0.1,
          margins = c(5,10))
dev.off()

write.csv(t(heatObj.2),file='crossSpecTfMatClytia.csv',row.names = F)
