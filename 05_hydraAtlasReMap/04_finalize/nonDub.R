library(Seurat)
library(tidyverse)
library(rstudioapi)
library(glmGamPoi)
library(plotly)
library(RColorBrewer)
library(patchwork)
library(ggplot2)
library(plyr)

#utility to convert transcript ID to gene ID
t2g <- function(x){
  vapply(x, function(y) gsub('HVAEP1_T(\\d+)[.]\\d','HVAEP1-G\\1',y),"")
}

setwd(dirname(getActiveDocumentContext()$path))

nonDub <- read.delim('../nmf/nonDubIDs.txt')[,1]

#get list of count matrix files (include path in name)
readMats <- list.files(pattern = 'dge.txt.gz', recursive = T, full.names = T)

#read in each individual read matrix, generate a separate seurat object for each, and do some initial filtering
inDs <- lapply(readMats, function(x) {
  y <- read.delim(x,stringsAsFactors = F,header = T, row.names = 1)
  
  #get batch name from parent folder of read matrix
  libName <- gsub('[.]\\/(D\\d+[^\\/]+)\\/.*','\\1',x)
  
  #apply batch names to cell barcodes to prevent redundant barcodes
  colnames(y) <- paste0(colnames(y),'-',libName)
  
  y <- y[,colnames(y) %in% nonDub]
  
  #initialize seurat object
  tmpDS <- CreateSeuratObject(counts = y, project = libName, min.cells = 3, min.features = 200)
  
  #perform preliminary filtering
  subset(tmpDS, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 & nCount_RNA > 500 & nCount_RNA < 75000)
})

#normalize data using SCTransform
inDs <- lapply(inDs, FUN = SCTransform, method = "glmGamPoi")

features <- SelectIntegrationFeatures(inDs,nfeatures = 2000)

inDs <- PrepSCTIntegration(object.list = inDs, anchor.features = features)

inDs <- lapply(X = inDs, FUN = RunPCA, features = features, verbose = F)

dsAnchors <- FindIntegrationAnchors(object.list = inDs,
                                    normalization.method = "SCT",
                                    anchor.features = features,
                                    dims = 1:45,
                                    reduction = "rpca",
                                    k.anchor = 5)

rm(inDs)
gc()

ds <- IntegrateData(anchorset = dsAnchors, normalization.method = "SCT", dims = 1:45)

rm(dsAnchors)

#determining the dimensionality of the data
ds <- RunPCA(ds, npcs = 80)

ElbowPlot(ds,ndims = 80)

#generate UMAP projection and perform louvain clustering

ndUse=55
ds <- RunUMAP(ds, reduction = "pca", dims = 1:ndUse, min.dist = 0.3, spread=0.38, seed.use = 600,n.neighbors = 50)
ds <- FindNeighbors(ds, reduction = "pca", dims = 1:ndUse)
ds <- FindClusters(ds, resolution = 0.7, graph.name = 'integrated_snn')

DimPlot(ds,label = T) + NoLegend() + NoAxes()



####look at markers####
DefaultAssay(ds) <- 'SCT'

markList <-read.csv('markerPanel.csv', header = F)

gg <- FeaturePlot(ds,t2g(markList[1,1]),order = T) + 
  NoAxes() + labs(subtitle = paste0(markList[1,2],'    ',markList[1,3])) + 
  theme(legend.key.height = unit(0.2, 'cm')) + 
  theme(legend.text = element_text(size=8,face = "bold")) + 
  theme(legend.position="bottom")

for(i in 2:nrow(markList)){
  subGG <- FeaturePlot(ds,t2g(markList[i,1]),order = T) + 
    NoAxes() + labs(subtitle = paste0(markList[i,2],'    ',markList[i,3])) + 
    theme(legend.key.height = unit(0.2, 'cm')) + 
    theme(legend.text = element_text(size=8,face = "bold")) + 
    theme(legend.position="bottom")
  gg <- gg + subGG
}

png('nonDubMarks.png',width = 20,height = 24,units = 'in',res = 300)
gg + plot_layout(ncol = 6)
dev.off()

#get markers for mystery endo clusters

#cluster 41
mark41 <- FindMarkers(ds,ident.1 = 41,ident.2 = c(5,7,0),only.pos = T,logfc.threshold = 1,assay = 'integrated')

FeaturePlot(ds,rownames(mark41)[1:9],order = T)

#cluster 37
mark37 <- FindMarkers(ds,ident.1 = 37, ident.2 = c(5,7,0,9),only.pos = T,logfc.threshold = 1,assay = 'integrated')

annots <- read.delim('../../Genome_annotation/functionalAnnotation/upBlast/uniprotBlast.txt', header = F)[,1:2]

annots$V1 <- t2g(annots$V1)

mark37$id <- rownames(mark37)

mark37 <- merge(mark37,annots, by.x = 'id', by.y = 'V1',all.x = T)

mark37 <- mark37[order(mark37$p_val),]

FeaturePlot(ds,mark37$id[1:9],order = T)

#cluster 37 is stress

#drop the clusters 41 and 37

finalNonDub <- rownames(ds@meta.data[!(ds@meta.data$seurat_clusters %in% c(37,41)),])

ds <- subset(ds,cells = finalNonDub)

DefaultAssay(ds) <- 'integrated'

#determining the dimensionality of the data
ds <- RunPCA(ds, npcs = 80)

#generate UMAP projection and perform louvain clustering

# 0,21 0,2 4321

ndUse=55
ds <- RunUMAP(ds, reduction = "pca", dims = 1:ndUse, min.dist = 0.2, spread=0.2, seed.use = 963140, n.neighbors = 45)
ds <- FindNeighbors(ds, reduction = "pca", dims = 1:ndUse)
ds <- FindClusters(ds, resolution = 0.7, graph.name = 'integrated_snn')

DefaultAssay(ds) <- 'SCT'

DimPlot(ds,label = T) + NoLegend() + NoAxes()
ggsave('unlabeledNonDubUmap.pdf',width=8,height=8)
ggsave('unlabeledNonDubUmap.png',width=8,height=8,dpi=300)

write.table(colnames(ds),'finalNonDub.txt', row.names = F, col.names = F, quote = F)

saveRDS(ds,file='nonDubSeurat.rds')
ds <- readRDS('nonDubSeurat.rds')

#markers again
gg <- FeaturePlot(ds,t2g(markList[1,1]),order = T) + NoLegend() + NoAxes() + labs(subtitle = paste0(markList[1,2],'    ',markList[1,3]))

for(i in 2:nrow(markList)){
  subGG <- FeaturePlot(ds,t2g(markList[i,1]),order = T) + NoLegend() + NoAxes() + labs(subtitle = paste0(markList[i,2],'    ',markList[i,3]))
  gg <- gg + subGG
}

png('nonDubMarks.png',width = 20,height = 25,units = 'in',res = 300)
gg + plot_layout(ncol = 6)
dev.off()

####Label Final UMAP####
clusterNames <- data.frame(clusterNumber=levels(ds$seurat_clusters))

clusterNames$names <- c('En_BodyCol/SC','Ec_BodyCol/SC','I_ISC','I_Neuro','I_EarlyNem',
                        'I_FemGC','En_Head','En_BodyCol/SC','I_MaleGC','En_Foot',
                        'I_StenoNB','Ec_Head','I_DesmoNB','I_GranGl','En_Tentacle',
                        'I_Ec2N','I_SpumMucGl','I_DesmoNC','I_DesmoNB','I_ZymoGl',
                        'I_SpumMucGl','Ec_BasalDisk','Ec_Peduncle','I_Ec1N','I_En1N',
                        'I_DesmoNC','I_StenoNB','I_StenoNC','I_En2N','I_IsoNB',
                        'I_ZymoGl','Ec_Tentacle','I_Ec4N','I_Ec1/5N','I_Ec3N',
                        'I_IsoNC','I_Ec3N','I_En3N','I_GranGl','I_GlProgen')

ds$curatedIdent <- as.character(ds$seurat_clusters)
ds$curatedIdent <- factor(mapvalues(ds$curatedIdent, from = clusterNames$clusterNumber, to = clusterNames$names))

ds@active.ident <- ds$curatedIdent

saveRDS(ds,'nonDubLabeledSeurat.rds')
ds <- readRDS('nonDubLabeledSeurat.rds')

set.seed(12345)
colsUse <- sample(c("#eba92c","#984eec","#79dc3a","#d248de",
                    "#58ce5f","#d84fc4","#d1d52e","#7766e8",
                    "#9ccd4f","#e4349d","#5bda8f","#ad69d7",
                    "#cbc553","#5581f2","#e64f22","#3be9ca",
                    "#e3407a","#53ad70","#be5aa7","#759027",
                    "#7477da","#89b65d","#9a7ad0","#bf9439",
                    "#3e76cc","#cb7329","#4b96eb","#e04946",
                    "#de8cd8","#d0764d","#d75e95","#da5666"))

DimPlot(ds,label = T,repel = T,cols = colsUse) + NoAxes() + NoLegend()
ggsave('labeledUMAP.pdf',width = 8, height = 8)
ggsave('labeledUMAP.png',width = 8, height = 8,dpi=300)


DimPlot(ds,label = F,repel = T,cols = colsUse) + NoAxes() + NoLegend()

ds$lineage <- gsub('_.*','',ds$curatedIdent)

colsUse <- c('#4E8ECC','#9FCC58','#F27755')
DimPlot(ds,label=F,cols=colsUse, group.by = 'lineage') + NoAxes() + NoLegend()
ggsave('unlabeledUMAPwLineage.pdf',width = 8, height = 8)


FeaturePlot(ds,t2g('HVAEP1_T003934.1'),order = T)


####Misc Marker Search####
#get putative isorhiza markers

annots <- read.csv('../../Orthofinder/HVAEP1_annotation.csv')

ipr <- read.delim('../../Genome_annotation/functionalAnnotation/HVAEP1.prot.longestIso.fa.tsv', header = F)

DimPlot(ds,label=T)

isoMark1 <- FindMarkers(ds,ident.1 = 'I_IsoNC',ident.2 = c('I_StenoNC','I_DesmoNC','I_StenoNB','I_DesmoNB'),only.pos = T,logfc.threshold = 1)
isoMark1$id <- rownames(isoMark1)

isoMark1 <- isoMark1[order(isoMark1$p_val_adj),]

gg <- FeaturePlot(ds,isoMark1$id[1],order = T) + NoAxes()

for(i in 2:30){
  subGG <- FeaturePlot(ds,isoMark1$id[i],order = T) + NoAxes()
  gg <- gg + subGG
}

png('iso1Marks.png',width = 25,height = 22,units = 'in',res = 300)
gg + plot_layout(ncol = 6)
dev.off()

FeaturePlot(ds,rownames(isoMark1[1:20,]))

isoMark2 <- FindMarkers(ds,ident.1 = 'I_IsoNC',only.pos = T,logfc.threshold = 1)
isoMark2$id <- rownames(isoMark2)

isoMark2 <- isoMark2[order(isoMark2$p_val_adj),]

gg <- FeaturePlot(ds,isoMark2$id[1],order = T) + NoAxes()

for(i in 2:30){
  subGG <- FeaturePlot(ds,isoMark2$id[i],order = T) + NoAxes()
  gg <- gg + subGG
}

png('iso2Marks.png',width = 25,height = 22,units = 'in',res = 300)
gg + plot_layout(ncol = 6)
dev.off()

#try to pick good in situ candidates

#mature isorhiza
FeaturePlot(ds,'HVAEP1-G027705')

FeaturePlot(ds,'HVAEP1-G010119')

FeaturePlot(ds,'HVAEP1-G017166')

FeaturePlot(ds,'HVAEP1-G008733',order=T,pt.size = 0.8) + NoAxes()
ggsave('isoIshGene.png',width=9,height = 8,dpi=300)

#developing isorhiza
FeaturePlot(ds,'HVAEP1-G015357')


FeaturePlot(ds,'HVAEP1-G005822')

#scleraxis
FeaturePlot(ds, 'HVAEP1-G017021',order=T,pt.size = 0.8) + NoAxes()
ggsave('scleraxis.png',width=9,height=8,dpi=300)

#iso ion channel
FeaturePlot(ds, 'HVAEP1-G002529',order=T,pt.size = 0.8) + NoAxes()
ggsave('isoIon.png',width=9,height=8,dpi=300)

#generate comprehensive marker list
dsMarks <- FindAllMarkers(ds,logfc.threshold = 1, only.pos = T)

write.csv(dsMarks,'nonDubMarkers.csv',row.names = F)
