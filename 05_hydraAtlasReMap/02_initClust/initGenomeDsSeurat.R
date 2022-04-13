####Initial R session Setup####
library(Seurat)
library(tidyverse)
library(rstudioapi)
library(glmGamPoi)
library(plotly)
library(RColorBrewer)
library(patchwork)
library(plyr)

#utility to convert transcript ID to gene ID
t2g <- function(x){
  vapply(x, function(y) gsub('HVAEP1_T(\\d+)[.]\\d','HVAEP1-G\\1',y),"")
}

setwd(dirname(getActiveDocumentContext()$path))

####Importing Read Tables####

#get list of count matrix files (include path in name)
readMats <- list.files(pattern = 'dge.txt.gz', recursive = T, full.names = T)

#read in each individual read matrix, generate a separate seurat object for each, and do some initial filtering
inDs <- lapply(readMats, function(x) {
  y <- read.delim(x,stringsAsFactors = F,header = T, row.names = 1)
  
  #get batch name from parent folder of read matrix
  libName <- gsub('[.]\\/(D\\d+[^\\/]+)\\/.*','\\1',x)
  
  #apply batch names to cell barcodes to prevent redundant barcodes
  colnames(y) <- paste0(colnames(y),'-',libName)
  
  #add mt IDs
  rownames(y) <- gsub('^NC','MT-NC',rownames(y))
  
  #initialize seurat object
  tmpDS <- CreateSeuratObject(counts = y, project = libName, min.cells = 3, min.features = 200)
  
  #perform preliminary filtering
  subset(tmpDS, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 & nCount_RNA > 500 & nCount_RNA < 75000)
})

####Batch Correction####
#In the next section we align the individual samples and perform batch correction using reciprocal PCA

#normalize data using SCTransform
inDs <- lapply(inDs, FUN = SCTransform, method = "glmGamPoi")

#use reciprocal PCA to perform batch effect correction and integrate the different libraries.
features <- SelectIntegrationFeatures(inDs)

inDs <- PrepSCTIntegration(object.list = inDs, anchor.features = features)

inDs <- lapply(X = inDs, FUN = RunPCA, features = features, verbose = F)

dsAnchors <- FindIntegrationAnchors(object.list = inDs,
                                    normalization.method = "SCT",
                                    anchor.features = features,
                                    dims = 1:35,
                                    reduction = "rpca",
                                    k.anchor = 5)

ds <- IntegrateData(anchorset = dsAnchors, normalization.method = "SCT", dims = 1:35)

rm(dsAnchors,inDs)

####Clustering and UMAP####

#determining the dimensionality of the data
ds <- RunPCA(ds, npcs = 80)

#VlnPlot(ds, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ElbowPlot(ds,ndims = 80)

ds <- RunUMAP(ds, reduction = "pca", dims = 1:60, min.dist = 0.2, spread=0.2, seed.use = 521, n.neighbors = 45)
ds <- FindNeighbors(ds, reduction = "pca", dims = 1:60)
ds <- FindClusters(ds, resolution = 1.25, graph.name = 'integrated_snn')

allMarks <- FindAllMarkers(ds,only.pos = T,logfc.threshold = 1)

DimPlot(ds,label = T) + NoLegend() + NoAxes()
ggsave('withDubUmap.pdf',width = 6, height=6)
ggsave('withDubUmap.png',width = 6, height=6, units = 'in', dpi = 300)

#plot markers
markList <-read.csv('markerPanel.csv', header = F)

gg <- FeaturePlot(ds,t2g(markList[1,1]),order = T) + NoLegend() + NoAxes() + labs(subtitle = paste0(markList[1,2],'    ',markList[1,3]))

for(i in 2:nrow(markList)){
  subGG <- FeaturePlot(ds,t2g(markList[i,1]),order = T) + NoLegend() + NoAxes() + labs(subtitle = paste0(markList[i,2],'    ',markList[i,3]))
  gg <- gg + subGG
}

png('initReclustMarks.png',width = 25,height = 25,units = 'in',res = 300)
gg + plot_layout(ncol = 6)
dev.off()

####save/load####
saveRDS(ds, file = 'initGenDs.rds')

ds <- readRDS('initGenDs.rds')

###label umap

clusterNames <- data.frame(clusterNumber=levels(ds$seurat_clusters))

clusterNames$names <- c('En_BodyCol/SC','Ec_BodyCol/SC','I_ISC','En_BodyCol/SC','I_FemGc',
                        'I_EarlyNem','Ec_Head','I_MaleGC','I_Neuro','En_Foot','I_StenoNB',
                        'I_SpumousGl','I_GranularGl','En_Head','En_Tent','I_DesmoNB',
                        'I_Neuro','I_EarlyNem','Ec_Peduncle','I_ZymogenGl','I_Ec2N',
                        'Ec_Battery','I_DesmoNC','I_StenoNB','Ec_NB_Dubs','I_DesmoNB',
                        'Ec_BasalDisk','En_Gl_Dubs','En_Hypo','I_En1N','I_DesmoNC',
                        'I_ZymogenGl','I_En2N','I_IsoNB','I_Ec4N','I_maturingNC',
                        'I_Ec1N','Ec_Battery','I_Ec3N','Ec_NB_Dubs','I_SpumousGl',
                        'I_Ec3N','En_BodyCol/SC','I_Ec1N','I_Ec5N','I_StenoNC',
                        'I_En3N','I_Ec1N','I_IsoNC','En_NC_Dubs','I_GlProgen',
                        'En_BodyCol/SC','I_GranularGl','Ec_SomaticGerm','Ec_NC_Dubs')

set.seed(128)
upal <- sample(c("#d69058","#4071f3","#94b92c","#6f46bf","#47b94b",
                 "#c257d1","#7dad46","#956bed","#c9a72e","#4c64cf",
                 "#e28a25","#5d84f4","#df5425","#4c93e9","#ce4230",
                 "#4dbadd","#dc3748","#55b974","#d248ad","#567725",
                 "#aa6fe2","#ae9338","#7747a9","#a9a558","#bc74d8",
                 "#4b9b72","#e03873","#3eb8ad","#d14086","#7b8c50",
                 "#8582e5","#ba5d1f","#4067be","#db8c43","#5959a3",
                 "#8d5d23","#9b87da","#cb5b46","#5790d3","#d9475d",
                 "#5d74a7","#da9d75","#8c509a","#965e3a","#d38edd",
                 "#a35a5a","#a59dda","#b5495a","#c07eb8","#e08584",
                 "#c15dab","#a14767","#d289b0","#9a4979","#e778a3"))

ds$curatedIdent <- as.character(ds$seurat_clusters)
ds$curatedIdent <- factor(mapvalues(ds$curatedIdent, from = clusterNames$clusterNumber, to = clusterNames$names))

ds@active.ident <- ds$curatedIdent

DimPlot(ds,label = T,repel = T,cols = upal) + NoAxes() + NoLegend()
ggsave('labeledDubUmap.pdf',width=8,height=8)
ggsave('labeledDubUmap.png',width=8,height=8,dpi=300)


saveRDS(ds,'labeledDubDs.rds')
ds <- readRDS('labeledDubDs.rds')

####doublet ID

#get ecto markers
ectoMark <- FindMarkers(ds,ident.1 = c(18,53,1,6,37,26),logfc.threshold = 1,only.pos = T)

ectoMark <- rownames(ectoMark[ectoMark$p_val_adj == 0,])

#get endo markers
endoMark <- FindMarkers(ds,ident.1 = c(28,13,0,9,14),logfc.threshold = 1,only.pos = T)

endoMark <- rownames(endoMark[endoMark$p_val_adj == 0,])

#get nemato markers

nemMark <- c(39,24,23,25,33,15,10)

nemMark <- lapply(nemMark, function(x){
  FindMarkers(ds,ident.1 = x,ident.2=1,logfc.threshold = 1,only.pos = T)
})

nemMark <- lapply(nemMark, function(x){
  rownames(x[x$p_val_adj < 1e-200,])
})

batMark <- FindMarkers(ds,ident.1 = 21,ident.2 = c(37,6),only.pos = T,logfc.threshold = 1)

batMark <- rownames(batMark[batMark$p_val_adj < 1e-180,])

#get neuro markers
ectoC <- c(18,53,1,6,37,26)
endoC <- c(28,13,0,9,14)

neuroSubC <- list(20,c(43,47,36),34,c(41,38),32,46,29,44)

neuroMark <- lapply(neuroSubC, function(x){
  FindMarkers(ds,ident.1 = x,ident.2=c(ectoC,endoC),logfc.threshold = 1,only.pos = T)
})

neuroMark <- lapply(neuroMark, function(x){
  rownames(x[x$p_val_adj == 0,])
})


#get sexy markers
sexyMark <- c(4,7)

sexyMark <- lapply(sexyMark, function(x){
  FindMarkers(ds,ident.1 = x,ident.2=2,logfc.threshold = 1,only.pos = T)
})

sexyMark <- lapply(sexyMark, function(x){
  rownames(x[x$p_val_adj == 0,])
})

#get gland markers
glandMark <- list(c(19,12,52),31,c(41,11))

glandMark <- lapply(glandMark, function(x){
  FindMarkers(ds,ident.1 = x, ident.2 = endoC, logfc.threshold = 1,only.pos = T)
})

glandMark <- lapply(glandMark, function(x){
  rownames(x[x$p_val_adj == 0,])
})

markList <- list(ectoMark,endoMark,sexyMark[[1]],sexyMark[[2]],batMark)

markList <- append(markList,values = c(neuroMark,glandMark,nemMark))

markNames <- c('ecto','endo','fem','male','nemBat',
               paste0('neur',LETTERS[1:length(neuroMark)]),
               paste0('gland',LETTERS[1:3]),
               paste0('nem',LETTERS[1:7]))

ds <- AddModuleScore(ds,features=markList,
                     name = markNames,
                     assay = 'SCT')


png('initDubSubScores.png',width = 20,height = 25,units = 'in',res = 300)
FeaturePlot(ds,paste0(markNames,1:length(markNames)), order=T)
dev.off()


dubs <- ds@meta.data[,paste0(markNames,1:length(markNames))]

dubs <- as.data.frame(dubs > 0.2)

dubs$ectoDub <- apply(dubs[,-2],1,sum)

dubs$ectoDub[!dubs$ecto1] <- 0

dubs$ectoDub <- dubs$ectoDub > 1

dubs$endoDub <- apply(dubs[,-1],1,sum)

dubs$endoDub[!dubs$endo2] <- 0

dubs$endoDub <- dubs$endoDub > 1

dubs$endoEcto <- (dubs$ecto1 + dubs$endo2) > 1

dubs$dubTot <- dubs$ectoDub | dubs$endoDub | dubs$endoEcto

ds@meta.data$dubTest <- as.numeric(dubs$dubTot)

png(filename = 'initDsClust.png',height = 4,width = 4,units = 'in',res = 300)
DimPlot(ds,label = T) + NoLegend() + NoAxes()
dev.off()

png(filename = 'initDsDubs.pdf',height = 6,width = 6,units = 'in',res = 300)
FeaturePlot(ds,c('dubTest')) + NoLegend() + NoAxes()
dev.off()

#export non-Doublets

write.table(rownames(ds@meta.data[(!ds@meta.data$dubTest),]), file = 'nondub.tsv',quote = F, row.names = F, col.names = F)


####NMF prep####
rawC <- t(as.matrix(ds@assays$RNA@counts))
write.table(rawC,file="unfilt.whole.raw.counts.tsv",sep = '\t', quote = F)
rm(rawC)

normC <- t(as.matrix(ds@assays$SCT@data))
write.table(normC,file="unfilt.whole.norm.counts.tsv",sep = '\t', quote = F)
rm(normC)

#list of variable genes to look at
write.table(ds@assays$integrated@var.features, file = 'unfilt.whole.genes.tsv', row.names = F, col.names = F, quote = F)

