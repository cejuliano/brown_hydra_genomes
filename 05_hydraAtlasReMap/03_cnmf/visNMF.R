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

#import labeled seurat object
ds <- readRDS('../dropSeqMapping/nonDubLabeledSeurat.rds')

#import NMF cell scores 
k56.usage <- read.delim('final/whole_unfilt_fine_narrow.usages.k_56.dt_0_13.consensus.txt',row.names = 1)

#normalize cell score values so that the sum of scores for each cell sums to 1
k56.usage <- t(apply(k56.usage,1,function(x) x/sum(x)))

#add the metagene scores to the Seurat object metadata table
md <- ds@meta.data[,1:8]

md$id <- rownames(md)

k56.usage <- as.data.frame(k56.usage)

colnames(k56.usage) <- gsub('X','mg',colnames(k56.usage))

#incorporate descriptions of each metagene into the metagene name
mgD <- read.csv('k56_mg_annot.csv',header = F)

colnames(k56.usage) <- paste(colnames(k56.usage), mgD$V2,sep=' ')

k56.usage$id <- rownames(k56.usage) 

md <- merge(md, k56.usage, by = 'id',all.x = T)

rownames(md) <- md$id

md$id <- NULL

ds@meta.data <- md

#get the names of the metagene columns in the metadata DF
plotMods <- grep('mg',colnames(ds@meta.data),value = T)

#plot all metagene cell scores on the atlas UMAP
gg <- FeaturePlot(ds,plotMods[1],order=T) + NoLegend() + NoAxes()

for(i in 2:length(plotMods)){
  subGG <- FeaturePlot(ds,t2g(plotMods[i]),order=T) + NoLegend() + NoAxes()
  gg <- gg + subGG
}

#how many rows/columns should the plot have?
dims <- sqrt(56)
dims <- c(ceiling(dims),floor(dims))

pdf('k56usage.pdf',width = 4*dims[1],height = 4*dims[2])
gg + plot_layout(ncol = dims[1])
dev.off()

png('k56usage.png',width = 4*dims[1],height = 4*dims[2],units='in',res=300)
gg + plot_layout(ncol = dims[1])
dev.off()


#sexy metagene batch composition

FeaturePlot(ds,'X2')

summary(k56.usage$X2)

# ds$sexyTest <- ds@meta.data$X2 > 0.15
# 
# FeaturePlot(ds,'sexyTest') + NoLegend() + NoAxes()

sexyMG <- data.frame(cell.id = colnames(ds),
                     mg.member = ds@meta.data$X2 > 0.15,
                     batch = ds@meta.data$orig.ident,
                     clust = ds@meta.data$curatedIdent)

sexyMG <- sexyMG[sexyMG$clust == 'Ec_BodyCol/SC',]

sexyMG <- split(sexyMG,sexyMG$mg.member)

sexyMG <- lapply(sexyMG,function(x) aggregate(x,list(x$batch),length))

sexyMG <- lapply(1:2,function(x) {
  newDF.name <- names(sexyMG)[x]
  newDF <- sexyMG[[x]]
  newDF$cell.id <- newDF$cell.id/sum(newDF$cell.id)
  newDF$mg.member <- as.logical(newDF.name)
  return(newDF)
})

sexyMG <- do.call(rbind,sexyMG)

sexyMG <- sexyMG[order(-sexyMG$cell.id),]

sexyMG$Group.1 <- factor(sexyMG$Group.1,unique(c(sexyMG[sexyMG$mg.member == TRUE,'Group.1'],sexyMG$Group.1)))

ggplot(sexyMG,aes(x=mg.member,y=cell.id)) + geom_col(aes(fill = Group.1)) + theme_bw()

ggsave('sexyMgBatch.pdf',width=6,height=8)
