library(Seurat)
library(rstudioapi)
library(patchwork)

setwd(dirname(getActiveDocumentContext()$path))

#import UMAP from original clytia publication
cl <- readRDS('annotatedCl.rds')

#k37 course NMF module score visualization

#import nmf cell scores
k37.usage <- read.delim('nmf/cl_fine.usages.k_37.dt_0_13.consensus.txt',row.names = 1)

#normalize cell scores so that all scores in a cell sum to 1
k37.usage <- t(apply(k37.usage,1,function(x) x/sum(x)))

#append metagene cell scores to seurat metadata
md <- cl@meta.data

md$id <- rownames(md)

k37.usage <- as.data.frame(k37.usage)

colnames(k37.usage) <- gsub('X','mg',colnames(k37.usage))

#import metagene descriptions
#(based on initial review of cell score plots)
mgD <- read.csv('nmf/clMetaAnnot.csv',header = F)

colnames(k37.usage) <- paste(colnames(k37.usage), mgD$V2,sep=' ')

k37.usage$id <- rownames(k37.usage) 

md <- merge(md, k37.usage, by = 'id',all.x = T)

md[is.na(md)] <- 0

rownames(md) <- md$id

md$id <- NULL

cl@meta.data <- md

#plot metagene cell scores on the clytia UMAP
plotMods <- grep('mg',colnames(cl@meta.data),value = T)

gg <- FeaturePlot(cl,plotMods[1],order=T) + NoLegend() + NoAxes()

for(i in 2:length(plotMods)){
  subGG <- FeaturePlot(cl,plotMods[i],order=T) + NoLegend() + NoAxes()
  gg <- gg + subGG
}

dims <- sqrt(37)
dims <- c(ceiling(dims),floor(dims))

png('k37usage.png',4*dims[1],height = 5*dims[2],units = 'in',res = 300)
gg + plot_layout(dims[1])
dev.off()

set.seed(12345)
umapPal <- sample(c("#b76749","#6d6be8","#94b92c","#a358d6","#52b648","#cf43a8",
             "#4ac38b","#ab52b7","#59a85a","#5c58c0","#c9a72e","#5581f2",
             "#e28624","#5693dd","#d6482c","#45c1b8","#d44248","#47afd4",
             "#c56633","#a187e3","#8b9e41","#df82da","#4c976f","#d14a87",
             "#238e7e","#d24461","#767d3e","#525ea7","#bc8d39","#9565ab",
             "#bc8f57","#8b8dc7","#c26b6c","#b96daa","#b96580","#b66c96"))

DimPlot(cl,label=T,repel=T,cols=umapPal,group.by = 'annosSub',pt.size=0.8) + NoAxes() + NoLegend()
ggsave('clUmapLabeled.pdf',width = 15, height = 15)


DimPlot(cl,cols=umapPal,group.by = 'annosSub',pt.size=0.8) + NoAxes()
ggsave('clUmapLabeledLegend.pdf',width = 25, height = 15)




