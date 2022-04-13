library(Seurat)
library(SeuratDisk)
library(plyr)
library(tidyverse)
library(rstudioapi)
library(glmGamPoi)
library(plotly)
library(RColorBrewer)
library(reticulate)
library(Matrix)

setwd(dirname(getActiveDocumentContext()$path))

#convert annData object from Clytia scRNA-seq publication to make it easier to import into Seurat
#only need to run once
#Convert("bus_fs_combo_raw.h5ad", dest = "h5seurat", overwrite = TRUE)

#load clustering analysis object from Clytia scRNA-seq paper
clPP <- LoadH5Seurat('bus_fs_combo_raw.h5seurat')

#import newly mapped Clytia gene expression matrix
inMat <- Read10X_h5('raw_feature_bc_matrix.h5')

#drop any cell IDs that aren't in the published clustering data object
inMat <- inMat[,colnames(inMat) %in% rownames(clPP@meta.data)]

#create seurat object
#drop genes expressed in fewer than 3 cells and cells with fewer than 200 genes
cl <- CreateSeuratObject(counts = inMat, project = 'clDS', min.cells = 3, min.features = 200)

VlnPlot(cl, features = c("nFeature_RNA", "nCount_RNA"))

cl <- subset(cl, nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 100000)

cl <- SCTransform(cl, method = "glmGamPoi")

cl <- RunPCA(cl, npcs = 60, verbose = FALSE)

ElbowPlot(cl,ndims = 60)

cl <- RunUMAP(cl, dims = 1:45, min.dist = 0.38,verbose = FALSE, seed.use = 12345)

cl <- FindNeighbors(cl, dims = 1:45, verbose = FALSE)
cl <- FindClusters(cl, resolution = 0.5, verbose = FALSE)

DimPlot(cl, label=T) + NoLegend() + NoAxes()
ggsave('clRemapUMAP.pdf',width=8,height=8)
ggsave('clRemapUMAP.png',width=8,height=8,dpi=300)

saveRDS(cl,file = 'initSeurat.rds')
cl <- readRDS('initSeurat.rds')


#import annotations from preprint
use_condaenv("rScanpy", required = T)

sc <- import("scanpy")

adata <- sc$read_h5ad('fedStarved_withUMAPPaga.h5ad')

exprs <- t(adata$X)
colnames(exprs) <- adata$obs_names$to_list()
rownames(exprs) <- adata$var_names$to_list()
# Create the clPP object
clPP <- CreateSeuratObject(exprs)
# Set the expression assay
clPP <- SetAssayData(clPP, "data", exprs)
# Add observation metadata
clPP <- AddMetaData(clPP, adata$obs)
# Add embedding
embedding <- adata$obsm["X_umap"]
rownames(embedding) <- adata$obs_names$to_list()
colnames(embedding) <- c("umap_1", "umap_2")
clPP[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")

#save reformatted object for use later
saveRDS(clPP,'annotatedCl.rds')

DimPlot(clPP,group.by = 'annosSub', label = T, repel = T) + NoLegend() + NoAxes()
ggsave('clOrigUMAP.pdf',width=8,height=8)
ggsave('clOrigUMAP.png',width=8,height=8,dpi=300)

saveRDS(clPP,file = 'processedCl.rds')

cl@meta.data$origAnnot <- mapvalues(rownames(cl@meta.data), from=rownames(clPP@meta.data), to = as.character(clPP@meta.data$annosSub), warn_missing = F)

set.seed(12345)
umapPal <- sample(c("#b76749","#6d6be8","#94b92c","#a358d6","#52b648","#cf43a8",
                    "#4ac38b","#ab52b7","#59a85a","#5c58c0","#c9a72e","#5581f2",
                    "#e28624","#5693dd","#d6482c","#45c1b8","#d44248","#47afd4",
                    "#c56633","#a187e3","#8b9e41","#df82da","#4c976f","#d14a87",
                    "#238e7e","#d24461","#767d3e","#525ea7","#bc8d39","#9565ab",
                    "#bc8f57","#8b8dc7","#c26b6c","#b96daa","#b96580","#b66c96"))

DimPlot(cl,group.by = 'origAnnot', label = T,repel = T,cols = umapPal,pt.size=0.8) + NoLegend() + NoAxes()
ggsave('clRemapUmapLabeled.pdf',width = 15, height = 15)
ggsave('clRemapUMAPorigAnnot.png',width=8,height=8,dpi=300)