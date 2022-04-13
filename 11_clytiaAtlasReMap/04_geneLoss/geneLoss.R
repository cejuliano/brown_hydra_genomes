library(rstudioapi)
library(Seurat)
library(plyr)
library(patchwork)
library(ggplot2)
library(reticulate)

setwd(dirname(getActiveDocumentContext()$path))

#import orthologs resolved at the level of hydrozoa
hydTab <- read.delim('Results_Sep15_1/Phylogenetic_Hierarchical_Orthogroups/N10Mod.tsv')

#drop any empty columns (species not in clade) 
hydTab <- hydTab[,!is.na(hydTab[1,])]

#drop AEP transcriptome (redundant)
hydTab <- hydTab[,-14]

#because cruxmelitensis has duplicated gene IDs we need to drop it from this analysis
hydTab <- hydTab[,-5]

#look just at orthologs within Hydra genus
hydraSpec <- c('H_circumcincta','H_oligactis','H_viridissima','H_vulgaris105','H_vulgarisAEP')

hydTab.hyd <- hydTab[,hydraSpec]

#find orthogroups that are completely absent from Hydra genus
hydTab.hyd.test <- apply(hydTab.hyd,1,function(x) length(which(x != "")) > 0)

hydTab.lost <- hydTab[!hydTab.hyd.test,]

#make sure the genes lost in Hydra are present in other hydrozoans (more than one)

hydTab.lost <- hydTab.lost[hydTab.lost$H_echinata != '' |
                             hydTab.lost$C_hemisphaerica != '',]

#specifically look at things lost in Hydra but retained in clytia

hydTab.lost.cl <- hydTab.lost[hydTab.lost$C_hemisphaerica != '',]

hydTab.lost.cl <- unlist(strsplit(hydTab.lost.cl$C_hemisphaerica, split = ", "))

write.table(hydTab.lost.cl,'clytiaLostInHydra.txt',row.names = F, col.names = F, quote = F)

####Lost Hydra gene expression in Clytia####

#import remapped clytia single cell data
cl <- readRDS('../ds/cl/remap/initSeurat.rds')

#convert neuronal subclustering object from clytia paper into seurat object
#only needs to be run once
use_condaenv("rScanpy", required = T)

sc <- import("scanpy")

adata <- sc$read_h5ad('../ds/cl/remap/neuron_subpops_fs.h5ad')

exprs <- t(adata$X)
colnames(exprs) <- adata$obs_names$to_list()
rownames(exprs) <- adata$var_names$to_list()
# Create the Seurat object
seurat <- CreateSeuratObject(exprs)
# Set the expression assay
seurat <- SetAssayData(seurat, "data", exprs)
# Add observation metadata
seurat <- AddMetaData(seurat, adata$obs)
# Add embedding
embedding <- adata$obsm["X_umap"]
rownames(embedding) <- adata$obs_names$to_list()
colnames(embedding) <- c("umap_1", "umap_2")
seurat[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")

#save converted object to save time later
saveRDS(seurat,'annotatedNeuroCl.rds')

#import annotations from preprint (also use their umap)
cl.annot <- readRDS('../ds/cl/remap/annotatedCl.rds')

cl.neuro <- readRDS('annotatedNeuroCl.rds')

#save the publication UMAP coords
cl.umap <- cl.annot@reductions$umap

#save publication metadata (includes cell identities)
cl.annot <- cl.annot@meta.data

#drop any cells not in the remapped dataset
cl.annot <- cl.annot[rownames(cl.annot) %in% colnames(cl),]

#make sure cell order is the same
cl.annot <- cl.annot[match(rownames(cl.annot),colnames(cl)),]

#only keep the cell annotations from the publication metadata
cl.annot <- cl.annot[,c('annos','annosSub')]

#drop whole-animal annotations of neuronal cells
cl.annot <- cl.annot[!(rownames(cl.annot) %in% colnames(cl.neuro)),]

#extract neuronal analysis metadata
cl.annot.n <- cl.neuro@meta.data

#drop any cells not in the remapped dataset
cl.annot.n <- cl.annot.n[rownames(cl.annot.n) %in% colnames(cl),]

cl.annot.n.rn <- rownames(cl.annot.n)

#subset metadata to just include neuron cell identities
cl.annot.n <- data.frame(annos='Neuron',annosSub=paste0('Neuro_',cl.annot.n$louvain_neur))

rownames(cl.annot.n) <- cl.annot.n.rn

#combine neuronal and whole-animal annotations
cl.annot <- rbind(cl.annot,cl.annot.n)

#make sure cell order in the new metadata table is the same as the seurat object 
cl.annot <- cl.annot[colnames(cl),]

#add cell type annotations to remapped dataset
cl@meta.data <- cbind(cl@meta.data,cl.annot[,c('annos','annosSub')])

#add publication umap to remapped data
cl@reductions$oldUmap <- cl.umap

DimPlot(cl,group.by='annosSub',reduction = 'oldUmap',label=T,repel = T) + NoLegend() + NoAxes()
ggsave('fullLabClUMAP.png',width = 9,height = 9,dpi = 300)

#convert underscores to dash for gene names (seurat formatting requirement)
hydTab.lost.cl.ds <- gsub('_','-',hydTab.lost.cl)

#drop genes that aren't in clytia seurat object
hydTab.lost.cl.ds <- hydTab.lost.cl.ds[hydTab.lost.cl.ds %in% rownames(cl)]

#calculate holistic score for expression of lost genes
cl <- AddModuleScore(cl,list(hydTab.lost.cl.ds),name='hyLost')

FeaturePlot(cl,'hyLost1',reduction = 'oldUmap') + NoAxes()
ggsave('clLostGeneScores.png',width = 9,height = 9,dpi = 300)

lostPlot <- cl@meta.data[,c(9,10)]

lostPlot$annosSub <- factor(lostPlot$annosSub,
                            levels = rev(c("i-Cells","Very Early Oocytes","Small Oocytes",
                                           "Medium Oocytes","Gland Cells-A","Gland Cells-B",
                                           "Gland Cells-C","Gland Cells-D","Gland Cells-E",
                                           "Neuro_0","Neuro_1","Neuro_2","Neuro_3","Neuro_4",
                                           "Neuro_5","Neuro_6","Neuro_7","Neuro_8","Neuro_9",
                                           "Neuro_10","Neuro_11","Neuro_12","Neuro_13","Neuro_14",
                                           "Nematocyte Precursors","Early Nematoblasts",
                                           "Mid Nematoblasts","Late Nematoblasts","Differentiating Nematocytes",
                                           "Terminal Differentiating Nematocytes","Mature Nematocytes",
                                           "Gonad Epidermis","Manubrium Epidermis","Exumbrella Epidermis",
                                           "Tentacle Epidermis","Striated Muscle of Subumbrella","Striated Muscle of Velum",
                                           "Radial Smooth Muscles","GastroDigestive-A","GastroDigestive-B",
                                           "GastroDigestive-C","GastroDigestive-D","GastroDigestive-E",
                                           "GastroDigestive-F","Endodermal Plate","Tentacle Bulb Distal Gastroderm",
                                           "Tentacle GFP Cells")))

ggplot(lostPlot,aes(y=annosSub,x=hyLost1,fill=annosSub)) + 
  geom_boxplot() + 
  #geom_jitter(height = 0, width = 0.1, size = 0.1,alpha=0.4) +
  theme_bw() +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 300, vjust = 0, hjust=0.05))
ggsave('lostScores.pdf',width = 6.5,height = 15)

saveRDS(cl,file='neuroLabGeneLossCl.rds')
