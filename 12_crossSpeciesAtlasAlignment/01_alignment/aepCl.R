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
library(fields)
library(BiocNeighbors)
library(networkD3)
library(reticulate)

#utility to convert transcript ID to gene ID
t2g <- function(x){
  vapply(x, function(y) gsub('HVAEP1_T(\\d+)[.]\\d','HVAEP1-G\\1',y),"")
}

setwd(dirname(getActiveDocumentContext()$path))

#import list of non-doublet Hydra cell IDs
nonDub <- read.delim('../dropSeqMapping/finalNonDub.txt')[,1]

#import orthofinder results to identify one-to-one orthologs between clytia and hydra
orthos <- read.delim('../../Orthofinder/Results_Sep15_1/Phylogenetic_Hierarchical_Orthogroups/N14.tsv')

orthos <- orthos[,c('C_hemisphaerica','H_vulgarisAEP')]

#drop genes without orthologs in the other species
orthos <- orthos[orthos$C_hemisphaerica != '' & orthos$H_vulgarisAEP != '',]

#drop genes that don't have one-to-one orthology
orthos <- orthos[!(grepl(',',orthos$C_hemisphaerica) | grepl(',',orthos$H_vulgarisAEP)),]

orthos$H_vulgarisAEP <- gsub('HVAEP1_T(\\d+)[.].*','HVAEP1_G\\1',orthos$H_vulgarisAEP)

#import data and convert
#get list of hydra count matrix files (include path in name)
readMats <- list.files(path = '../dropSeqMapping', pattern = 'dge.txt.gz', recursive = T, full.names = T)

#read in each individual read matrix, generate a separate seurat object for each, and do some initial filtering
inAep <- lapply(readMats, function(x) {
  y <- read.delim(x,stringsAsFactors = F,header = T, row.names = 1)
  
  #get batch name from parent folder of read matrix
  libName <- gsub('[.]+\\/dropSeqMapping\\/(D\\d+[^\\/]+)\\/.*','\\1',x)
  
  #apply batch names to cell barcodes to prevent redundant barcodes
  colnames(y) <- paste0(colnames(y),'-',libName)
  
  #drop doublet cells
  y <- y[,colnames(y) %in% nonDub]
  
  #only keep genes that have one-to-one orthology with clytia genes
  y <- y[rownames(y) %in% t2g(orthos$H_vulgarisAEP),]
  
  #initialize seurat object
  tmpDS <- CreateSeuratObject(counts = y, project = libName, min.cells = 3, min.features = 200)
  
  #perform preliminary filtering
  subset(tmpDS, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 & nCount_RNA > 500 & nCount_RNA < 75000)
})

#import raw count matrix for remapped clytia data
cl <- Read10X_h5('../cl/remap/raw_feature_bc_matrix.h5')

#drop genes without one-to-one orthology with Hydra genes
cl <- cl[rownames(cl) %in% orthos$C_hemisphaerica,]

#convert clytia gene names to their hydra equivalents
rownames(cl) <- mapvalues(rownames(cl),from = orthos$C_hemisphaerica, to = orthos$H_vulgarisAEP, warn_missing = F)

#import the original annotated clytia dataset
clPP <- readRDS('../cl/remap/annotatedCl.rds')

#remove any clytia cells that weren't in the original publication
cl <- cl[,colnames(cl) %in% rownames(clPP@meta.data)]

rm(clPP)

cl <- CreateSeuratObject(counts = cl, project = 'clDS', min.cells = 3, min.features = 200)

cl <- subset(cl, nFeature_RNA < 4000 & nCount_RNA > 500 & nCount_RNA < 100000)

#combine clytia object with list of hydra objects
aepCl <- inAep

aepCl <- append(aepCl,cl)

rm(inAep)


#normalize data using SCTransform
aepCl <- lapply(aepCl, FUN = SCTransform, method = "glmGamPoi")

#integrate data (both across batches and across species) using reciprocal pca
features <- SelectIntegrationFeatures(aepCl,assay=rep('SCT',length(aepCl)))

aepCl <- PrepSCTIntegration(object.list = aepCl, anchor.features = features,assay=rep('SCT',length(aepCl)))

aepCl <- lapply(X = aepCl, FUN = RunPCA, features = features, verbose = F, assay = 'SCT', npcs = 50)

dsAnchors <- FindIntegrationAnchors(object.list = aepCl,
                                    normalization.method = "SCT",
                                    anchor.features = features,
                                    dims = 1:50,
                                    reduction = "rpca",
                                    k.anchor = 20)

rm(aepCl)
gc()

aepCl.int <- IntegrateData(anchorset = dsAnchors, normalization.method = "SCT", dims = 1:50)

rm(dsAnchors)
gc()

#determining the dimensionality of the data
aepCl.int <- RunPCA(aepCl.int, npcs = 80,verbose=F)

ElbowPlot(aepCl.int,ndims = 80)

#generate UMAP projection and perform louvain clustering
ndUse=1:30
aepCl.int <- RunUMAP(aepCl.int, reduction = "pca", dims = ndUse, min.dist = 0.18,spread = 0.2, seed.use = 4321)
aepCl.int <- FindNeighbors(aepCl.int, reduction = "pca", dims = ndUse)
aepCl.int <- FindClusters(aepCl.int, resolution = 0.4, graph.name = 'integrated_snn')

DimPlot(aepCl.int) + NoAxes() + NoLegend()

aepCl.int@meta.data$species <- grepl('^D',aepCl.int@meta.data$orig.ident)

aepCl.int@meta.data$species[aepCl.int@meta.data$species == T] <- 'AEP'

aepCl.int@meta.data$species[aepCl.int@meta.data$species == 'FALSE'] <- 'Cl'

DimPlot(aepCl.int,group.by = 'species') + NoAxes()

DimPlot(aepCl.int,split.by = 'species') + NoAxes() + NoLegend()

####Bring in curated cell identity annotations####
clPP <- readRDS('../cl/remap/annotatedCl.rds')

ds <- readRDS('../dropSeqMapping/nonDubLabeledSeurat.rds')

annotTab <- rbind(data.frame(ID=colnames(clPP), curatedIdent = as.character(clPP@meta.data$annos),curatedIdentSub=as.character(clPP@meta.data$annosSub)),
                  data.frame(ID=colnames(ds), curatedIdent = as.character(ds@meta.data$curatedIdent),curatedIdentSub = as.character(ds@meta.data$curatedIdent)))

#only keep cells that made it into the integrated object
annotTab <- annotTab[annotTab$ID %in% colnames(aepCl.int),]

#reorder cells so they are the same order as the integrated metadata table
rownames(annotTab) <- annotTab$ID
annotTab <- annotTab[colnames(aepCl.int),]

aepCl.int@meta.data$curatedIdent <- annotTab$curatedIdent

aepCl.int@meta.data$curatedIdentSub <- annotTab$curatedIdentSub

#colors for cluster plot
clustPal <- c("#45adec","#cbee2f","#d565f5","#85e538","#f452cd",
              "#57e958","#997cfb","#eadc25","#528efb","#abd533",
              "#d07be8","#65c737","#a585f0","#54e170","#e27cd1",
              "#4eb13a","#f2699f","#3ee698","#f76748","#3be6ea",
              "#f57417","#3295e9","#d0cc36","#6e8de9","#def469",
              "#7d9af7","#95ca42","#47a2f7","#e8b027","#2499d7",
              "#e8802b","#9d94e5","#a1f078","#db93d2","#83ce5a",
              "#f2686c","#51dbb4","#e77f4a","#56c88d","#e87c86",
              "#6bc565","#e8886e","#82e596","#c98b23","#a6f4b4",
              "#c78b3c","#5ab267","#eaaa43","#99cd86","#df9f64",
              "#ade486","#e8aa59","#91a824","#d6e58f","#b39420",
              "#cce679","#b0953a","#efda63","#8eac57","#d8b644",
              "#94af45","#dcbd61","#aaa126","#ceca7b","#c9ca4d",
              "#a4a154","#e8de7a","#afa746")

set.seed(5321)
clustPal <- clustPal[sample(length(clustPal))]

#plot using broad cluster labels from clytia paper
DimPlot(aepCl.int,group.by = 'curatedIdent',split.by = 'species', label = T, repel = T, cols = clustPal,pt.size=0.5) + NoAxes() + NoLegend()
ggsave('curatedIdentBySpec.pdf',width = 22,height = 10)

#plot using more granular cluster labels from clytia paper
#(a bit crowded)
DimPlot(aepCl.int,group.by = 'curatedIdentSub',split.by = 'species', label = T, repel = F, cols = clustPal,pt.size=0.5) + NoAxes() + NoLegend()
ggsave('curatedIdentBySpecSub.pdf',width = 30,height = 10)

#a very course level of clustering for the integrated data just to group together broad cell types
aepCl.int <- FindClusters(aepCl.int, resolution = 0.05, graph.name = 'integrated_snn')
DimPlot(aepCl.int,group.by = 'integrated_snn_res.0.05',split.by = 'species', repel = F,pt.size=0.5) + NoAxes() + NoLegend()
ggsave('intClust.pdf',width = 10,height = 10)

#save object for later
saveRDS(aepCl.int,'aepClInt.rds')
aepCl.int <- readRDS('aepClInt.rds')

#####cross-species cell type alignment quantification####

#alignment score
#the average number of mutual nearest cross-species neighbors of 
#each cell relative to the maximum possible number of neighbors.

#convert neuronal subclustering object from clytia paper into seurat object
#only needs to be run once
use_condaenv("rScanpy", required = T)

sc <- import("scanpy")

adata <- sc$read_h5ad('../cl/remap/neuron_subpops_fs.h5ad')

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

clN <- seurat
rm(seurat)

clN <- readRDS('annotatedNeuroCl.rds')

DimPlot(clN, group.by = 'louvain_neur', label = T) + NoAxes() + NoLegend()

#for neuronal cells, replace the cell labels from the
#from the whole-animal analysis with the labels
#from the neuronal subcluster analysis
annotTab <- annotTab[!(annotTab$ID %in% colnames(clN)),]

annotTab <- rbind(data.frame(ID=colnames(clN), curatedIdent = 'Neuron',curatedIdentSub=paste0('neuro_',as.character(clN$louvain_neur))),annotTab)

annotTab <- annotTab[annotTab$ID %in% colnames(aepCl.int),]

rownames(annotTab) <- annotTab$ID

annotTab <- annotTab[colnames(aepCl.int),]

aepCl.int@meta.data$curatedIdent <- annotTab$curatedIdent

aepCl.int@meta.data$curatedIdentSub <- annotTab$curatedIdentSub

DefaultAssay(aepCl.int) <- 'integrated'

#extract cell score for principal components from seurat object
pcSpace <- aepCl.int@reductions$pca@cell.embeddings[,1:80]

#split cell PC scores by species
pcSpace.a <- pcSpace[aepCl.int$species == 'AEP',]

pcSpace.c <- pcSpace[aepCl.int$species == 'Cl',]

#specify the number of nearest neighbors to use for MNN
kUse <- 30

#find nearest neighbors for all cells in aligned PC space
mnnRes <- findMutualNN(pcSpace.a,pcSpace.c, k1 = kUse)

#get the cross-species neighbor pair cell IDs
mnnRes <- data.frame(AEP=mnnRes[[1]],Cl=mnnRes[[2]])

mnnRes$AEP <- rownames(pcSpace.a)[mnnRes$AEP]
mnnRes$Cl <- rownames(pcSpace.c)[mnnRes$Cl]

mnnRes.orig <- mnnRes

#convert neighbor pair cell IDs to cell type labels
mnnRes$Aepclust <- mapvalues(mnnRes$AEP,from = rownames(aepCl.int@meta.data), to = aepCl.int$curatedIdentSub, warn_missing = F)
mnnRes$ClClust <- mapvalues(mnnRes$Cl,from = rownames(aepCl.int@meta.data), to = aepCl.int$curatedIdentSub, warn_missing = F)

write.csv(mnnRes,'mnnRes.csv')

#split pair assignments by AEP cell
mnnRes.a <- split(mnnRes,  mnnRes$Aepclust)

#for a given hydra cell type, calculate the portion of
#it's neighbors made up by each clytia cell type
mnnRes.a <- lapply(mnnRes.a, function(x) {
  res <- as.data.frame(table(x$ClClust))
  res$score <- (res$Freq/length(unique(x$AEP)))/kUse
  return(res)
})

#set a minimum score cutoff of 0.05
mnnRes.a <- lapply(mnnRes.a,function(x){
  x[x$score >= 0.05,]
})

#drop cell types without any meaningful alignment
mnnRes.a <- mnnRes.a[sapply(mnnRes.a,nrow) > 0]

#add a column with the hydra cell ID
#(to keep track of stuff when the list is collapsed)
mnnRes.a <- lapply(1:length(mnnRes.a), function(x){
  df <- mnnRes.a[[x]]
  df$aepClust <- names(mnnRes.a)[x]
  return(df)
})

mnnRes.a <- do.call(rbind,mnnRes.a)

#these objects are formatted in a particular way to fit the input requirements
#of the sankey plotting funciton
#this object simply lists all of the cell types in the analysis (both clytia and hydra)
nodeDf <- data.frame(name=unique(c(as.character(mnnRes.a$Var1),mnnRes.a$aepClust)))

#need to manually set the order so that the resulting plot has a logical order
nodeDf <- nodeDf[c(1,3,4,5,2,12,13,14,11,16,15,17,22,10,6,20,8,7,9,21,18,19,23,
                   24,27,26,28,30,29,46,25,34,36,35,37,39,38,47,41,33,43,48,44,31,32,45,49,40,42,50),,drop=F]

#this object specifies the cell type pairs and the alignment score between them
linkDf <- data.frame(source = mnnRes.a$Var1, target=mnnRes.a$aepClust, value=mnnRes.a$score)

#these encode the index for the cell type in the nodeDF object (zero indexed)
linkDf$IDsource <- match(linkDf$source, nodeDf$name)-1 
linkDf$IDtarget <- match(linkDf$target, nodeDf$name)-1

#plot using plotly sankey plot
fig <- plot_ly(
  type = "sankey",
  orientation = "h",
  node = list(
    label = nodeDf$name,
    color = clustPal[1:48],
    pad = 20,
    thickness = 20,
    line = list(
      color = "black",
      width = 0.5
    )
  ),
  link = list(
    source = linkDf$IDsource,
    target = linkDf$IDtarget,
    value =  linkDf$value
  )
)

fig

orca(fig, "crossSankey.pdf",width=500,height=800)

write.csv(nodeDf,'mnnSankeyNode.csv')
write.csv(linkDf,'mnnSankeyLink.csv')

saveRDS(aepCl.int,'aepClIntNeuroLab.rds')
aepCl.int <- readRDS('aepClIntNeuroLab.rds')

#get distance between mnn pairs

#make distance matrix
pcSpace.a <- pcSpace[aepCl.int$species == 'AEP',1:80]

pcSpace.c <- pcSpace[aepCl.int$species == 'Cl',1:80]

cDist <- rdist(pcSpace.a,pcSpace.c)

rownames(cDist) <- rownames(pcSpace.a)

colnames(cDist) <- rownames(pcSpace.c)

#for every hydra cell, get the average distance for the 30 closest clytia cells
aepDist <- vapply(rownames(pcSpace.a), function(x){
  dists <- cDist[x,]
  dists <- dists[order(dists)]
  dists <- mean(dists[1:30])
  return(dists)
}, numeric(1))

#for every clytia cell, get the average distance for the 30 closest hydra cells
clDist <- vapply(rownames(pcSpace.c), function(x){
  dists <- cDist[,x]
  dists <- dists[order(dists)]
  dists <- mean(dists[1:30])
  return(dists)
}, numeric(1))

dist.df <- data.frame(cellID = c(rownames(pcSpace.a),rownames(pcSpace.c)), dist = c(aepDist,clDist))

dist.df <- dist.df[match(colnames(aepCl.int),dist.df$cellID),]

dist.df$curatedIdent <- aepCl.int$curatedIdentSub

write.csv(dist.df, file = 'crossSpecDist.csv',row.names = F)

aepCl.int$dist <- dist.df$dist

FeaturePlot(aepCl.int, 'dist',order = T)
ggsave('crossSpecDistance.pdf',width = 5.5, height = 5)
FeaturePlot(aepCl.int, 'dist', order = T, split.by = 'species', pt.size = 0.4)
ggsave('crossSpecDistanceSplit.pdf',width = 10, height = 5)

plotDf <- aepCl.int@meta.data

plotDf$curatedIdentSub <- factor(plotDf$curatedIdentSub)

plotFacts <- levels(plotDf$curatedIdentSub)

factOrder <- c("I_ISC","I_FemGC","I_MaleGC",
  "I_GlProgen","I_ZymoGl","I_GranGl",
  "I_SpumMucGl","I_Neuro","I_Ec1/5N",
  "I_Ec1N","I_Ec2N","I_Ec3N",
  "I_Ec4N","I_En1N","I_En2N",
  "I_En3N","I_EarlyNem","I_DesmoNB",
  "I_StenoNB","I_IsoNB","I_DesmoNC",
  "I_StenoNC","I_IsoNC","Ec_BasalDisk",
  "Ec_Peduncle","Ec_BodyCol/SC","Ec_Head",
  "Ec_Tentacle","En_Foot","En_BodyCol/SC",
  "En_Head","En_Tentacle",
  "i-Cells","Very Early Oocytes","Small Oocytes",
  "Medium Oocytes","Gland Cells-A","Gland Cells-B",
  "Gland Cells-C","Gland Cells-D","Gland Cells-E",
  "neuro_0","neuro_1","neuro_2",
  "neuro_3","neuro_4","neuro_5",
  "neuro_6","neuro_7","neuro_8",
  "neuro_9","neuro_10","neuro_11",
  "neuro_12","neuro_13","neuro_14",
  "Nematocyte Precursors","Early Nematoblasts","Mid Nematoblasts",
  "Late Nematoblasts","Differentiating Nematocytes","Terminal Differentiating Nematocytes",
  "Mature Nematocytes","Exumbrella Epidermis","Manubrium Epidermis",
  "Gonad Epidermis","Tentacle Epidermis","Tentacle GFP Cells",
  "Radial Smooth Muscles","Striated Muscle of Subumbrella","Striated Muscle of Velum",
  "Endodermal Plate","GastroDigestive-A","GastroDigestive-B",
  "GastroDigestive-C","GastroDigestive-D","GastroDigestive-E",
  "GastroDigestive-F","Tentacle Bulb Distal Gastroderm")

plotDf$curatedIdentSub <- factor(plotDf$curatedIdentSub, levels = factOrder)

ggplot(plotDf[plotDf$species == 'AEP',],aes(x=curatedIdentSub,y=dist,fill=curatedIdentSub)) + 
  geom_boxplot() + 
  #geom_jitter(height = 0, width = 0.1, size = 0.1,alpha=0.4) +
  theme_bw() +
  theme(legend.position="none") +
  coord_flip()
ggsave('hvDistBox.pdf',width = 5, height = 15)

ggplot(plotDf[plotDf$species == 'Cl',],aes(x=curatedIdentSub,y=dist,fill=curatedIdentSub)) + 
  geom_boxplot() + 
  #geom_jitter(height = 0, width = 0.1, size = 0.1,alpha=0.4) +
  theme_bw() +
  theme(legend.position="none") +
  coord_flip()
ggsave('clDistBox.pdf',width = 6, height = 15)



