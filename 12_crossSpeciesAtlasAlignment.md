# Identifying Conserved Cell-Type-Specific Transcriptional Programs in *Hydra* and *Clytia*

This document covers our comparative analysis of the *Hydra* and *Clytia* single-cell atlases. This analysis entailed using reciprocal principal component analysis to align the two datasets, quantifying overall transcriptional similarities between cell types using mutual nearest neighbor analysis, and identifying orthologs with conserved expression patterns.

[TOC]



## Aligning the *Clytia* and *Hydra* Single-Cell Atlases

The process of aligning single cell datasets from different species uses the same approach as correcting for inter-sample batch effects. Essentially, we simply treat inter-species differences as (dramatic) batch effects. For this to be possible, the two datasets need to use the same gene IDs. In order to unambiguously link orthologous genes between *Hydra* and *Clytia*, we limited our alignment analysis to only genes that had a one-to-one orthology between the two species, as determined through our Orthofinder analysis. We then simply converted the gene IDs for the *Clytia* data to their *Hydra* equivalents. 

In the following code chunk, we identify one-to-one ortholog pairs, import and perform initial filtering on the *Clytia* and *Hydra* gene expression matrices, convert the *Clytia* gene IDs to their *Hydra* equivalents, and concatenate all of the separate samples (15 *Hydra* libraries and one *Clytia* library) into a list of Seurat objects. 

(*snippet from 01_alignment/aepCl.R*)

```R
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
```

We then performed standard dataset integration using reciprocal PCA analysis as implemented through Seurat. Following integration we performed PCA on the aligned dataset. 

(*snippet from 01_alignment/aepCl.R*)

```R
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
```

![integratedElbow](resources/integratedElbow.png)

We then performed Louvain clustering and UMAP dimensional reduction using 30 PCs. We then visualized the clustering results using the UMAP plot.

(*snippet from 01_alignment/aepCl.R*)

```R
#generate UMAP projection and perform louvain clustering
ndUse=1:30
aepCl.int <- RunUMAP(aepCl.int, reduction = "pca", dims = ndUse, min.dist = 0.18,spread = 0.2, seed.use = 4321)
aepCl.int <- FindNeighbors(aepCl.int, reduction = "pca", dims = ndUse)
aepCl.int <- FindClusters(aepCl.int, resolution = 0.4, graph.name = 'integrated_snn')

DimPlot(aepCl.int) + NoAxes() + NoLegend()
```

![initIntUMAP](resources/initIntegratedUMAP.png)

We also split the plot by species to see the contribution of each species to different cell populations

(*snippet from 01_alignment/aepCl.R*)

```R
aepCl.int@meta.data$species <- grepl('^D',aepCl.int@meta.data$orig.ident)

aepCl.int@meta.data$species[aepCl.int@meta.data$species == T] <- 'AEP'

aepCl.int@meta.data$species[aepCl.int@meta.data$species == 'FALSE'] <- 'Cl'

DimPlot(aepCl.int,group.by = 'species') + NoAxes()
```

![initUmapBySpec](resources/initIntUMAPbySpec.png)

(*snippet from 01_alignment/aepCl.R*)

```R
DimPlot(aepCl.int,split.by = 'species') + NoAxes() + NoLegend()
```

![initIntUMAPsplitSpec](resources/initIntUMAPsplitSpec.png)



To determine if our integrated dataset was properly aligning homologous cell types in the two species, we imported the cell annotations from the unintegrated atlases for each species. Specifically, we used the cell type annotations from our analysis of the AEP-mapped *Hydra* atlas (described in `05_hydraAtlasReMap.md`) and the annotations from the original *Clytia* atlas publication.

(*snippet from 01_alignment/aepCl.R*)

```R
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
```

![intUMAPcuratedIDbroad](resources/intUMAPcuratedIdentBroad.png)

(*snippet from 01_alignment/aepCl.R*)

```R
#plot using more granular cluster labels from clytia paper
#(a bit crowded)
DimPlot(aepCl.int,group.by = 'curatedIdentSub',split.by = 'species', label = T, repel = F, cols = clustPal,pt.size=0.5) + NoAxes() + NoLegend()
```

![intUMAPcuratedIDbroad](resources/intUMAPcuratedIdentFine.png)

(*snippet from 01_alignment/aepCl.R*)

```R
#a very course level of clustering for the integrated data just to group together broad cell types
aepCl.int <- FindClusters(aepCl.int, resolution = 0.05, graph.name = 'integrated_snn')
DimPlot(aepCl.int,group.by = 'integrated_snn_res.0.05',split.by = 'species', repel = F,pt.size=0.5) + NoAxes() + NoLegend()

#save object for later
saveRDS(aepCl.int,'aepClInt.rds')
```

![intUMAPbroad](resources/intUMAPbroad.png)

Based on these plots, our integration was successful, as homologous cell types colocalize in the UMAP and are clustered together by the Louvain algorithm.

## Quantifying Transcriptional Similarities Between *Clytia* and *Hydra* Cell Types

Although UMAPs can give an indication of similarity between different cells, it is not a quantitative way of assessing similarity. To assess similarity quantitatively, we adopted the alignment score metric proposed by Tarashansky et al. (2021). 

To generate this score, the 30 nearest inter-species neighbors for each cell in the aligned principal component space are identified. Then, for a given cell type in one species, the number of nearest cross-species neighbors that were from each cell type in the other species is tabulated. The final alignment score is defined as the portion of total cross-species pairs for a cell type in one species that belonged to a particular cell type in another species. Thus, a higher alignment score indicates that the two cell types being compared shared a higher number of neighbors in the aligned PC space. 

Because we were interested in identifying possible similarities between neuronal subtypes in *Clytia* and *Hydra*, and because the whole-animal version of the *Clytia* atlas did not resolve individual neuronal subtypes, we incorporated cell type labels from the *Clytia* neuronal subclustering analysis before performing the alignment quantification. 

To do this we first had to download the neuronal sub clustering object (available [here](https://data.caltech.edu/records/1804)) and convert it to a Seurat object.

(*snippet from 01_alignment/aepCl.R*)

```R
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
```

![clNeuroUMAP](resources/clNeuroUMAP.png)

We then updated the cell labels for *Clytia* neurons using the higher resolution cluster labels from the neuron sub-clustering

(*snippet from 01_alignment/aepCl.R*)

```R
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
```

We then extracted the principal component cell scores from the integrated seurat object and separated them by species. We then performed a mutual nearest neighbor analysis to identify the *Clytia* cells that was most similar to each *Hydra* cell. We then converted these cell IDs to cell type labels and tabulated the number of cells from each *Clytia* cell type associated with each *Hydra* cell type. Finally, we used these values to generate a Sankey diagram showing all alignment scores greater than 0.05 (indicating that greater than 5% of all nearest neighbors for a *Hydra* cell type were made up of a given *Clytia* cell type).

(*snippet from 01_alignment/aepCl.R*)

```R
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
```

<img src="resources/crossSankey.png" alt="crossSankey" style="zoom:33%;" />

As an alternative approach to visualize and evaluate the degree of transcriptional similarity between *Hydra* and *Clytia* cells, we also calculated a distance metric that captured how different the overall transcriptional profiles of cells from one species were when compared to their most similar transcriptional neighbors from the other species. To generate this metric, we calculated the average distance to the 30 nearest cross-species nearest neighbors in aligned principal component space for each cell. We then plotted these values on the cross-species UMAP to visually link these distance values to cell type annotations.

```R
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
```

<img src="resources/crossSpecDistance.png" alt="crossSpecDistance" style="zoom:33%;" />

```R
FeaturePlot(aepCl.int, 'dist', order = T, split.by = 'species', pt.size = 0.4)
```

![crossSpecDistanceSplit](resources/crossSpecDistanceSplit.png)

We also summarized the distance metric results using a boxplot that grouped distance scores by cell type in each species. The first boxplot presents distances for *Hydra* cell types.

```R
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
```

<img src="resources/hvDistBox.png" alt="hvDistBox" style="zoom:25%;" />

The second boxplot presents distances for *Clytia* cell types.

```R
ggplot(plotDf[plotDf$species == 'Cl',],aes(x=curatedIdentSub,y=dist,fill=curatedIdentSub)) + 
  geom_boxplot() + 
  #geom_jitter(height = 0, width = 0.1, size = 0.1,alpha=0.4) +
  theme_bw() +
  theme(legend.position="none") +
  coord_flip()
```

<img src="resources/clDistBox.png" alt="clDistBox" style="zoom:25%;" />

## Identifying Orthologous Genes with Conserved Expression Patterns

The alignment score we calculated above could be used to holistically examine transcriptional similarities, but it did not give us access to the genes whose expression patterns were conserved in *Hydra* and *Clytia*. To identify genes with conserved expression patterns, we needed a way to correlate the expression patterns in *Hydra* and *Clytia* cells. 

To do this, we used a very high resolution Louvain clustering analysis to generate ad-hoc 'pseudo-cells'. These pseudo-cells grouped together cells (regardless of species) that were close together in the aligned principal component space. We calculated average gene expression values in each species for each pseudo-cell and then, by matching pseudo-cell labels from each species, identified genes with similar expression patterns. 

We started by first importing the Seurat object containing the aligned *Clytia* and *Hydra* data and performing a high resolution Louvain clustering analysis. This generated a total of 132 'pseudo-cell' clusters.

(*03_expressionConservation/aepCl_Cor.R*)

```R
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
```

![pseudoCellUMAP](resources/pseudoCellUMAP.png)

We then extracted the normalized read matrix from the Seurat object, split cells both by species and pseudo-cell ID, and calculated average read counts. We then calculated correlation scores for pseudo-cell expression when comparing the two species. We used the resulting scores as a readout of how conserved the cell-type-specificity was for a pair of *Clytia* and *Hydra* orthologs.

(*03_expressionConservation/aepCl_Cor.R*)

```R
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
```

For identifying conserved expression patterns, we used high resolution clusters, which allowed us to perform a relatively more refined comparison of expression in the two species. To visualize these conserved expression patterns, we opted to use a gene-by-cell-type heatmap. In order to avoid having a high number of columns (cell types) in the heatmap, thus improving overall readability, we used a lower resolution clustering to summarize the data. We generated descriptive labels for these lower resolution clusters based on the identities of both the *Hydra* and *Clytia* cells making up each cluster, although we showed some preference for labels based on *Hydra* annotations as *Hydra* cell types are generally better characterized.

(*03_expressionConservation/aepCl_Cor.R*)

```R
#after identifying genes with similar expression patterns
#convey the data using a heatmap based on a lower resolution clustering
#otherwise it's just an insane number of columns

#import more descriptive names for lower resolution clusters
corPlotLabels <- read.csv('corPlotClustsLables.csv',header=F)

#make umap with descriptive names for clusters used for heatmap plots
aepCl.int$renameClust <- mapvalues(aepCl.int$integrated_snn_res.0.4,from = corPlotLabels$V1, corPlotLabels$V2, warn_missing = F)

DimPlot(aepCl.int,group.by = 'renameClust',label = T, repel = T) + NoLegend() + NoAxes()
```

![corPlotClusts](resources/corPlotClusts.png)

Using these broader clusters we again averaged gene expression by clusters (split by species) to generate the values that would populate the heatmap plot. We restricted these heatmaps to only those ortholog pairs with high correlation values (correlation score > 0.65)

We first plotted the *Hydra* data:

(*03_expressionConservation/aepCl_Cor.R*)

```R
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
```

<img src="resources/aepAllHeat.png" alt="aepAllHeat" style="zoom:25%;" />

We then plotted the *Clytia* data

(*03_expressionConservation/aepCl_Cor.R*)

```R
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
```

<img src="resources/clAllHeat.png" alt="clAllHeat" style="zoom:25%;" />

We also generated a table including correlation scores for all high scoring (correlation score > 0.65) ortholog pairs as well as functional annotations (putative vertebrate orthologs, protein domains, etc.)

(*03_expressionConservation/aepCl_Cor.R*)

````R
#generate table of correlation scores for export (only high scoring genes)

#bring in functional annotation information to make browsing results easier to interpret
annots <- read.csv('../../Orthofinder/HVAEP1_annotation.csv')

annots$H_vulgarisAEP <- t2g(annots$H_vulgarisAEP)

corRes.annot <- merge(corRes,annots,by.x = 'id', by.y = 'H_vulgarisAEP', all.x = T)

corRes.annot <- corRes.annot[order(-corRes.annot$cor),]

corRes.annot <- corRes.annot[corRes.annot$cor > 0.65,]

corRes.annot$clOrtho <- mapvalues(corRes.annot$id, from = orthos$H_vulgarisAEP, to = orthos$C_hemisphaerica, warn_missing = F)
write.csv(corRes.annot,'crossSpecExpCor.csv',row.names = F)
````

Finally, we generated another set of heat maps that included only putative transcription factors (predictions were based on interproscan results as described in `03_aepGenomeAnnotation.md`). First for AEP:

(*03_expressionConservation/aepCl_Cor.R*)

```R
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
```

![aepTfHeat](resources/aepTfHeat.png)

Then for *Clytia*:

(*03_expressionConservation/aepCl_Cor.R*)

```R
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
```

![clTfHeat](resources/clTfHeat.png)



## Files Associated with This Document



```
12_crossSpeciesAtlasAlignment/
├── 01_alignment
│   ├── aepClInt.rds
				R binary file containing the aligned Hydra and Clytia single-cell atlases
				formatted as a Seurat object.
│   └── aepCl.R
				R script that aligns the Clytia and Hydra single-cell atlases and quantifies
				and visualizes transcriptional similarities across homologous cell types using 
				a mutual nearest neighbor-based alignment score.
├── 02_similarityQuant
│   ├── aepClIntNeuroLab.rds
				R binary file containing a modified version of aepClInt.rds that includes
				neuronal subtype annotations from the original Clytia atlas publication.
│   ├── crossSpecDist.csv
				Table that provides the average distance to the 30 nearest cross-species
				nearest neighbors for every cell in the cross-species atlas. The table
				also includes the barcode and the curated cell type annotation for each
				cell.
│   ├── mnnRes.csv
				Table that lists the mutual nearest cross-species neighbors for all cells
				in the aligned Clytia and Hydra atlas.
│   ├── mnnSankeyLink.csv
				Table used for generating a sankey plot visualizing the alignment scores generated
				by systematically comparing Clytia and Hydra cell types in the aligned cross-species
				single-cell atlas. The table includes all cross-species cell type pairs that had an
				alignment score ≥ 0.05. In addition to specifying the two cell types, it also includes
				the alignment score under the 'value' column.
│   └── mnnSankeyNode.csv
				Table used for generating a sankey plot visualizing the alignment scores generated
				by systematically comparing Clytia and Hydra cell types in the aligned cross-species
				single-cell atlas. This table simply lists all cell types included in the alignment
        score sankey plot (required by the function that made the plot).
└── 03_expressionConservation
    ├── aepCl_Cor.R
    		R script that identifies genes with conserved expression patterns in the cross-
    		species atlas and visualizes the results using heatmaps.
    ├── corPlotClustsLables.csv
    		Table containing descriptive names for the Louvain clusters in the cross-species
        atlas used for generating the heatmaps visualizing conserved expression patterns.
    ├── corPlotClustsOrder.csv
    		Table that specifies the order that the clusters listed in corPlotClustsLables.csv
    		will be displayed in the heatmaps visualizing conserved expression patterns.
    ├── crossSpecExpCor.csv
    		Table listing all genes with a cross-species expression correlation score above 0.65.
    		The table also includes functional annotation data for each gene.
    ├── crossSpecMatClytia.csv
    		Table containing the values used to generate the heatmap depicting the expression 
    		patterns of all Clytia genes with a cross-species expression correlation score 
    		above 0.65.
    ├── crossSpecMatHydra.csv
    		Table containing the values used to generate the heatmap depicting the expression 
    		patterns of all Hydra genes with a cross-species expression correlation score 
    		above 0.65.
    ├── crossSpecTfMatClytia.csv
    		Table containing the values used to generate the heatmap depicting the expression 
    		patterns of all Clytia transcription factors with a cross-species expression 
    		correlation score above 0.65.
    ├── crossSpecTfMatHydra.csv
    		Table containing the values used to generate the heatmap depicting the expression
        patterns of all Hydra transcription factors with a cross-species expression 
        correlation score above 0.65.
    ├── geneCoreRes.csv
    		Table containing the cross-species expression correlation scores for all genes in 
    		the cross-species atlas.
    └── pcCellIDs.csv
    		Table containing the pseudo-cell assignments for all cells in the cross-species 
    		atlas.
```

