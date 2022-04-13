# Identifying Candidate Regulators of Cell-Type-Specific Transcription in *Hydra*

This document covers our analysis aimed at identifying transcriptional regulators of gene co-expression modules in the *Hydra* single-cell atlas. This entailed integrating our non-negative matrix factorization (NMF) results (described in `05_hydraAtlasReMap.md`) with our transcription factor binding site conservation analysis (described in `07_genomeConservation.md`) to identify motifs enriched in the putative regulatory sequences of sets of co-expressed genes (i.e., metagenes). Following the motif enrichment analysis, we then looked for transcription factors whose expression pattern and sequence binding preferences correlated with motif enrichement patterns. 

[TOC]

## Preparatory Motif and Transcription Factor Analyses

Prior to the main analysis, several supplemental files needed to be generated. These included annotations linking conserved motifs to their most proximal gene, information on motif sequence similarity (to reduce redundancy in the enrichment results), the names of the transcription factor associated with each binding motif, and predictions of which gene models in the AEP assembly are transcription factors based on their pfam domain annotations. 

### Linking Predicted Transcription Factor Sites to Their Nearest Genes 

First, we used UROPA to link conserved transcription factor binding sites to their putative target genes based on proximity. Specifically, we linked a motif to the nearest transcription start site as long as it fell within 30 Kb. We did this for both the binding site predictions based on the bona fide JASPAR motif sequences as well as the predictions based on shuffled motif sequences (both generated in `07_genomeConservation.md`).

We used the following config file to annotate the binding site predictions for the bona fide motif sequences (`conMotsATAC.bed`):

(*01_motifPrep/conMotATACAnnot.json*)

```json
{
"queries":[
    {"distance": [30000, 30000], "feature": ["gene"], "feature_anchor": ["start"], "name": "query_1", "relative_location": ["Downstream", "Upstream", "OverlapStart", "FeatureInsidePeak", "PeakInsideFeature", "OverlapEnd"], "strand": "ignore"}
          ],
"show_attributes": ["gene_id"],
"priority": "False",
"gtf": "HVAEP1.GeneModels.longestIso.gtf",
"bed": "conMotsATAC.bed",
"prefix": "conMotATAC",
"outdir": "motifAnnot",
"threads": "1",
"output_by_query": "False"
}
```

We used the following config file to annotate the binding site predictions for the shuffled motif predictions (`conShufMotsATAC.bed`)

(*01_motifPrep/conShufMotATACAnnot.json*)

```json
{
"queries":[
    {"distance": [30000, 30000], "feature": ["gene"], "feature_anchor": ["start"], "name": "query_1", "relative_location": ["Downstream", "Upstream", "OverlapStart", "FeatureInsidePeak", "PeakInsideFeature", "OverlapEnd"], "strand": "ignore"}
          ],
"show_attributes": ["gene_id"],
"priority": "False",
"gtf": "HVAEP1.GeneModels.longestIso.gtf",
"bed": "conShufMotsATAC.bed",
"prefix": "conShufMotsATAC",
"outdir": "motifAnnot",
"threads": "1",
"output_by_query": "False"
}
```

We then ran the two annotation analyses with the following commands:

```bash
uropa -i conMotATACAnnot.json 
uropa -i conShufMotATACAnnot.json
```

### Clustering Motifs to Reduce Redundancy

Closely related transcription factors that use the same DNA binding domain tend to have highly similar (or virtually identical) binding preferences. Thus, motif enrichment analyses tend to have many redundant motifs in their results. To identify redundant motifs in our compiled JASPAR database, we used two clustering analyses: one that grouped motifs purely based on their sequence composition and another that grouped motifs based on their enrichment patterns in the single-cell atlas. In this section, we describe the former clustering approach (the latter approach is described later). 

To perform the sequence-based clustering, we chose to make use of the `compareMotifs.pl` utility script provided as part of the HOMER suite of motif enrichment tools. In order to use this script, we needed to reformat our JASPAR-formated motifs into the HOMER motif format. We did this using a utility script also provided as part of HOMER:

```bash
perl parseJasparMatrix.pl pooledJasparNR.txt > pooledJasparNR.homer.motifs
```

We then used the `compareMotifs.pl` to generate a matrix of pairwise similarity scores for all motifs in our custom JASPAR database

(*01_motifPrep/homerMotCompare.sh*)

```bash
#!/bin/bash

perl compareMotifs.pl pooledJasparNR.homer.motifs motCompare -matrix motComp.txt
```

Within a custom R script, we used these similarity scores for a hierarchical clustering analysis that grouped together motifs with similar sequence composition. These results were later combined with a second hierarchical clustering analysis (described in the "Single-Cell Motif Enrichment Analysis" section) to arrive at the final motif cluster assignments.

(*01_motifPrep/clusterMots.R*)

```R
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

#import similarity matrix for JASPAR animal PWMs
similarityMat <- read.delim("motComp.txt", row.names = 1, header = T, sep = "\t")

#reformat motif names (rows and columns)
newNames <- rownames(similarityMat)
newNames <- gsub('/Jaspar','',newNames)
newNames <- gsub('/','_',newNames)

rownames(similarityMat) <- newNames
colnames(similarityMat) <- newNames

#generate distance matrix from correlation coefficient
d <- as.dist(1 - similarityMat)

#use hierarchical clustering to group motifs based on correlation
hc1 <- hclust(d, method = "average")

pdf('motDistClust.pdf',width=150,height=20)
plot(hc1)
rect.hclust(hc1 , h = 0.3)
abline(h = 0.3, col = 'red')
dev.off()

#cut the branches at 0.3 and group accordingly
clusts <- cutree(hc1, h = 0.3)
length(unique(clusts))

#save cluster information
repMotifs <- data.frame(ID = names(clusts), clust = clusts)

write.csv(repMotifs, file = "motif_clusters.csv")
```

Below is the tree generated by the hierarchical clustering analysis used to generate the `motif_clusters.csv` file:

![motClust](resources/motDistClust.png)



### Annotate Transcription Factors in the Genome Gene Models

Ultimately, the goal of the analysis was to assign transcription factors as candidate regulators of co-expressed genes. This required that we determine which gene models in the AEP assembly are likely to be transcription factors. To do this, we made use of our InterProScan results (described in `03_aepGenomeAnnotation.md`). Specifically, we used a manually curated list of protein domains and gene ontology terms to generate a candidate list of transcription factors.

(*01_motifPrep/tfList.R*)

```R
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

#import ipr results
ipr <- read.delim('HVAEP1.prot.longestIso.fa.tsv', header = F)

#pull initial list of candidate TFs based on either domain composition or GO annotation
ipr.tf <- ipr[grepl('GO:0003700|GO:0003677|GO:0006355|GO:0043565',ipr$V14) | grepl('helix-loop-helix|hmg|sox|winged|c2h2',ipr$V6,ignore.case = T),]

#drop genes associated with non-TF GO terms
exclIDs <- ipr[grepl('GO:0000814|GO:0008180|GO:0005840|GO:0031213|GO:0032021|GO:0006281|GO:0000150|GO:0003721|GO:0030337|GO:0006298|GO:0000723|GO:0006302|GO:0003899|GO:0006351|GO:0006313|GO:0003887|GO:0000786|GO:0003697|GO:0006260|GO:0016592|GO:0003755|GO:0006265|GO:0003684|GO:0006338|GO:0006384|GO:0006383|GO:0003723|GO:0006270|GO:0003910|GO:0006367|GO:0032508|GO:0006289',ipr$V14),1]

#drop genes associated with non-TF annotations
moreEclIDs <- ipr[grepl('Oxygenase|transferase|Transposase|phosphodiesterase|peptidase|DNA REPLICATION|Tetratricopeptide|membrane-bound|hydrolase|proteinase|mitochondrial|reductase|synthase|autointegration|nuclease|TRANSFERASE|histone|helicase|topoisomerase|DNA polymerase|rejoining|rad51|Membrane-anchored|Translin|BESS|seet_6|LIPOMA|QSOX',ipr$V6, ignore.case = T),1]

ipr.tf <- ipr.tf[!(ipr.tf$V1 %in% c(exclIDs,moreEclIDs)),]

#subset to just the IDs
tfIDs <- unique(ipr.tf$V1)

#convert transcript IDs to gene IDs
tfIDs <- gsub('HVAEP1_T(\\d+)[.]\\d+','HVAEP1_G\\1',tfIDs)

write.csv(ipr.tf,'fullTF.csv',row.names=F)
write.table(tfIDs,'tfIDs.txt',row.names = F,col.names = F,quote = F)
```

Based on this analysis, we identified a total of 811 candidate transcription factors in the AEP gene models. 

### Download Motif Metadata

In order to determine which transcription factors could bind which JASPAR motifs, we made use of the uniprot database protein domain annotations available for each JASPAR motif. JASPAR provides uniprot IDs for each motif in its database, and each uniprot entry contains Pfam domain annotations. This allowed us to determine the pfam domains associated with each motif in our motif database. These domains could then be linked to the domain predictions from our InterProScan analysis to associate specific gene models with JASPAR motifs. These gene/motif links provided the basis for predicting regulators of gene co-expression.

To access the metadata for each JASPAR motif, we first had to compile a list of all JASPAR motif IDs of interest:

```bash
grep '>' pooledJasparNR.txt | sed 's/>\([^\t]\+\).*/\1/g' > jasparIDs.txt
```

We then used both the JASPAR and uniprot REST APIs to download the Pfam domains associated with each JASPAR motif. We then cross-referenced these domains with the domain composition of our list of predicted transcription factor gene models to generate lists of genes that could plausibly bind to each binding motif. 

(*01_motifPrep/protDBLink.R*)

```bash
library(rstudioapi)
library(jsonlite)
library(httr)
library(RCurl)
setwd(dirname(getActiveDocumentContext()$path))

#get full list of jaspar motif IDs
motif.IDs <- read.delim('jasparIDs.txt')[,1]

results <- list(length(motif.IDs))

#get TF name and family for each motif ID
for (i in 1:length(motif.IDs)) {
  print(i)
  link <- paste0("http://jaspar2020.genereg.net/api/v1/matrix/",as.character(motif.IDs[i]),"/")
  result <- fromJSON(url(link))
  result.n <- paste(result[[c("name")]],collapse = ', ')
  result.f <- paste(result[[c("family")]],collapse = ', ')
  results[[i]] <- c(result.n,result.f)
}

#collapse results into DF
motInfo <- data.frame(ID=motif.IDs,
                      name=vapply(results,function(x) x[1],""),
                      family=vapply(results,function(x) x[2],""))

write.csv(motInfo,'motifInfo.csv')

#get uniprot IDs for each motif from JASPAR
results <- list(length(motif.IDs))

for (i in 1:length(motif.IDs)) {
  print(i)
  link <- paste0("http://jaspar2020.genereg.net/api/v1/matrix/",as.character(motif.IDs[i]),"/")
  result <- fromJSON(url(link))
  result <- result[["uniprot_ids"]]
  if(length(result != 0)) {
    results[[i]] <- result
  } else {
    print("Empty")
  }
}

#pull names for each result
names(results) <- motif.IDs

pfamAnnot <- list()
#Pull the pfam domains assocaited with each swissprot entry
for (i in 1:length(results)) {
  print(i)
  subRes <- results[[i]]
  enGene <- c()
  for(ID in subRes){
    print(ID)
    link <- paste0("https://www.uniprot.org/uniprot/",ID,".txt")
    result <- GET(link)
    result <- rawToChar(result$content)
    result <- strsplit(result, "\n")
    result <- result[[1]]
    result <- result[grepl("Pfam", result)]
    result <- strsplit(result,'; ')
    result <- lapply(result,function(x) x[2])
    result <- unlist(result)
    pfamAnnot[[i]] <- unique(c(enGene,result))
  }
}

names(pfamAnnot) <- names(results)

#pull ipr results
ipr <- read.delim('../../Genome_annotation/functionalAnnotation/HVAEP1.prot.longestIso.fa.tsv', header = F)
ipr$V1 <- gsub('_T(\\d+)[.]\\d+','_G\\1',ipr$V1)

#import TF IDs
tfID <- read.delim('../../Genome_annotation/functionalAnnotation/tfIDs.txt',header=F)[,1]

#subset ipr table to include only TFs
ipr.tf <- ipr[ipr$V1 %in% tfID,]

#get PF domains associated with all the TF entries in the ipr table
ipr.tf.pf <- unique(ipr.tf[ipr.tf$V4 == 'Pfam',5])

#drop pfam domains from swissprot that didn't show up in our ipr domain list from our hydra TFs
pfamAnnot <- lapply(pfamAnnot,function(x) x[x %in% ipr.tf.pf])

#get the Hydra genes that contain the domains linked to each motif
pfamAnnot.genes <- lapply(pfamAnnot,function(x) unique(ipr[ipr$V5 %in% x,1]))

#drop motifs without a linked gene
pfamAnnot.genes <- pfamAnnot.genes[sapply(pfamAnnot.genes,length) > 0]

saveRDS(pfamAnnot.genes,'motifBindPfam.rds')

```

## Single-Cell Motif Enrichment Analysis

### Analysis Using Bona Fide Binding Motifs

We performed the main motif enrichment and gene co-expression regulator predictions within a single script (`hydraRegulators.R`). We have opted to break this script into smaller code blocks in this document to facilitate explanation.

#### Performing Gene Set Enrichment Analyses

After the initial steps of setting up the environment for the analysis, our first step was to generate gene sets. These gene sets, which group genes by the motifs present in their presumptive regulatory regions, served as the basis for a gene set enrichment analysis (GSEA; described in the next section). Essentially, the GSEA analysis determined if the presence of a particular binding motif in a gene's regulatory sequence was positively correlated with the gene having a higher NMF gene score for a particular metagene (the metagene gene score reflects how well a gene's expression pattern mimics the metagene expression pattern).

(*snippet from 02_enrichmentAnalysis/hydraRegulators.R*)

```R
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
```

To perform the motif enrichment analysis we used the `fgseaMultilevel` function from the `fgsea` package. This function requires a list of gene sets and a named vector of scores. For our analysis, the gene sets list was a list made up of vectors of gene name, with the name for each vector being the motif ID associated with the gene IDs contained in that vector. The named vector of scores correspond to the gene scores for a specific metagene. 

For our analysis, we performed gsea iteratively for each of the 56 metagenes we identified in the *Hydra* single-cell atlas. We then dropped (i.e., converted to 0) any enrichment score that failed to pass our significance cutoff (adjusted p-value ≤ 0.01). We then combined all the enrichment results into a single data frame.

(*snippet from 02_enrichmentAnalysis/hydraRegulators.R*)

```R
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
```

#### Generating Single-Cell Motif Enrichment Scores

To map our enrichment results onto our Seurat object (for visualization and to correlate the results with gene expression data), we used the NMF cell scores. NMF cell scores are weights that specify how strongly each metagene contributes to a cell's overall transcriptional profile. We therefore used these cell scores to generate a weighted average of enrichment scores for each motif at a single-cell level, such that the enrichment score from a high scoring metagene contributed more strongly to a cell's enrichment score than a lowly scoring metagene.

(*snippet from 02_enrichmentAnalysis/hydraRegulators.R*)

```R
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
```

We then used the enrichment scores we calculated for our Hydra atlas to generate a heat map of single-cell enrichment scores averaged across different cell types. For this plot, and for other downstream analyses, we dropped any motifs that were not found to show signs of conservation in our cross-species whole-genome alignment (specified in the `motifConservationStats.csv` file generated in `07_genomeConservation.md`). This section of code also makes use of the file `k56_mg_annot.csv`, which was manually generated and includes descriptive names for each metagene. There's also the file `mg_order.csv` which specifies the order in which to metagenes appear in the heat map. 

(*snippet from 02_enrichmentAnalysis/hydraRegulators.R*)

```R
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
```

<img src="resources/motifHeatmap.png" alt="motifHeat" style="zoom:10%;" />

#### Removing Redundant Motifs from the Enrichment Results

This initial heat map had over 300 rows, and wasn't very compact. To reduce the number of rows in the plot, we removed redundant motif entries, which we defined as motifs with highly similar sequence composition and enrichment patterns. We had already grouped motifs according to sequence similarity. In the following section of the analysis, we performed another clustering analysis where we clustered motifs based on correlation scores calculated by comparing motif enrichment results. 

We grouped the motifs in a way that integrated both the sequence-based and the enrichment-based clustering analyses, such that motifs were only considered redundant if they had similar sequence composition **and** enrichment results. To collapse a group of motifs that were flagged as redundant, we presented the averaged enrichment profile of all motifs within that group. 

(*snippet from 02_enrichmentAnalysis/hydraRegulators.R*)

```R
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
```

![motDistClustEn](resources/motDistClustEn.png)

We then generated a new, 'low-redundancy' version of the heat map. 

(*snippet from 02_enrichmentAnalysis/hydraRegulators.R*)

```R
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
```







<img src="resources/motifHeatmapLR.png" alt="motHeatLr" style="zoom:10%;" />

#### Identifying Candidate Regulators of Enriched Motifs

We next linked enriched motifs to putative regulators. To do this we needed to correlate a motif's enrichment pattern in the single-cell atlas to the expression patterns of predicted transcription factors. One potential challenge for performing this type of correlation analysis involving single-cell gene expression data is that the low depth of scRNA-seq makes it susceptible to frequent 'drop-outs', where a gene with moderate expression has zero counts in a non-trivial number of cells. To mitigate the issues this might cause for the correlation analysis between motif enrichment and gene expression, we generated a matrix of imputed read counts based on our NMF analysis. 

The goal of NMF is to generate two matrices (the cell score matrix and the gene score matrix) that, when multiplied together, create an approximation of the original data matrix (the single-cell gene expression matrix). We took advantage of this to create an expression matrix that roughly recapitulated the original expression data, but without any drop-outs, smoothing out the data and making it more suitable for correlation analyses. This was the matrix that we used to determine the correlation score between predicted transcription factors and enriched motifs.

(*snippet from 02_enrichmentAnalysis/hydraRegulators.R*)

```R
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
```

From this motif/transcription factor correlation matrix, we generated a summary table of all putative co-expression regulators. A candidate regulator was defined as a transcription factor that had a correlation score of 0.5 or greater with a motif that it could plausibly bind based on its Pfam domain content. This table included the gene ID of the candidate regulator as well as a list of all the enriched motifs that transcription factor could potentially be regulating. We also added some functional annotation data for the candidate regulator, such as predicted orthologs and Pfam domains.

(*snippet from 02_enrichmentAnalysis/hydraRegulators.R*)

```R
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
```

To generate the motif enrichment and gene expression plots provided in the paper we used the template below. The motif IDs and gene IDs were swapped to create the other plots presented throughout the text.

(*snippet from 02_enrichmentAnalysis/hydraRegulators.R*)

```R
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
```

![ebfM](resources/ebfM.png)

![ebfGexp](resources/ebf.png)

### Analysis Using Shuffled Binding Motifs

We next wanted to determine the extent to which the enrichment patterns we observed above were due to true biological signal, as opposed to being artifacts driven by chance. To test this, we repeated the single-cell motif enrichment analysis using shuffled instead of bona fide transcription factor binding motifs. We had already generated a list of shuffled motifs for the analysis described in `07_genomeConservation.md`. In addition we had subjected predicted instances of these shuffled motifs to the same filtering criteria (i.e., requiring conservation across multiple *Hydra* genomes and localization to open chromatin) we used for generating binding site predictions for our database of bona fide binding sites. 

Because the shuffled motifs should be random, non-functional sequences, any enrichment patterns we identify using the shuffled motifs will likely be artifacts. If we get similar enrichment results using both the shuffled and bona fide versions of a motif, we can determine that our enrichment results are likely not biologically meaningful; however, if the bona fide motif has a dramatically different enrichment pattern than the shuffled motif, then we can conclude that the enrichment results are not purely the result of sequence bias.

#### Performing Gene Set Enrichment Analyses

To repeat our enrichment analysis using shuffled motifs, we used the `conShufMotsATAC_finalhits.txt` table to generate gene sets. We then performed gsea for each metagene in the *Hydra* single cell atlas using these shuffled motif gene sets:

(*snippet from 02_enrichmentAnalysis/hydraRegulatorsShuf.R*)

```R
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
mots <- read.delim('../../alignment_conservation/motifAnnot/conShufMotsATAC_finalhits.txt')

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
saveRDS(enrichDF,'enrichDFnCDSATAC_shuf.rds')
```

#### Generating Single-Cell Motif Enrichment Scores

We then used NMF cell scores to translate the metagene enrichment results into single-cell enrichment scores for each motif:

(*snippet from 02_enrichmentAnalysis/hydraRegulatorsShuf.R*)

```R
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

saveRDS(motScores,'motScores_shuf.rds')
```

We then visualized these motif enrichment scores using a motif by cell type heat map:

(*snippet from 02_enrichmentAnalysis/hydraRegulatorsShuf.R*)

```R
####motif enrichment heatmap####

#load metagene descriptions
mgD <- read.csv('k56_mg_annot.csv',header = F)

mgD$V1 <- paste0('mg',mgD$V1)

#initialize plotting data object
nesPlot <- enrichDF

#rename the columns to have the more descriptive metagene name
colnames(nesPlot) <- mgD$V2

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
pdf('motifHeatmap_shuf.pdf',width = 10,height = 60)
heatmap.2(as.matrix(nesPlot),
          Colv = F,
          scale = 'none',
          dendrogram = 'none',
          col = viridis(30),
          trace='none',
          key = F,
          keysize = 0.1,
          margins = c(10,10),
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="average"))
dev.off()
write.csv(nesPlot,'hydraEnrichmentMatrix_shuf.csv')
```

<img src="resources/motifHeatmap_shuf.png" alt="motifHeatmap_shuf" style="zoom:10%;" />

This heat map revealed a large number of enriched motifs when using the shuffled motif set, although markedly fewer metagenes had enriched motifs when compared to the bona fide motif set. 

### Comparing Enrichment Results for Genuine and Shuffled Motifs

To determine if the enrichment patterns we observed were similar for both the shuffled and bona fide versions of each motif (indicating that the enrichment patterns were driven by sequence bias), we calculated correlation scores for the enrichment results for the bona fide and shuffled motif results. We found that the enrichment results were highly different for the vast majority of motifs, indicating that the enrichment patterns we observed were not primarily being driven by sequence bias artifacts.

(*snippet from 02_enrichmentAnalysis/hydraRegulatorsShuf.R*)

```R
####compare shuffled results to the results from bona fide motifs####

#first just compare raw enrichment results

#import the bona fide enrichment results
enrichDF.bf <- readRDS('enrichDFnCDSATAC.rds')

#need to make sure the rows are identical between the two dfs
enrichDF.bf.comp <- lapply(rownames(enrichDF),function(x){
  if(x %in% rownames(enrichDF.bf)){
    return(enrichDF.bf[x,])
  } else {
    return(rep(0,ncol(enrichDF.bf)))
  }
})

enrichDF.bf.comp <- do.call(rbind,enrichDF.bf.comp)

rownames(enrichDF.bf.comp) <- rownames(enrichDF)

#measure similarity between the shuffled and bona fide motif enrichment results
corRes <- vapply(1:nrow(enrichDF), function(x) cor(as.numeric(enrichDF[x,]),as.numeric(enrichDF.bf.comp[x,])),numeric(1))

corRes <- data.frame(mot = rownames(enrichDF),score=corRes)

#need to fix cases where the cor score couldn't be calculated because of invariance 
#in one of the two DFs
corRes[is.na(corRes$score),'score'] <- 0

summary(corRes$score)
```

```
    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
-0.10616  0.00000  0.00000  0.03438  0.00000  1.00000 
```

```R
ggplot(corRes,aes(x='motifCor',y=score)) + geom_violin() + geom_jitter(width = 0.3, height = 0.02) + theme_bw()
```

![enCorScor](resources/enCorScor.png)

One of the criteria we used for identifying transcriptional regulators was that their transcriptional expression pattern should overlap with the enrichment pattern of their target motif. We therefore sought to determine if the enrichment patterns we calculated using bona fide motifs tended to correlate more closely with the expression patterns of their candidate regulators than the enrichment patterns calculated using shuffled motifs. 

To do this we calculated a correlation score for each possible gene/motif pair (determined above based on pfam domains) by comparing the NMF-imputed gene expression values for the gene and the corresponding motif enrichment scores. We calculated this correlation score for both the bona fide and shuffled versions of the motif. We then compared the two correlation scores (one for bona fide and one for shuffled) using a scatter plot, including only on those gene/motif pairs that had a high correlation score (above 0.4) for either motif version.

(*snippet from 02_enrichmentAnalysis/hydraRegulatorsShuf.R*)

```R
#next check to see if the protein/motif correspondence is better with the non-shuffled motif than it is with the
#shuffled motif

#load imputed gene expression matrix
gscores <- readRDS('nmfNormExp.rds')

#get list of genes linked to motifs by pfam domain
pfamMot <- readRDS('../../alignment_conservation/motifDB/motifBindPfam.rds')

#check similarity between motif enrichment and gene expression
corRes <- apply(gscores.tf,2,function(x) apply(motScores, 2, function(y) cor(x,y)))

saveRDS(corRes,'tfMotifCor_shuf.rds')

corRes <- readRDS('tfMotifCor_shuf.rds')

#load gene-motif correlation scores for bona fide motifs
corRes.bf <- readRDS('tfMotifCor.rds')

#generate table of possible motif/gene pairs
corComp <- data.frame(motID = unlist(sapply(1:length(pfamMot), function(x) rep(names(pfamMot)[x],length(pfamMot[[x]])))),
                      geneID = unlist(pfamMot))

#pull correlation score for each gene/motif pair for the shuffled motif data
corComp$cor <- apply(corComp,1,function(x){
  if(x[1] %in% rownames(corRes) & x[2] %in% colnames(corRes)){
    return(corRes[x[1],x[2]])
  } else {
    return(0)
  }
})

#pull correlation score for each gene/motif pair for the bona fide motif data
corComp$corBF <- apply(corComp,1,function(x){
  if(x[1] %in% rownames(corRes.bf) & x[2] %in% colnames(corRes.bf)){
    return(corRes.bf[x[1],x[2]])
  } else {
    return(0)
  }
})

#drop gene/motif pairs that had zeroes for both bona fide and
#shuffled motifs
corComp <- corComp[corComp$cor != 0 | corComp$corBF != 0,]

#for plot, drop gene/motif pairs that did not high scores for either 
#shuffled or bona fide motif sequences
corComp.p <- corComp[corComp$cor > 0.4 | corComp$corBF > 0.4,]

#visualize correspondence between shuffled and bona fide motifs
ggplot(corComp.p,aes(x=cor,y=corBF)) + geom_point() + theme_bw()
```

![motEnCorComp](resources/motEnCorComp.png)

(Note: the x-axis corresponds to scores calculated using shuffled motif enrichment patterns and the y-axis corresponds to scores calculated using bona fide motif enrichment patterns)

The plot above shows that for any given gene/motif pair, the degree of correspondence between enrichment and expression for the shuffled motif has virtually no correlation to the degree of correspondence for the bona fide motif. This can also be demonstrated more concisely by just calculating the correlation between the two scores:

```R
cor(corComp$cor,corComp$corBF)
```

```
[1] 0.09411901
```

In addition, the plot reveals a trend where unshuffled motifs tend to have higher correlation scores than the shuffled motifs, as there are more points with high values on the y-axis than the x-axis. This can also be demonstrated using a box plot showing the distribution of correlation scores for bona fide and shuffled motifs:

(*snippet from 02_enrichmentAnalysis/hydraRegulatorsShuf.R*)

```R
#is the correspondence between motif enrichment and gene expression
#stronger for bona fide motifs than for shuffled motifs?
corComp.vln <- data.frame(corType=rep(c('shuf','bf'),c(nrow(corComp),nrow(corComp))),
                          corScore=c(corComp$cor,corComp$corBF))

ggplot(corComp.vln,aes(x=corType,y=corScore)) + geom_boxplot() + theme_bw()
```

![geneMotCorComp](resources/geneMotCorComp.png)

Although most gene/motif pairs showed no correlation (even for bona fide motif sequences) the overall correspondence between motif enrichment and gene expression was significantly higher for bona fide motifs than for shuffled motifs. 

(*snippet from 02_enrichmentAnalysis/hydraRegulatorsShuf.R*)

```R
#is the difference significant?
t.test(x=corComp$cor,y=corComp$corBF)
```

```
	Welch Two Sample t-test

data:  corComp$cor and corComp$corBF
t = -21.551, df = 36786, p-value < 2.2e-16
alternative hypothesis: true difference in means is not equal to 0
95 percent confidence interval:
 -0.03683734 -0.03069529
sample estimates:
   mean of x    mean of y 
-0.008893753  0.024872563 
```

Overall, these results strongly suggest that the enrichment patterns we observed using our bona fide motif set are not purely driven by sequence bias artifacts. In addition, the enrichment patterns of unshuffled motifs more closely matched the expression patterns of their putative regulators than would be expected by chance.



## Files Associated with This Document



```
10_hydraRegulators/
├── 01_motifPrep
│   ├── clusterMots.R
				R script that uses a similarity matrix generated by homerMotCompare.sh
				to perform hierarchical clustering of the motifs in pooledJasparNR.txt
				based on sequence composition. 
│   ├── conMotATACAnnot.json
				JSON config file used by the UROPA pipeline to link predicted 
				transcription factor binding sites to their target genes based on 
				proximity.
│   ├── conMotATAC_finalhits.txt
				Output generated by UROPA linking predicted transcription factor
				binding sites to their target genes based on proximity.
│   ├── conShufMotATACAnnot.json
				JSON config file used by the UROPA pipeline to link conserved instances 
				of shuffled transcription factor binding motifs to their target genes 
				based on proximity.
│   ├── conShufMotsATAC_finalhits.txt
				Output generated by UROPA linking conserved instances of shuffled 
				transcription factor binding motifs to their target genes based on 
				proximity.
│   ├── fullTF.csv
				Full interproscan protein domain annotations for all predicted 
				transcription factors in the AEP assembly gene models.
│   ├── homerMotCompare.sh
				Shell script that uses the compareMotifs.pl script from HOMER to 
				generate a similarity matrix comparing the sequence composition
				of all motifs in pooledJasparNR.txt
│   ├── jasparIDs.txt
				Text file that lists all JASPAR motif IDs in pooledJasparNR.txt. Used
				by protDBLink.R.
│   ├── motComp.txt
				Text file containing a similarity matrix comparing the sequence 
				composition of all motifs in pooledJasparNR.txt. Generated by 
				homerMotCompare.sh
│   ├── motifBindPfam.rds
				Binary R object file containing a list object that links individual
				JASPAR motifs to AEP gene models that could potential bind those 
				motifs based on pfam protein domain annotations.
│   ├── motif_clusters.csv
				Table containing the clustering results from clusterMots.R. Motifs
        were grouped using a hierarchical clustering analysis based on
        sequence similarities.
│   ├── motifInfo.csv
				Table that provides the name and transcription factor family associated
				with each motif in pooledJasparNR.txt.
│   ├── pooledJasparNR.homer.motifs
				HOMER-formatted version of the pooledJasparNR.txt motif database.
│   ├── protDBLink.R
				R script that links AEP gene models to putative target binding motifs
				based on pfam domain annotations. Specifically, it identifies the
				DNA binding domain associated with each JASPAR motif and then links
				that motif to all AEP gene models that also contain that binding
				domain.
│   ├── tfIDs.txt
				Text file listing the gene ID of all putative transcription factors
				in the AEP gene models.
│   └── tfList.R
				R script that uses inteproscan protein domain annotations to identify
				putative transcription factors in the AEP gene models. It does this
				using a combination of GO terms and lists of known DNA binding domains.
└── 02_enrichmentAnalysis
    ├── enrichDFnCDSATAC.rds
    		Binary R object file containing a motif by metagene matrix of normalized
    		enrichment scores (NES) for the Hydra single cell atlas. Higher values
    		indicate stronger enrichment. NES values that did not pass our significance
    		theshold (P. adj < 0.01) were set to zero.
    ├── enrichDFnCDSATAC_shuf.rds
    		Binary R object file containing a motif by metagene matrix of normalized
    		enrichment scores (NES) for the Hydra single cell atlas calculated using 
    		shuffled control versions of each motif. Higher values indicate stronger 
    		enrichment. NES values that did not pass our significance theshold (P. 
    		adj < 0.01) were set to zero.
    ├── finalMotClust.csv
    		Table containing the final motif clustering results that grouped JASPAR
    		motifs based on both sequence composition (results taken from 
    		motif_clusters.csv) and motif enrichment patterns. This was used to reduce
    		the redundancy of the final enrichment results matrix 
    		(hydraEnrichmentMatrixLR.csv)
    ├── hydraEnrichmentMatrix.csv
    		Motif by cell type matrix of enrichment scores in the Hydra single cell
    		atlas. These values were used to generate the heatmaps presented in the
    		paper. 
    ├── hydraEnrichmentMatrixLR.csv
    		Low redundancy version of hydraEnrichmentMatrix.csv where rows (i.e., motifs)
    		have been collapsed according to the clustering in finalMotClust.csv.
    ├── hydraEnrichmentMatrix_shuf.csv
    		Motif by cell type matrix of enrichment scores in the Hydra single cell
    		atlas calculated using shuffled control versions of each motif.
    ├── hydraRegulators.R
    		R script that uses gene set enrichment analysis to identify motifs enriched
    		in each metagene in the Hydra single cell atlas. It then maps these values
    		onto individual cells to generate single-cell enrichment scores. Finally,
    		it compares motif enrichment patterns to the expression patterns of putative
    		transcription factors that could plausibly bind each motif to identify
    		condidate regulators of cell-type-specific expression.
    ├── hydraRegulatorsShuf.R
    		R script that repeats the enrichment analysis from hydraRegulators.R, but 
    		instead uses shuffled control versions of each motif. It then compares
    		the enrichment results for the shuffled and non-shuffled versions to 
    		evaluate the degree to which the bona fide motif enrichment results
    		are meaningfully different from the results generated using random noise
    		(simulated using the shuffled motifs).
    ├── k56_mg_annot.csv
    		Table that provides descriptive names for each metagene in the Hydra single
    		cell atlas.
    ├── mg_order.csv
    		Table that specifies the column order of the motif enrichment heatmaps 
    		(e.g., hydraEnrichmentMatrixLR.csv)
    ├── motScores.rds
    		Binary R object file containing a motif by cell matrix of motif enrichment 
    		scores for the Hydra single cell atlas.
    ├── motScores_shuf.rds
    		Binary R object file containing a motif by cell matrix of motif enrichment 
    		scores for the Hydra single cell atlas calculated using shuffled control 
    		versions of each motif.
    ├── nmfNormExp.rds
    		Binary R object file containing a gene by cell matrix of imputed/smoothened
    		expression values for the Hydra single cell atlas calculated using NMF. Used
    		to facilitate single cell gene expression correlation analyses.
    ├── pfamMotCor.csv
    		Table containing correlation scores for every motif/gene pair specified in 
    		motifBindPfam.rds. The score reflects how well the enrichment pattern for the
    		motif matched the expression pattern of the putative regulator.
    ├── regulatorResults.csv
    		Table listing all candidate regulators of cell type specific transcription
    		in the Hydra single cell atlas. All putative transcription factors whose 
    		expression pattern was correlated (correlation score threshold of ≥ 0.5)
    		with the enrichment pattern of at least one motif were included in the table.
    		The table also lists all motifs that passed the correlation threshold for each
    		transcription factor, as well as some functional annotation information for each 
    		candidate regulator.
    ├── tfMotifCor.rds
    		Motif-by-gene correlation matrix comparing single-cell GSEA enrichment scores to 
    		the NMF-imputed expression patterns of putative transcription factors in the 
    		Hydra atlas.
		└── tfMotifCor_shuf.rds
				Motif-by-gene correlation matrix comparing single-cell GSEA enrichment scores of 
				shuffled control motifs to the NMF-imputed expression patterns of putative 
				transcription factors in the Hydra atlas.
				
```



