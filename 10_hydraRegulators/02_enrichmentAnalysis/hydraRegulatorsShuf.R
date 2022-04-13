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
enrichDF <- readRDS('enrichDFnCDSATAC_shuf.rds')

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
motScores <- readRDS('motScores_shuf.rds')

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

ggplot(corRes,aes(x='motifCor',y=score)) + geom_violin() + geom_jitter(width = 0.3, height = 0.02) + theme_bw()

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

cor(corComp$cor,corComp$corBF)

#is the correspondence between motif enrichment and gene expression
#stronger for bona fide motifs than for shuffled motifs?
corComp.vln <- data.frame(corType=rep(c('shuf','bf'),c(nrow(corComp),nrow(corComp))),
                          corScore=c(corComp$cor,corComp$corBF))

ggplot(corComp.vln,aes(x=corType,y=corScore)) + geom_boxplot() + theme_bw()

#is the difference significant?
t.test(x=corComp$cor,y=corComp$corBF)
