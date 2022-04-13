library(rstudioapi)
library(plyr)
library(gplots)
library(viridis)

setwd(dirname(getActiveDocumentContext()$path))

#list clytia enrichment result files
clEn.files <- list.files('enOutShuf', pattern='ame.tsv',recursive = T, full.names = T)

#import clytia enrichment tables
clEn <- lapply(clEn.files, read.delim)

#drop comment lines from results
clEn <- lapply(clEn,function(x) x[!grepl('^#',x[,1]),])

#attach metagene IDs to objects in list
names(clEn) <- gsub('.*(mg\\d+)/ame.tsv','\\1',clEn.files)

#drop empty results
clEn <- clEn[sapply(clEn,length) > 0]

#add in column listing the metagene used to 
#calculate the states for that row
clEn.DF <- lapply(1:length(clEn),function(x){
  newDf <- clEn[[x]]
  newDf$mg <- names(clEn)[x]
  return(newDf)
})

#collapse into DF
clEn.DF <- do.call(rbind,clEn.DF)

####plotting####

#generate clytia enrichment scores
#will just base it on fold enrichment
clEn.DF$fc <- clEn.DF$X.TP/clEn.DF$X.FP

write.csv(clEn.DF,'clEnrichmentResShuf.csv',row.names = F)

#get list of all enriched motifs
motsUse <- unique(clEn.DF$motif_ID)

#load cl metagene cell scores
clNMF <- read.delim('../nmf/cl_fine.usages.k_37.dt_0_13.consensus.txt',row.names = 1)

#normalize cell scores so that the scores for a single cell sum to 1
clNMF <- t(apply(clNMF,1,function(x) x/sum(x)))

#for each enriched motif, go through each metagene that shows enrichment and multiply the fold enrichment by the cell score
enScores <- lapply(motsUse,function(x){
  motRes <- clEn.DF[clEn.DF$motif_ID == x,]
  motScores <- apply(motRes,1,function(y){
    mgUse <- as.numeric(gsub('mg','',y[18]))
    return(clNMF[,mgUse] * as.numeric(y[19]))
  })
  motScores <- apply(motScores,1,sum)
  return(motScores)
})

#get motif name (as opposed to JASPAR ID)
names(enScores) <- mapvalues(motsUse,from=clEn.DF$motif_ID,to=clEn.DF$motif_alt_ID,warn_missing = F)

#collapse to matrix
#this can be used to generate UMAP plots of enrichment scores
enScores <- do.call(cbind,enScores)

saveRDS(enScores,'clEnScoresShuf.rds')

#make matrix for heatmap plotting (motif by metagene enrichment score heatmap)
enScores.hm <- lapply(motsUse,function(x){
  motRes <- clEn.DF[clEn.DF$motif_ID == x,]
  mgIndex <- as.numeric(gsub('mg','',motRes$mg))
  enRow <- rep(0,37)
  enRow[mgIndex] <- motRes$fc
  enRow <- enRow/max(enRow)
  names(enRow) <- paste0('mg',1:37)
  return(enRow)
})

#collapse to matrix
enScores.hm <- do.call(rbind,enScores.hm)

#add in motif name to rownames
motsUse.names <- unique(clEn.DF$motif_alt_ID)

rownames(enScores.hm) <- paste(motsUse,motsUse.names,sep=' ')

#bring in metagene descriptions for colnames
mgAnnot <- read.csv('../nmf/clMetaAnnot.csv',header=F)

colnames(enScores.hm) <- paste(colnames(enScores.hm),mgAnnot$V2,sep=' ')

#reorder the columns so their grouped more logically
colOrder <- read.csv('../nmf/clMgOrder.csv',header=F)

enScores.hm <- enScores.hm[,colOrder$V1]

#plot enrichment heatmap
pdf('clMotifHeatmapShuf.pdf',width = 15,height = 60)
heatmap.2(enScores.hm,
          Colv = F,
          scale = 'none',
          dendrogram = 'none',
          col = viridis(30),
          trace='none',
          key = F,
          keysize = 0.1,
          margins = c(10,10),
          colsep = c(3,9,15,24,33,35),
          distfun = function(x) as.dist(1-cor(t(x))),
          hclustfun = function(x) hclust(x, method="average"))
dev.off()

write.csv(enScores.hm,file='clEnHeatmapShuf.csv')
