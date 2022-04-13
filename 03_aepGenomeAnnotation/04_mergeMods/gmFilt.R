library(rstudioapi)
library(plyr)

#set the working directory to be the folder in which this script is located
setwd(dirname(getActiveDocumentContext()$path))

bE.in <- read.delim("brakerExoOlap.bed", header = F)

#name reformatting
bE.in$V10 <- gsub(';.*','',bE.in$V10)
bE.in$V10 <- gsub('ID=','',bE.in$V10)
bE.in$V4 <- gsub('ID=','',bE.in$V4)
bE.in$V4 <- gsub(';.*','',bE.in$V4)
bE.in$V4 <- gsub('[.]t.*','',bE.in$V4)

#drop identical overlap pairs (caused by isoforms)
bE.in <- bE.in[!duplicated(paste(bE.in$V4,bE.in$V10)),]

#self olap for exonerate 
e.in <- read.delim("exoExoOlap.bed", header = F)

e.in$V10 <- gsub(';.*','',e.in$V10)
e.in$V10 <- gsub('ID=','',e.in$V10)

e.in$V4 <- gsub(';.*','',e.in$V4)
e.in$V4 <- gsub('ID=','',e.in$V4)

e.in <- e.in[e.in$V4 != e.in$V10,]

blastHits <- read.delim("blastpPrimProt.outfmt6", header = F)
blastHits$V1 <- gsub("-mRNA.*","",blastHits$V1)

#get sum of bit score for all blast hits for a given gene model
blastHits.score <- aggregate(blastHits[,12], list(blastHits$V1), FUN = sum)

#only keep the isoform with the highest score
blastHits.score$Group.1 <- gsub('[.]t\\d*','',blastHits.score$Group.1)
blastHits.score <- blastHits.score[order(-blastHits.score$x),]
blastHits.score <- blastHits.score[!duplicated(blastHits.score$Group.1),]

#add braker gene model blast score to overlap table
bE.in$bScore <- mapvalues(bE.in$V4,from = blastHits.score$Group.1, to = blastHits.score$x, warn_missing = F)

bE.in[grepl('^file',bE.in$bScore),'bScore'] <- 0

bE.in$bScore <- as.numeric(bE.in$bScore)

#add exonerate gene model blast score to overlap table
bE.in$eScore <- mapvalues(bE.in$V10,from = blastHits.score$Group.1, to = blastHits.score$x, warn_missing = F)

bE.in[!grepl('^\\d',bE.in$eScore),'eScore'] <- 0

bE.in$eScore <- as.numeric(bE.in$eScore)

#get the id of the gene model that has the higher score
bE.in$better <- bE.in$bScore >= bE.in$eScore
bE.in$worse <- ''

bE.in[bE.in$better == FALSE,'worse'] <- bE.in[bE.in$better == FALSE,4]
bE.in$better[bE.in$better == FALSE] <- bE.in[bE.in$better == FALSE,10]
bE.in[bE.in$better == TRUE,'worse'] <- bE.in[bE.in$better == TRUE,10]
bE.in$better[bE.in$better == TRUE] <- bE.in[bE.in$better == TRUE,4]

bE.in <- bE.in[!(bE.in$better %in% bE.in$worse),]

e.in$e1Score <- mapvalues(e.in$V4,from = blastHits.score$Group.1, to = blastHits.score$x, warn_missing = F)
e.in[!grepl('^\\d',e.in$e1Score),'e1Score'] <- 0
e.in$e1Score <- as.numeric(e.in$e1Score)

e.in$e2Score <- mapvalues(e.in$V10,from = blastHits.score$Group.1, to = blastHits.score$x, warn_missing = F)
e.in[!grepl('^\\d',e.in$e2Score),'e2Score'] <- 0
e.in$e2Score <- as.numeric(e.in$e2Score)

#get the id of the gene model that has the higher score
e.in$better <- e.in$e1Score >= e.in$e2Score
e.in$worse <- ''

e.in[e.in$better == FALSE,'worse'] <- e.in[e.in$better == FALSE,4]
e.in$better[e.in$better == FALSE] <- e.in[e.in$better == FALSE,10]
e.in[e.in$better == TRUE,'worse'] <- e.in[e.in$better == TRUE,10]
e.in$better[e.in$better == TRUE] <- e.in[e.in$better == TRUE,4]

e.in <- e.in[!(e.in$e1Score == 0 & e.in$e2Score == 0),]

e.in$sPair <- apply(e.in[,c(4,10)],1, function(x) {
  vIn <- as.vector(x)
  vIn <- vIn[order(vIn)]
  vOut <- paste(vIn,collapse = '_')
  return(vOut)
})

e.in <- e.in[!duplicated(e.in$sPair),]

e.exclude <- unique(e.in$worse)

#pull in all gene ID names
exoGMs <- read.delim("exoGenes.bed", header = F)
bGMs <- read.delim("brakerGenes.bed", header = F)

#name reformating
exoGMs$V4 <- gsub(';.*','',exoGMs$V4)
exoGMs$V4 <- gsub('ID=','',exoGMs$V4)

bGMs$V4 <- gsub(';.*','',bGMs$V4)
bGMs$V4 <- gsub('ID=','',bGMs$V4)

#keep any gene models that had a blast hit and didn't have an overlap
bGMs <- bGMs[!(bGMs$V4 %in% bE.in$V4),]
bGMs <- bGMs[bGMs$V4 %in% blastHits.score$Group.1,]

exoGMs <- exoGMs[!(exoGMs$V4 %in% bE.in$V10),]
exoGMs <- exoGMs[exoGMs$V4 %in% blastHits.score$Group.1,]

bE.in <- bE.in[!(bE.in$bScore == 0 & bE.in$eScore == 0),]

#get list of gene models to keep
gm.keep <- c(bGMs$V4,exoGMs$V4,bE.in$better)
gm.keep <- unique(gm.keep[!(gm.keep %in% e.exclude)])

#filter braker gene models
bGff <- read.delim("braker.fix.gff3",header = F)

bGff$gID <- gsub('.*gene_id=','',bGff$V9)
bGff$gID <- gsub(';','',bGff$gID)

bGff <- bGff[bGff$gID %in% gm.keep,]

bGff$gID <- NULL

write.table(bGff, file = 'braker.keep.gff3', row.names = F, col.names = F, sep = '\t', quote = F)

#filter exonerate gene models
eGff <- read.delim('exoCat.complete.gff3', header = F, skip = 1)

eGff$gID <- gsub('.*gene_id=','',eGff$V9)
eGff$gID <- gsub(';.*','',eGff$gID)

eGff <- eGff[eGff$gID %in% gm.keep,]

eGff$gID <- NULL

write.table(eGff, file = 'exo.keep.gff3', row.names = F, col.names = F, sep = '\t', quote = F)


