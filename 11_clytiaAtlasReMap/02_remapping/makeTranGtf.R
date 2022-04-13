library(rstudioapi)
library(Biostrings)
library(plyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tranSeq <- readDNAStringSet('clFinal.longestIso.tran.fa')

tranSeq.dim <- data.frame(ID = names(tranSeq), length = nchar(tranSeq))

tranSeq.dim$ID <- gsub(' .*','',tranSeq.dim$ID)

inGtf <- read.delim('clFinal.merge.longestIso.rename.gtf',header = F, skip = 1)

inGtf <- inGtf[inGtf$V3 == 'mRNA',]

inGtf$V1 <- gsub('.*gene_id ([^ ]+).*','\\1',inGtf$V9)

inGtf$V5 <- mapvalues(inGtf$V1, from = tranSeq.dim$ID, to = tranSeq.dim$length)

inGtf$V4 <- 1

inGtf$V7 <- '+'

inGtf$V3 <- 'exon'

inGtf$new9 <- inGtf$V9


inGtf$new9 <- gsub('([^;]) ([^;])','\\1 "\\2',inGtf$new9)
inGtf$new9 <- gsub(' ;','" ;',inGtf$new9)
inGtf$new9 <- gsub(' $','"',inGtf$new9)

inGtf$V9 <- inGtf$new9

write.table(inGtf[,1:9], file = 'clFinal.transcriptome.gtf', row.names = F, col.names = F, quote = F, sep = '\t')
