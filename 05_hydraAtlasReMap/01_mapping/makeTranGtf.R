library(rstudioapi)
library(Biostrings)
library(plyr)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

tranSeq <- readDNAStringSet('HVAEP1.tran.longestIso.fa')

#get the lengths of all sequences in the fasta file
tranSeq.dim <- data.frame(ID = names(tranSeq), length = nchar(tranSeq))

#get transcript names
tranSeq.dim$ID <- gsub(' .*','',tranSeq.dim$ID)

inGtf <- read.delim('HVAEP1.GeneModels.longestIso.gtf',header = F, skip = 1)

#pull only the transcript rows from the gtf
inGtf <- inGtf[inGtf$V3 == 'mRNA',]

#set the chrom column to be the transcript ID
inGtf$V1 <- gsub('.*(HVAEP1_T\\d+[.]\\d+).*','\\1',inGtf$V9)

#use the length info to set the end coordinate for each transcript
inGtf$V5 <- mapvalues(inGtf$V1, from = tranSeq.dim$ID, to = tranSeq.dim$length)

#set all start coordinates to one
inGtf$V4 <- 1

#set all strand info to +
inGtf$V7 <- '+'

#make all features exons
inGtf$V3 <- 'exon'

inGtf$new9 <- inGtf$V9

#some tag reformatting to fix the specific requirements of the Drop-seq pipeline
inGtf$new9 <- gsub('ID[^;]+;','',inGtf$new9)
inGtf$new9 <- gsub('Parent[^;]+; ','',inGtf$new9)

inGtf$new9 <- gsub('HVAEP','"HVAEP',inGtf$new9)
inGtf$new9 <- gsub(' ;','";',inGtf$new9)
inGtf$new9 <- gsub(' $','"',inGtf$new9)

inGtf$new9 <- paste0(inGtf$new9,';',gsub('_id','_name',inGtf$new9))

inGtf$V9 <- inGtf$new9

write.table(inGtf[,1:9], file = 'HVAEP1.transcriptome.gtf', row.names = F, col.names = F, quote = F, sep = '\t')