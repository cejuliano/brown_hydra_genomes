library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

inGff <- read.delim('clFinal.merge.longestIso.rfmt.tag.gtf',sep = '\t', header = F, skip = 1)

#drop non-coding gene entries
inGff <- inGff[inGff$V3 != 'region',]
inGff <- inGff[inGff$V3 != 'tRNA',]
inGff <- inGff[inGff$V3 != 'ncRNA_gene',]

#give the gene ID a simpler format
inGff$gID <- gsub('.*gene_id ([^;]*).*','\\1',inGff$V9)

inGff$gID <- gsub('GENE.(.*)~~.*','\\1',inGff$gID)
inGff$gID <- gsub('gene:','',inGff$gID)

inGff$gID <- gsub(' +$','',inGff$gID)

inGff <- inGff[!grepl('rank \\d+$',inGff$gID),]

inGff <- inGff[inGff$V3 %in% c('gene','mRNA','exon','CDS','five_prime_UTR','three_prime_UTR'),]


inGff$gID <- paste0('gene_id "',inGff$gID,'" ; gene_name "',inGff$gID,'" ; transcript_id "',inGff$gID,'" ; transcript_name "',inGff$gID,'"')

write.table(inGff[,c(1:8,10)],file = 'clFinal.merge.longestIso.rename.gtf',sep = '\t',quote = F, row.names = F, col.names = F)
