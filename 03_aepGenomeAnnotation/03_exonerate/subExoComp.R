library(rstudioapi)

#set the working directory to be the folder in which this script is located
setwd(dirname(getActiveDocumentContext()$path))

e.full <- read.delim("exoCat.gff3",header = F)

e.keep <- read.delim("exoHeaders.txt", header = F)
e.keep$V1 <- gsub('>','',e.keep$V1)
e.keep$V1 <- gsub('.*gene=','',e.keep$V1)
e.keep$V1 <- gsub(' .*','',e.keep$V1)

e.full$gID <- e.full$V9
e.full$gID <- gsub('.*gene_id=','',e.full$gID)
e.full$gID <- gsub(';.*','',e.full$gID)

e.sub <- e.full[e.full$gID %in% e.keep$V1,]
e.sub$gID <- NULL

write.table(e.sub,file = 'exoCat.complete.gff3',row.names = F, col.names = F, sep = '\t',quote = F)