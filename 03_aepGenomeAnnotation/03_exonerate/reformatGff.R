args <- commandArgs(trailingOnly = T)

inGff <- read.delim(file = args[1], sep = "\t", stringsAsFactors = F, header = F)

inGff <- inGff[!grepl("^#",inGff$V1),]
inGff <- inGff[!grepl("^-",inGff$V1),]

stCoord <- gsub(".*:","",inGff[1,1])
stCoord <- as.numeric(gsub("-.*","",stCoord))

inGff$V4 <- inGff$V4 + stCoord
inGff$V5 <- inGff$V5 + stCoord

inGff$V1 <- gsub(":.*","",inGff$V1)

gName <- inGff[grepl("gene",inGff$V3),9]
gName <- gsub(".*sequence ","",gName)
gName <- gsub(" ;.*","",gName)

inGff$V9 <- gsub("gene_id 0 ; ","",inGff$V9)

gName <- paste0(" ; gene_id ",gName)


inGff$V9 <- paste0(inGff$V9,gName)
inGff$V9 <- gsub("^ ; ","",inGff$V9)

end5 <- inGff[inGff$V3 == 'utr5',5]
end5 <- max(end5)

if(inGff[1,7] == '+') {
  inGff[inGff$V3 == 'cds' & inGff$V4 <= end5,4] <- end5 + 1
}

end5 <- inGff[inGff$V3 == 'utr5',4]
end5 <- min(end5)

if(inGff[1,7] == '-') {
  inGff[inGff$V3 == 'cds' & inGff$V5 >= end5,5] <- end5 - 1
}

write.table(file = paste0(args[1],".gff"), inGff, sep = "\t", quote = F, row.names = F, col.names = F)