options(stringsAsFactors = F)
args <- commandArgs(trailingOnly = T)

inGff <- read.delim(args[1], sep = "\t", skip = 1, header = F)

gName <- inGff[inGff$V3 == "gene"[1],9]

gName <- gsub("ID=","",gName)
gName <- gsub(";.*","",gName)

inGff$V3 <- gsub("^RNA$","mRNA",inGff$V3)
badRNA <- which(inGff$V3 == "mRNA")

for(i in badRNA) {
  badFormat <- gsub('.*Parent=([^;]*);.*',"\\1",inGff[i,9])
  goodFormat <- toupper(badFormat)
  goodFormat <- paste0("Parent=",goodFormat)
  inGff[i,9] <- gsub(paste0("Parent=",badFormat),goodFormat,inGff[i,9])
}

inGff$V9 <- gsub("nbisL\\d-cds",paste0(gName,"-mRNA"),inGff$V9)
inGff$V9 <- gsub("nbisL\\d",gName,inGff$V9)
inGff$V9 <- gsub("nbis",gName,inGff$V9)
inGff$V9 <- gsub("ID=exon",paste0("ID=",gName,"-exon"),inGff$V9)
inGff$V9 <- gsub("ID=cds",paste0("ID=",gName,"-cds"),inGff$V9)
inGff$V9 <- gsub("ID=(\\d+);",paste0("ID=",gName,"-intron","-\\1;"),inGff$V9)
inGff$V9 <- gsub('T(\\d+)AEP','t\\1aep',inGff$V9)

inGff[inGff$V3 == 'mRNA',9] <- gsub('^(ID=[^;]+)exon','\\1mRNA',inGff[inGff$V3 == 'mRNA',9])
inGff[inGff$V3 %in% c('intron','exon','cds','three_prime_UTR','five_prime_UTR'),9] <- gsub('(Parent=[^-;]+-)[^-;]+','\\1mRNA',inGff[inGff$V3 %in% c('intron','exon','cds','three_prime_UTR','five_prime_UTR'),9])



write.table(inGff, file = gsub(".gff3",".pfix.gff3",args[1]), quote = F, sep = "\t", row.names = F, col.names = F)  
