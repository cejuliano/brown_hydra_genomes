setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

inG <- read.delim("HVAEP1.geneModels.pUpdate1.filt.gff3", header = F, skip = 1)

inG <- inG[!grepl("^#",inG$V1),]

#PASA added some of its own models, which we don't really want
#so we just drop them
inG <- inG[!grepl('novel_model|temp_model',inG[,9]),]

#get the gene IDs for each row
inG$gID <- gsub(".*HVAEP1_[TG](\\d+).*","\\1",inG$V9)

#partition rows by gene ID
inG.List <- split(inG, inG$gID)

#sort rows by coords, from 5' to 3' relative to the gene in question
inG.List <- lapply(inG.List, function(x) if(x$V7[1] == '-'){
  return(x[order(-x$V4),])} else {
    return(x[order(x$V4),])
  })

#initialize list of problematic short 3' UTRs
short3 <- list()

#for loop to check all gene models for problematic UTRs
for(i in 1:length(inG.List)){
  subG <- inG.List[[i]]
  #if the gene doesn't even have 3' UTR we can skip it
  if(nrow(subG[subG$V3 == 'three_prime_UTR',]) == 0) {
    next
  }
  #pull 3' UTR rows and calculate each UTR segments length
  prime3 <- subG[subG$V3 == 'three_prime_UTR',]
  prime3L <- prime3[nrow(prime3),5] - prime3[nrow(prime3),4]
  
  #next look for exons at the end of genes that are exclusively made of UTR sequence
  exTest <- subG[subG$V3 == 'exon',]
  exTest <- (exTest[nrow(exTest),4] == prime3[nrow(prime3),4]) & (exTest[nrow(exTest),5] == prime3[nrow(prime3),5])
  
  #flag any 3' UTR-only exons if they are shorter than 20 bp
  if(prime3L <= 20 & exTest){
    short3[[as.character(subG[1,'gID'])]] <- subG
  }
}

#initialize a list of fixed UTRs
short3.fix <- list()

# go through the list of problematic 3' UTRs and drop them from the gene model
for(i in 1:length(short3)) {
  subG <- short3[[i]]
  if(subG$V7[1] == '-') {
    #get current boundary from faulty utr
    oldB <- subG[nrow(subG),4]
    
    #get the new boundary from the next leftmost thing
    newB <- unique(subG$V4)
    newB <- newB[length(newB) - 1]
    
    #delete the bad UTR
    subG <- subG[-which(subG$V3 == 'three_prime_UTR' & subG$V4 == oldB),]
    subG <- subG[-which(subG$V3 == 'exon' & subG$V4 == oldB),]
    
    #update the new Boundary for other rows
    subG[subG$V4 == oldB,'V4'] <- newB
    
    short3.fix[[as.character(subG$gID[1])]] <- subG
    
  } else {
    #get current boundary from faulty utr
    oldB <- subG[nrow(subG),5]
    
    #get the new boundary from the next leftmost thing
    newB <- unique(subG$V5)
    newB <- newB[length(newB) - 1]
    
    #delete the short UTR and it's exon
    subG <- subG[-which(subG$V3 == 'three_prime_UTR' & subG$V5 == oldB),]
    subG <- subG[-which(subG$V3 == 'exon' & subG$V5 == oldB),]
    
    #update the new Boundary for other rows
    subG[subG$V5 == oldB,5] <- newB
    
    short3.fix[[as.character(subG$gID[1])]] <- subG
  }
}

#merged fixed gene models with remaining entries
outG <- inG.List[!(names(inG.List) %in% names(short3.fix))]

outG <- c(outG, short3.fix)

outG <- do.call(rbind,outG)

#make sure output is coordinate sorted
chrNum <- as.numeric(gsub('chr-','',outG$V1))

outG <- outG[order(chrNum,outG$V4),]

#give every row the appropriate gene ID tag
outG$gID <- paste0(';gene_id=HVAEP1_G',outG$gID)

outG$V9 <- paste0(outG$V9,outG$gID)

outG$V9 <- gsub('[.]\\d[.](\\d)[.]','_\\1_',outG$V9)

outG$V9 <- gsub('Name=[^;]+;','',outG$V9)

write.table(outG[,1:9], file = "HVAEP1.geneModels.pUpdate1.filt.uFix.gff3", quote = F, row.names = F, col.names = F, sep = '\t')
