library(rstudioapi)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

inGff <- read.delim('clPasa.merge.longestIso.gff3',sep = '\t', header = F, skip = 1)

#pull out the root assembly ID (and isoform number) for each gene model
inGff$aID <- gsub('.*(asmbl_\\d+[.]p\\d+).*','\\1',inGff$V9)

#make a new object, where we'll find the problematic overlapping genes
idProb <- inGff[inGff$V3 == 'gene',]

#drop the p# isoform suffix
idProb$aIDsub <- gsub('[.]p\\d+','',idProb$aID)

#find all cases where an assembly name is duplicated
idProb <- idProb[idProb$aIDsub %in% idProb[duplicated(idProb$aIDsub),'aIDsub'],]

#for all assembly #'s that are duplicated, calculate the total
#length of coding sequence for each isoform 
idProb$len <- vapply(idProb$aID, function(x) {
  codRow <- inGff[inGff$aID == x & inGff$V3 == 'CDS',4:5]
  codRow$len <- codRow$V5 - codRow$V4
  return(sum(codRow$len))
},numeric(1))

#break up DF into sub DFs grouped by assembly #
idProb.list <- split(idProb,idProb$aIDsub)

#flag any isoform that isn't the longest for that assembly #
idProb.list <- lapply(idProb.list, function(x){
  x$drop <- (x$len != max(x$len)) | (x$len == 0)
  return(x)
})

#sometimes there are ties, so just arbitrarily pick
#one of the longest ones to keep
idProb.list <- lapply(idProb.list, function(x){
  if(length(which(x$drop == F)) > 1){
    x <- x[order(-x$len),]
    x$drop <- T
    x$drop[1] <- F
  }
  return(x)
})

idProb.df <- do.call(rbind,idProb.list)

#get the ID for all isoforms that we'll be dropping
dropThese <- idProb.df[idProb.df$drop,'aID']

#remove all redundant isoform rows
inGff <- inGff[!(inGff$aID %in% dropThese),]

write.table(inGff[,c(1:9)],file = 'clPasa.merge.longestIso.lr.gff3',sep = '\t',quote = F, row.names = F, col.names = F)
