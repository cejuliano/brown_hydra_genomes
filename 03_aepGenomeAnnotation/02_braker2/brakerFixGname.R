library(rstudioapi)

#set the working directory to be the folder in which this script is located
setwd(dirname(getActiveDocumentContext()$path))

gffIn <- read.delim("braker.gff3", header = F)

#pull out only those rose that are from genemark predictions (these are the problem ones)
gffIn.GM <- gffIn[grepl('GeneMark.hmm',gffIn$V2),]

#id the parent gene for each row
gffIn.GM$parent <- gsub('.*Parent=','',gffIn.GM$V9)

#remove gene models with incomplete ORFs that lack a start codon
gffIn.GM <- gffIn.GM[gffIn.GM$parent %in% gffIn.GM[gffIn.GM$V3 == 'start_codon','parent'],] 

#split dataframe by parent ID (groups rows by gene model)
gffIn.GM <- split(gffIn.GM,gffIn.GM$parent)

#now we need to add an mRNA row for each gene model
gffIn.GM <- lapply(gffIn.GM, function(x) {
  newDF <- x
  #pick an arbitrary row that we will remake into an mRNA row
  mRow <- newDF[1,]
  #extract parent ID for gene model in this DF
  #this is the basis for the gene name
  pID <- gsub(';','',newDF[1,'parent'])
  #set the start and end to encompass the full span of all rows for this gene prediction
  mRow[,4] <- min(newDF[,4])
  mRow[,5] <- max(newDF[,5])
  #rename feature type
  mRow[,3] <- 'mRNA'
  mRow[,6] <- '.'
  mRow[,8] <- '.'
  #give the mRNA row the proper mRNA ID in the tags column
  mRow[,9] <- paste0('ID=',pID,'-mRNA-1;Parent=',pID)
  #we also need to create a gene row
  #(the GM predictions only have one isoform, so the gene entry is identical to the mRNA entry)
  gRow <- mRow
  gRow[,9] <- paste0('ID=',pID)
  gRow[,3] <- 'gene'
  #make sure all the rows for CDS, exons, etc. have the new transcript ID for their parent tag
  newDF[,9] <- gsub(paste0('Parent=',pID),paste0('Parent=',pID,'-mRNA-1'),newDF[,9])
  #combine everything
  newDF <- rbind(gRow,mRow,newDF)
  #link all the rows with the gene ID
  newDF[,9] <- paste0(newDF[,9],'gene_id=',pID)
  #delete temp parent column
  newDF$parent <- NULL
  return(newDF)
})

#compile the fixed gene models into a DF
gffIn.GM <- do.call(rbind,gffIn.GM)

#drop the old, improperly formated versions of the GM models
gffIn <- gffIn[!grepl('GeneMark.hmm',gffIn$V2),]

#all the Augustus models lack a gene row, and some of them also lack a mRNA row
#so we need to fix that too

#extract the gene IDs for all rows
gffIn$gID <- gffIn$V9

gffIn$gID <- gsub(';.*','',gffIn$gID)
gffIn$gID <- gsub('[.]t.*','',gffIn$gID)
gffIn$gID <- gsub('ID=','',gffIn$gID)

#there are some weird rows that I don't understand, so let's just drop them
gffIn <- gffIn[!(gffIn$V3 %in% c('initial','terminal','internal')),]

#split augustus models into a list of DFs grouped by gene ID
gffInList <- split(gffIn, gffIn$gID)

gffInList <- lapply(gffInList, function(x) {
  old.df <- x
  #check if the augustus models have an mRNA row
  #if they don't, add them (same approach as above)
  if(length(old.df[old.df$V3 == 'mRNA',1]) == 0) {
    mRow <- old.df[1,]
    mRow[,3] <- 'mRNA'
    mRow[,4] <- min(old.df$V4)
    mRow[,5] <- max(old.df$V5)
    mRow[,6] <- '.'
    mRow[,8] <- '.'
    mRow[,9] <- paste0('ID=',old.df[1,10],'.t1;Parent=',old.df[1,10],';')
    
    #also add a gene row
    gRow <- mRow
    gRow[,9] <- paste0('ID=',old.df[1,10],';')
    gRow[,3] <- 'gene'
    
    new.df <- rbind(gRow,mRow,old.df)
  } else {  #if there's already an mRNA row, then just add the gene row
    gRow <- old.df[old.df$V3 == 'mRNA',][1,]
    gRow[,9] <- paste0('ID=',old.df[1,10],';')
    gRow[,3] <- 'gene'
    gRow[,4] <- min(old.df$V4)
    gRow[,5] <- max(old.df$V5)
    
    new.df <- rbind(gRow,old.df)
  }
  new.df$V9 <- paste0(new.df$V9,'gene_id=',old.df[1,10],';')
  return(new.df)
})

#re-order the gene models so they are consecutive
gffInList <- gffInList[order(as.numeric(gsub('.*g','',names(gffInList))))]

#combine list of DFs into single DF
newGff <- do.call(rbind,gffInList)
#drop temp gene ID row
newGff$gID <- NULL

#combine augustus and GM gene models
newGff <- rbind(newGff,gffIn.GM)

#export reformated GFF3
write.table(newGff, file = 'braker.fix.gff3', row.names = F, col.names = F, sep = '\t', quote = F)

