library(rstudioapi)
library(biomaRt)
library(Biostrings)
library(plyr)

setwd(dirname(getActiveDocumentContext()$path))

#get the list of all individual proteomes to be used for orthofinder analysis
seqList <- list.files("individual/primary_transcripts/", full.names = T)

#import AA sequences
seqs <- lapply(seqList, function(x) readAAStringSet(x))

#extract species names from filenames
seqNames <- gsub(".*/","",seqList)
seqNames <- gsub(".fa","",seqNames)
seqNames <- tolower(seqNames)
seqNames <- gsub("_","",seqNames)

names(seqs) <- seqNames

#generate ensembl database name from species name
setName <- vapply(seqNames, function(x) paste0(x,"_gene_ensembl"), character(1))

#check to find which species have databases available
ensembList <- lapply(setName, function(x) try(useEnsembl(biomart = "genes", dataset = x), silent = T))

#subset to only include species with an ensembl db hit
enSubset <- vapply(ensembList, function(x) !is.character(x), logical(1))

ensembList <- ensembList[enSubset]

seqs <- seqs[enSubset]

#initialize empty results object
annots <- list()

#for each proteome, download ensembl annotation data
#gene name, go terms, description, uniparc ID
for (i in 1:length(seqs)){

    res <- getBM(attributes = c("ensembl_gene_id","external_gene_name","external_gene_source","go_id","entrezgene_description","uniparc"),
                 filters = "ensembl_gene_id",
                 values = substr(gsub("[.].*","",seqs[[i]]@ranges@NAMES),1,18),
                 mart = ensembList[[i]])
    
  #drop any characters after a space in gene name
  res[,2] <- gsub(" .*$","",res[,2])
  annots[[i]] <- res
}

#generate ID that combines ensembl ID and gene name
#this will be the ID that's used to replace the AA fasta header
annots <- lapply(annots, function(x) cbind(x,finAnnot = paste(x[,1],x$external_gene_name, sep = "_")))

#subset to annots to just be a conversion table from old to new IDs
finAnnots <- lapply(annots, function(x) unique(x[,c(1,7)]))

#make sure to fix any instances where there was no gene name to append
finAnnots <- lapply(finAnnots, function(x) {
  x[,2] <- gsub('_$','',x[,2])
  return(x)
})

#replace the old names on the AA with the new ones that have the gene name included
newSeqs <- lapply(1:length(seqs), function(x) {
  newSeqObj <- seqs[[x]]
  newSeqObj@ranges@NAMES <- mapvalues(substr(gsub("[.].*","",newSeqObj@ranges@NAMES),1,18), 
                                      from = finAnnots[[x]][,1],
                                      to = gsub("_$","",finAnnots[[x]][,2]))
  return(newSeqObj)
})

#generate modified filenames for the output so as not to overwrite the original files
newFileNames <- vapply(names(seqs), function(x) {
  fChar <- toupper(substr(x,1,1))
  lChar <- substr(x,2,nchar(x))
  return(paste0(fChar,'_',lChar))
}, "")

#export new AA fastas with updated IDs
lapply(1:length(newSeqs), function(i) {
  writeXStringSet(newSeqs[[i]],
                  paste0("individual/primary_transcripts/",newFileNames[i],".fa"))
})

#subset annots to just be ensembl ID and GO terms
goCollapse <- lapply(annots, function(x) x[,c(1,4)])

#collapse GO terms by gene ID
goCollapse <- lapply(goCollapse, function(x) aggregate(x[,2], by = list(x[,1]), paste, collapse = ";"))

#fix cases where empty results were aggregated (creating things like ';;')
goCollapse <- lapply(goCollapse, function(x) {
  newDF <- x
  newDF[,2] <- gsub(';$|;;+|^;','',newDF[,2])
  return(newDF)
})

annots.rfmt <- annots

#reduce annotation tables to one row per gene
annots.rfmt <- lapply(annots.rfmt, function(x) unique(x[,-4]))

#replace old GO column with new, collapsed set of all GO terms for each gene
for(i in 1:length(annots.rfmt)) {
  tmpDF <- annots.rfmt[[i]]
  tmpDF$go <- mapvalues(tmpDF[,1], from = goCollapse[[i]][,1], to = goCollapse[[i]][,2], warn_missing = F)
  annots.rfmt[[i]] <- tmpDF[,-6]
}

#export table on gene functional data for each species
dir.create('ensemblAnnotation',showWarnings = F)

lapply(1:length(annots.rfmt), function(i) {
  write.csv(annots.rfmt[[i]],
                  paste0("ensemblAnnotation/",names(seqs)[i],".csv"), row.names = F)
})
