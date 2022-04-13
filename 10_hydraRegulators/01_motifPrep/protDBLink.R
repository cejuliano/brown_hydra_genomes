library(rstudioapi)
library(jsonlite)
library(httr)
library(RCurl)
setwd(dirname(getActiveDocumentContext()$path))

#get full list of jaspar motif IDs
motif.IDs <- read.delim('jasparIDs.txt')[,1]

results <- list(length(motif.IDs))

#get TF name and family for each motif ID
for (i in 1:length(motif.IDs)) {
  print(i)
  link <- paste0("http://jaspar2020.genereg.net/api/v1/matrix/",as.character(motif.IDs[i]),"/")
  result <- fromJSON(url(link))
  result.n <- paste(result[[c("name")]],collapse = ', ')
  result.f <- paste(result[[c("family")]],collapse = ', ')
  results[[i]] <- c(result.n,result.f)
}

#collapse results into DF
motInfo <- data.frame(ID=motif.IDs,
                      name=vapply(results,function(x) x[1],""),
                      family=vapply(results,function(x) x[2],""))

write.csv(motInfo,'motifInfo.csv')

#get uniprot IDs for each motif from JASPAR
results <- list(length(motif.IDs))

for (i in 1:length(motif.IDs)) {
  print(i)
  link <- paste0("http://jaspar2020.genereg.net/api/v1/matrix/",as.character(motif.IDs[i]),"/")
  result <- fromJSON(url(link))
  result <- result[["uniprot_ids"]]
  if(length(result != 0)) {
    results[[i]] <- result
  } else {
    print("Empty")
  }
}

#pull names for each result
names(results) <- motif.IDs

pfamAnnot <- list()
#Pull the pfam domains assocaited with each swissprot entry
for (i in 1:length(results)) {
  print(i)
  subRes <- results[[i]]
  enGene <- c()
  for(ID in subRes){
    print(ID)
    link <- paste0("https://www.uniprot.org/uniprot/",ID,".txt")
    result <- GET(link)
    result <- rawToChar(result$content)
    result <- strsplit(result, "\n")
    result <- result[[1]]
    result <- result[grepl("Pfam", result)]
    result <- strsplit(result,'; ')
    result <- lapply(result,function(x) x[2])
    result <- unlist(result)
    pfamAnnot[[i]] <- unique(c(enGene,result))
  }
}

names(pfamAnnot) <- names(results)

#pull ipr results
ipr <- read.delim('../../Genome_annotation/functionalAnnotation/HVAEP1.prot.longestIso.fa.tsv', header = F)
ipr$V1 <- gsub('_T(\\d+)[.]\\d+','_G\\1',ipr$V1)

#import TF IDs
tfID <- read.delim('../../Genome_annotation/functionalAnnotation/tfIDs.txt',header=F)[,1]

#subset ipr table to include only TFs
ipr.tf <- ipr[ipr$V1 %in% tfID,]

#get PF domains associated with all the TF entries in the ipr table
ipr.tf.pf <- unique(ipr.tf[ipr.tf$V4 == 'Pfam',5])

#drop pfam domains from swissprot that didn't show up in our ipr domain list from our hydra TFs
pfamAnnot <- lapply(pfamAnnot,function(x) x[x %in% ipr.tf.pf])

#get the Hydra genes that contain the domains linked to each motif
pfamAnnot.genes <- lapply(pfamAnnot,function(x) unique(ipr[ipr$V5 %in% x,1]))

#drop motifs without a linked gene
pfamAnnot.genes <- pfamAnnot.genes[sapply(pfamAnnot.genes,length) > 0]

saveRDS(pfamAnnot.genes,'motifBindPfam.rds')
