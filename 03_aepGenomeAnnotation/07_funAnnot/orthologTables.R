library(Biostrings)
library(plyr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#function to expand orthology table so that each AEP gene gets its own row
hExpand <- function(x){
  x <- paste0("Results_Sep15_1/Orthologues/Orthologues_H_vulgarisAEP/H_vulgarisAEP__v__",x,".tsv")
  specOrthos <- read.delim(x, sep = "\t", stringsAsFactors = F, header = T)
  s <- strsplit(specOrthos$H_vulgarisAEP, split = ", ")
  specOrthos <- data.frame(Orthogroup = rep(specOrthos$Orthogroup, sapply(s, length)), H_vulgarisAEP = unlist(s), Ortholog = rep(specOrthos[,3], sapply(s, length)))
  return(specOrthos)
}

#define the species that we want to pull functional annotations from based on orthology (picked well annotated systems)
#the order of these species also defines the priority that they receive when assigning function
#as in: if there's an ortholog from both humans and flies, we'll pick the human one for prediction function
orthoLists <- lapply(c("H_sapiens", "M_musculus","X_tropicalis", "D_melanogaster", "C_elegans"), function(x) hExpand(x))

#collapse all orthology assignments to a single table
orthoLists <- do.call(rbind,orthoLists)

#needed to add this because of one weirdly named drosophila gene
orthoLists$Ortholog <- gsub('_(.)_(.)_','_\\1\\2',orthoLists$Ortholog)

#drop redundant orthology assignments using the species prioritization defined above so that we only use
#one orthology assignment per Hydra gene
orthoLists <- orthoLists[!duplicated(orthoLists$H_vulgarisAEP),]

#only keep orthologs that had a gene name attached
orthoLists <- orthoLists[grepl("_",orthoLists$Ortholog),]

#split gene names from ensembl IDs
orthoLists.rfmt <- strsplit(orthoLists$Ortholog,', ')

#catch any remaining orthologs without a name
orthoLists.rfmt <- lapply(orthoLists.rfmt,function(x) x[grepl('_',x)])

orthoLists.rfmt.IDs <- orthoLists.rfmt

#isolate just the gene names
orthoLists.rfmt <- lapply(orthoLists.rfmt,function(x) gsub('.*_','',x))

#isolate just the ensembl IDs
orthoLists.rfmt.IDs <- lapply(orthoLists.rfmt.IDs,function(x) gsub('_.*','',x))

#bring in long gene names as well (more readable)
descTab <- c('hsapiens.csv','mmusculus.csv','xtropicalis.csv','dmelanogaster.csv','celegans.csv')

descTab <- vapply(descTab,function(x) paste0('../proteomes/ensemblAnnotation/',x),'')

#import ensembl gene info tables
descTab <- lapply(descTab,read.csv)

#collapse into a single table
descTab <- do.call(rbind,descTab)

#pull the long names associated with each Hydra gene in our ortholog table
descList <- lapply(orthoLists.rfmt.IDs,function(x) descTab[descTab$ensembl_gene_id %in% x,4])

#drop redundant long names
descList <- lapply(descList,unique)

#collapse short ortholog names into single string per Hydra gene
orthoLists.rfmt <- lapply(orthoLists.rfmt,function(x) paste(x,collapse = ', '))

#collapse ensembl IDs into single string for each Hydra gene
orthoLists.rfmt.IDs <- lapply(orthoLists.rfmt.IDs,function(x) paste(x,collapse = ', '))

#collapse long ortholog names into single string per Hydra gene
descList <- lapply(descList,function(x) paste(x,collapse = ', '))

#add separate ensembl ID, short name, and long name columns to ortholog table
orthoLists$Ortholog <- unlist(orthoLists.rfmt)
orthoLists$EnsemblID <- unlist(orthoLists.rfmt.IDs)
orthoLists$EnsemblLongName <- unlist(descList)

#this function attempts to automate the process of collapsing gene IDs from the same gene family into a single ID
#that retains accurate information about orthology
#it sort of works...
collapseOrthos <- function(x){
  
  #initialize vector where collapsed ortholog names will be placed
  ortho.annot <- character(0)
  
  #iterate through each set of orthologs (one set per Hydra gene)
  for (batch in x){
    
    #split each name into it's own character string
    batch <- strsplit(batch,", ")[[1]]
    
    #if there's only a single ortholog, there's nothing to collapse,
    #so just add it to the results and move on
    if (length(batch) == 1){
      ortho.annot <- c(ortho.annot,batch)
      next
    }
  
    #try to find the longest common prefixes shared across orthologs
    #this is aimed at identifying common gene names to collapse
    #into a single more compact name
    
    #initialize empty results vector
    ind.lcp <- character(0)
    
    #loop through each ortholog name and find
    #the longest common prefix (lcp) shared between
    #that gene and any of the other remaining
    #gene names
    for (i in 1:length(batch)){
      
      #current gene name being considered
      batchword <- batch[i]
      
      #all other ortholog names
      restword <- batch[-i]
      
      #variable to hold lcp candidates
      temphit <- ""
      
      #loop through other gene names and pull the lcp for each
      #if it beats the current lcp, assign it as the new lcp
      for (j in 1:length(restword)){
        subbatchword <- restword[j]
        lcpres <- substr(subbatchword, start = 1, stop = lcprefix(batchword,subbatchword))
        if(nchar(lcpres) > nchar(temphit)){
          temphit <- lcpres
        }
      }
      
      #after looping through all names, report the longest lcp
      ind.lcp <- c(ind.lcp, temphit)
    }
    
    #if the lcp is short or empty it was likely a singleton
    #meaning it can't be collapsed, so just go with the original name
    for (i in 1:length(ind.lcp)) {
      if(nchar(ind.lcp[i]) < 3){
        ind.lcp[i] <- batch[i]
      }
      if(ind.lcp[i] == ""){
        ind.lcp[i] <- batch[i]
      }
    }
    
    #collapse repeated lcps
    ind.lcp <- unique(ind.lcp)
    
    #drop any numbers that are likely used to refer to orthologs within an annotation set
    #(e.g., extract wnt from wnt8a and wnt8b)
    ind.lcp.pref <- unique(gsub("\\d+$|\\d+.$","",ind.lcp))
    
    #initialize empty results object
    ortho.collapse <- character(0)
    
    #for each gene prefix, collapse all suffixes into a single label
    #(e.g., for prefix DKK, collapse DKK1, DKK2, DKK3, and DKK4 into DKK1/2/3/4)
    for (pref in ind.lcp.pref) {
      #delete the prefix from the gene names
      #this should give all of the ortholog specific numbers/letters
      #(e.g., 8a and 8b for wnt8a and wnt8b)
      ind.lcp.pref.n <- gsub(pref,"",batch[grepl(pref,batch)])
      
      #drop cases where the prefix encompassed the entire gene name
      ind.lcp.pref.n <- ind.lcp.pref.n[ind.lcp.pref.n != ""]
      
      #drop letter suffixes in cases where there is a number then a letter
      #(e.g., convert 8a and 8b to just 8)
      ind.lcp.pref.n <- unique(gsub("(\\d+)\\D$","\\1",ind.lcp.pref.n))
      
      #order suffixes
      ind.lcp.pref.n <- ind.lcp.pref.n[order(ind.lcp.pref.n)]
      
      #collapse into a single string
      ind.lcp.pref.n <- paste(ind.lcp.pref.n, collapse = "/")
      
      #add on prefix
      ind.lcp.pref.n <- paste0(pref,ind.lcp.pref.n)
      
      #add to results list
      ortho.collapse <- c(ortho.collapse,ind.lcp.pref.n)
    }
    
    #after all prefixes have been processed, collapse compact names into single character string
    ortho.collapse <- paste(ortho.collapse, collapse = "; ")
    
    #append final collapsed string to results vector
    ortho.annot <- c(ortho.annot,ortho.collapse)
  }
  return(ortho.annot)
}

orthoLists$OrthologCollapse <- collapseOrthos(orthoLists$Ortholog)

#also bring in blast hits for manually deposited Hydra genbank sequences
gb <- read.delim("/Users/Jcazet/Google_Drive/Juliano_lab/References/genbank/gb2AEP.txt", sep = "\t", stringsAsFactors = F, header  = F)

#drop tailing pipe from accession numbers
gb$V2 <- gsub("[|]","",gb$V2)

#in cases where a genbank entry had multiple hits in the AEP gene models
#just go with the one with the best alignment score
gb <- gb[order(-gb$V12),]

gb <- gb[!duplicated(gb$V2),]

#also just pick a single genbank entry (best hit) for each AEP gene
gb <- gb[!duplicated(gb$V1),]

#subset to just include IDs to queries and hits
gb <- gb[,1:2]

#pull in the full fasta headers from the genbank entries
gb.head <- read.delim("/Users/Jcazet/Google_Drive/Juliano_lab/References/genbank/headers.txt", stringsAsFactors = F, header  = F)

#reformat to get just the gene name information from the header
gb.head$name <- gb.head[,1]
gb.head$name <- gsub(">.*?[|] ","",gb.head$name)
gb.head$name <- gsub("[|].*$","",gb.head$name)
gb.head$name <- gsub("Hydra \\D+? ","",gb.head$name)
gb.head$name <- gsub("strain \\w+? ","",gb.head$name)
gb.head$name <- gsub("isolate \\D+? ","",gb.head$name)
gb.head$name <- gsub(" mRNA","",gb.head$name)
gb.head$name <- gsub("\\D+? for ","",gb.head$name)

#get column of just the accession
gb.head$V1 <- gsub("[|].*","",gb.head$V1)
gb.head$V1 <- gsub(">","",gb.head$V1)

#add header info to blast hit table
gb <- merge(gb, gb.head, by.x = "V2" ,by.y = "V1", all.x = T)

#combine accession and gene name
gb$V2 <- paste0(gb$V2,'; ',gb$name)

gb <- gb[,1:2]

colnames(gb) <- c("genBankAnnotation","ID")

#add genbank hits to ortholog table
orthoLists <- merge(orthoLists, gb, by.x = "H_vulgarisAEP",by.y = "ID",all = T)

#bring in interpro protein domain predictions
ipr <- read.delim("../Genome_annotation/functionalAnnotation/HVAEP1.prot.longestIso.fa.tsv", header = F)

#extract PANTHER annotations
ipr.p <- ipr[ipr$V4 == "PANTHER",]

#drop entries without description
ipr.p <- ipr.p[ipr.p$V6 != "-",]

#remove redundant rows
ipr.p <- ipr.p[!duplicated(paste(ipr.p$V1,ipr.p$V6)),]

ipr.p <- ipr.p[!(ipr.p$V5 %in% gsub(":.*","",ipr.p[grepl(":",ipr.p$V5),"V5"])),]

ipr.p <- ipr.p[!duplicated(ipr.p$V1),]

#subset to just include the most useful columns
ipr.p <- ipr.p[,c(1,5,6)]

colnames(ipr.p) <- c("ID","PANTHER_ID","PANTHER_NAME")

#add PANTHER annotations to ortholog table
orthoLists <- merge(orthoLists, ipr.p, by.x = "H_vulgarisAEP", by.y = "ID", all = T)

#extract pfam predictions from interpro results
ipr.pf <- ipr[ipr$V4 == "Pfam",]

#drop redundant entries (cases where a domain shows up multiple times in a gene)
ipr.pf <- ipr.pf[!duplicated(paste(ipr.pf$V1,ipr.pf$V6)),]

#split table into subtables grouped by gene
ipr.pf.ls <- split(ipr.pf,ipr.pf$V1)

#set aside gene IDs
ipr.pf.ls.nm <- names(ipr.pf.ls)

#collapse the list of pfam domain names into a single character string
ipr.pf.ls.ds <- vapply(ipr.pf.ls, function(x) paste(x[,6],collapse = "; "),"")

#collapse the list of pfam domain IDs into a single character string
ipr.pf.ls.id <- vapply(ipr.pf.ls, function(x) paste(x[,5],collapse = "; "),"")

#generate pfam annotation table consisting of Hydra gene name, all
#associated pfam domain IDs, and corresponding pfam names
ipr.pf <- data.frame(ID = ipr.pf.ls.nm, PFAM_ID = ipr.pf.ls.id, PFAM_NAME = ipr.pf.ls.ds)

#add pfam annotations to ortholog list
orthoLists <- merge(orthoLists, ipr.pf, by.x = "H_vulgarisAEP", by.y = "ID", all = T)

#import uniprot blast hits
upHits <- read.delim('../proteomes/upBlast/uniprotBlast.txt',header = F)[,1:2]

colnames(upHits) <- c("ID",'UniprotHit')

#add uniprot hits to ortholog table
orthoLists <- merge(orthoLists, upHits, by.x = "H_vulgarisAEP", by.y = "ID", all = T)

write.csv(orthoLists,"HVAEP1_annotation.csv", row.names = F)