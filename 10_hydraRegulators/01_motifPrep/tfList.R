library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

#import ipr results
ipr <- read.delim('HVAEP1.prot.longestIso.fa.tsv', header = F)

#pull initial list of candidate TFs based on either domain composition or GO annotation
ipr.tf <- ipr[grepl('GO:0003700|GO:0003677|GO:0006355|GO:0043565',ipr$V14) | grepl('helix-loop-helix|hmg|sox|winged|c2h2',ipr$V6,ignore.case = T),]

#drop genes associated with non-TF GO terms
exclIDs <- ipr[grepl('GO:0000814|GO:0008180|GO:0005840|GO:0031213|GO:0032021|GO:0006281|GO:0000150|GO:0003721|GO:0030337|GO:0006298|GO:0000723|GO:0006302|GO:0003899|GO:0006351|GO:0006313|GO:0003887|GO:0000786|GO:0003697|GO:0006260|GO:0016592|GO:0003755|GO:0006265|GO:0003684|GO:0006338|GO:0006384|GO:0006383|GO:0003723|GO:0006270|GO:0003910|GO:0006367|GO:0032508|GO:0006289',ipr$V14),1]

#drop genes associated with non-TF annotations
moreEclIDs <- ipr[grepl('Oxygenase|transferase|Transposase|phosphodiesterase|peptidase|DNA REPLICATION|Tetratricopeptide|membrane-bound|hydrolase|proteinase|mitochondrial|reductase|synthase|autointegration|nuclease|TRANSFERASE|histone|helicase|topoisomerase|DNA polymerase|rejoining|rad51|Membrane-anchored|Translin|BESS|seet_6|LIPOMA|QSOX',ipr$V6, ignore.case = T),1]

ipr.tf <- ipr.tf[!(ipr.tf$V1 %in% c(exclIDs,moreEclIDs)),]

#subset to just the IDs
tfIDs <- unique(ipr.tf$V1)

#convert transcript IDs to gene IDs
tfIDs <- gsub('HVAEP1_T(\\d+)[.]\\d+','HVAEP1_G\\1',tfIDs)

write.csv(ipr.tf,'fullTF.csv',row.names=F)
write.table(tfIDs,'tfIDs.txt',row.names = F,col.names = F,quote = F)
