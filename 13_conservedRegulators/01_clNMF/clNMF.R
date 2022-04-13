library(Seurat)
library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

cl <- readRDS('initSeurat.rds')

####NMF prep####
rawC <- t(as.matrix(cl@assays$RNA@counts))
write.table(rawC,file="cl.raw.counts.tsv",sep = '\t', quote = F)
rm(rawC)

normC <- t(as.matrix(cl@assays$SCT@data))
write.table(normC,file="cl.norm.counts.tsv",sep = '\t', quote = F)
rm(normC)

#list of variable genes to look at
write.table(cl@assays$SCT@var.features, file = 'cl.genes.tsv', row.names = F, col.names = F, quote = F)
