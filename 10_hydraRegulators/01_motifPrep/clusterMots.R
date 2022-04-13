library(rstudioapi)

setwd(dirname(getActiveDocumentContext()$path))

#import similarity matrix for JASPAR animal PWMs
similarityMat <- read.delim("motComp.txt", row.names = 1, header = T, sep = "\t")

#reformat motif names (rows and columns)
newNames <- rownames(similarityMat)
newNames <- gsub('/Jaspar','',newNames)
newNames <- gsub('/','_',newNames)

rownames(similarityMat) <- newNames
colnames(similarityMat) <- newNames

#generate distance matrix from correlation coefficient
d <- as.dist(1 - similarityMat)

#use hierarchical clustering to group motifs based on correlation
hc1 <- hclust(d, method = "average")

pdf('motDistClust.pdf',width=150,height=20)
plot(hc1)
rect.hclust(hc1 , h = 0.3)
abline(h = 0.3, col = 'red')
dev.off()

#cut the branches at 0.3 and group accordingly
clusts <- cutree(hc1, h = 0.3)
length(unique(clusts))

#save cluster information
repMotifs <- data.frame(ID = names(clusts), clust = clusts)

write.csv(repMotifs, file = "motif_clusters.csv")
