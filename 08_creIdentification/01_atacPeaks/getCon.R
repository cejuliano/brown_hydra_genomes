int12 <- read.delim("12Int.bed",header = F, sep = "\t")
int13 <- read.delim("13Int.bed",header = F, sep = "\t")
int23 <- read.delim("23Int.bed",header = F, sep = "\t")
consensus <- int12[(int12[,11] + int13[,11] + int23[,11]) > 1,1:6]
write.table(consensus, file = "consensusAEP.bed", sep = '\t', quote = F, col.names = F, row.names = F)
