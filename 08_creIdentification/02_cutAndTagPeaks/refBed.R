args = commandArgs(trailingOnly=TRUE)

inBed <- read.delim(args[1],sep = "\t",header = F)

inBed[,6] <- "."

inBed[,5] <- inBed[,4]

inBed[,4] <- 1:nrow(inBed)

newCols <- data.frame(V1 = rep(0,nrow(inBed)), V2 = rep(0,nrow(inBed)), V3 = rep(0,nrow(inBed)))

outBed <- cbind(inBed,newCols)

write.table(outBed, file = gsub(".bed",".rfmt.bed",args[1]), quote = F, row.names = F, col.names = F, sep = "\t")
