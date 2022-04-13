args <- commandArgs(trailingOnly = T)

inCoords <- read.delim(args[1], header = F, sep = " ", stringsAsFactors = F)

inCoords$V3 <- NULL

inCoords$V2 <- inCoords$V2 + 3

inCoords$V2 <- inCoords$V2 - inCoords$V1

write.table(inCoords, file = args[1], sep = " ", quote = F, row.names = F, col.names = F)