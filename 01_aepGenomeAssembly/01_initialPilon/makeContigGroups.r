scaffs <- read.table("contLengths.txt", stringsAsFactors = F)

scaffs <- scaffs[order(-scaffs$V2),]

scaffV <- c()
scaffL <- c()

totalSeq <- 0
for(i in 1:nrow(scaffs)){
  if (scaffs[i,2] >= 50000000) {
    scaffL <- c(scaffL, scaffs[i,1])
  } else if (totalSeq + scaffs[i,2] < 50000000 & i != nrow(scaffs)) {
    totalSeq <- totalSeq + scaffs[i,2]
    scaffV <- c(scaffV, scaffs[i,1])
  } else if (i == nrow(scaffs)) {
    scaffV <- c(scaffV, scaffs[i,1])
    scaffL <- c(scaffL, paste(scaffV, collapse = " "))
  } else {
    totalSeq <- scaffs[i,2]
    scaffL <- c(scaffL, paste(scaffV, collapse = " "))
    scaffV <- c(scaffs[i,1])
  }
}

write.table(scaffL, file = "contGroups.txt", quote = F, sep = " ", row.names = F, col.names = F)
