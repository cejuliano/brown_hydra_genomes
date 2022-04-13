args <-  commandArgs(trailingOnly=TRUE)

prefix <- args[1]

intRes <- lapply(args[-1],read.delim,header=F,sep='\t')

intRes.sub <- lapply(intRes,function(x) x[,ncol(x)])

intRes.sub <- do.call(cbind,intRes.sub)

intRes.sub <- apply(intRes.sub,1,sum)

consensus <- intRes[[1]][intRes.sub > 1,1:6]

# consensus$V4 <- gsub('.*_MG_','',consensus$V4)

write.table(consensus, file = paste0("consensus",prefix,".bed"), sep = '\t', quote = F, col.names = F, row.names = F)
