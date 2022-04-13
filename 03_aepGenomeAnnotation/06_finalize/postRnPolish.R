setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

inG <- read.delim('HVAEP1.PU.RN.gff3', skip = 1, header = F)

inG$V9 <- gsub('Name=;','',inG$V9)

inG$V9 <- gsub('Alias=[^;]+;','',inG$V9)

inG$V9 <- gsub('deletions=[^;]+;','',inG$V9)

inG$gID <- gsub('.*HVAEP1_[TG](\\d+).*','\\1',inG$V9)

inG$V9 <- gsub('gene_id=.*$','',inG$V9)

inG$V9 <- paste0(inG$V9,'gene_id=HVAEP1_G',inG$gID)

inG$tID <- gsub('.*(HVAEP1_T\\d+[.]\\d+);.*','\\1',inG$V9)

inG$euFix <- inG$V9

inG$euFix <- gsub(';.*','',inG$euFix)

inG$euFix <- gsub('ID=T\\d+[.]\\d+(\\D+)','\\1',inG$euFix)
inG$euFix <- gsub('ID=T\\d+_\\d+_[^.]+[.](\\D+)','.\\1',inG$euFix)
inG$euFix <- gsub('ID=split[^.]+[.][^.]+','',inG$euFix)
inG$euFix <- gsub('_T\\d+[.]\\d+[.](\\D+)','.\\1',inG$euFix)
inG$euFix <- gsub('ID=file_1_file_1_jg\\d+[.]t\\d+(\\D+)','.\\1',inG$euFix)
inG$euFix <- gsub('ID=t\\d+aep-(\\D+)-(\\d+)','.\\1.\\2',inG$euFix)
inG$euFix <- gsub('ID=nbis-','.',inG$euFix)
inG$euFix <- gsub('five_prime_utr[\\.-]','utr5p',inG$euFix)
inG$euFix <- gsub('three_prime_utr[\\.-]','utr3p',inG$euFix)

inG$euFix <- paste0('ID=',inG$tID,inG$euFix)

targRow <- grepl('UTR',inG$V3) | grepl('exon',inG$V3)

inG[!targRow,'euFix'] <- ''

inG[targRow,9] <- gsub('ID=[^;]+','',inG[targRow,9])

inG[targRow,9] <- paste0(inG[targRow,'euFix'],inG[targRow,9])

inG$euFix <- NULL


inG$cFix <- inG$tID

inG[inG$V3 != 'CDS','cFix'] <- ''

rleRes <- data.frame(lengths = rle(inG$cFix)$lengths, values = rle(inG$cFix)$values)

rleRes <- lapply(rleRes$lengths, function(x) seq(from=1, to=x, by = 1))

rleRes <- do.call(c,rleRes)

inG$cFix <- paste0('.',rleRes,'.',inG$cFix)

inG[inG$V3 != 'CDS','cFix'] <- ''

inG[inG$V3 == 'CDS' ,9] <- gsub('ID=[^;]+','',inG[inG$V3 == 'CDS' ,9])
inG[inG$V3 == 'CDS' ,9] <- paste0('ID=cds',inG[inG$V3 == 'CDS' ,'cFix'],inG[inG$V3 == 'CDS' ,9])

inG$cFix <- NULL

inG$tID <- paste0('transcript_id=',inG$tID)

inG[inG$V3 == 'gene','tID'] <- ''

inG$V9 <- paste0(inG$V9,';',inG$tID)

inG$V9 <- gsub(';$','',inG$V9)

write.table(inG[,1:9],file = 'HVAEP1.PU.RN.pol.gff3', sep = '\t', row.names = F, col.names = F, quote = F)
