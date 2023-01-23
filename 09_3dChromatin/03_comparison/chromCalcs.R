library(ggplot2)
library(zoo)

#####cross chrom cent search#####

centCrossCheck <- function(specCheck,chrCheck){
  print(chrCheck)
  
  #get chromosome length for target chromosome
  chrSize <- cSizes[cSizes$V1 == chrCheck,'V2']
  
  #read in all interaction matrices for all cross-chromosome pairs that include target chr
  crossChrs <- list.files(path=specCheck,pattern = paste0(chrCheck,'_.*txt'), full.names = T)
  crossChrs <- lapply(crossChrs, read.delim, header=F, skip=1)
  
  #add in names so that each list entry is labeled with which chr pair is being quantified
  names(crossChrs) <- list.files(path=specCheck,pattern = paste0(chrCheck,'_.*txt'))
  
  #the target chromosome isn't always assigned to the x-axis (column V1)
  #sometimes it's the y-axis (column V2)
  
  #Here I'm looking at which column has a maximum coordinate value that is closest to
  #the length of the target chromosome in order to determine which column contains
  #the coordinates for the target chromosome
  colUse <- lapply(crossChrs, function(x) sapply(x, function(y) abs(max(y) - chrSize)))
  colUse <- sapply(colUse, function(x) which.min(x))
  
  # split the interaction frequency values according to their position along the target
  #chromosome
  crossChrs.sp <- lapply(1:length(crossChrs), function(x) {
    split(crossChrs[[x]],crossChrs[[x]][,colUse[x]])
  })
  
  #get the total number of interactions for each position along the target chromosome
  crossChrs.sp <- lapply(crossChrs.sp, function(x) sapply(x, function(y) sum(y$V3)))
  
  #check and make sure each cross-chromosome interaction pair has the same length
  #(sometimes it isn't)
  print(sapply(crossChrs.sp,length))
  
  crossChrs.df <- do.call(cbind,crossChrs.sp)
  
  #in cases 
  crossChrs.df <- crossChrs.df[rowSums(is.na(crossChrs.df)) != ncol(crossChrs.df),]
  
  crossChrs.ave <- apply(crossChrs.df,1,median, na.rm = T)
  
  #exclude the tips of each chromosome (to exclude telomeres)
  crossChrs.ave.trim <- crossChrs.ave[floor(length(crossChrs.ave)/10):floor(length(crossChrs.ave) - length(crossChrs.ave)/10)]
  
  crossChrs.zs <- (crossChrs.ave.trim - mean(crossChrs.ave.trim))/sd(crossChrs.ave.trim)
  
  plot(1:length(crossChrs.ave.trim), crossChrs.zs, type = 'l')
  
  print(max(crossChrs.zs))
}

specCheck <- 'aep'

cSizes <- read.delim(paste0(specCheck,'/', specCheck,'.genome'), header=F)

aepScores <- sapply(cSizes$V1, function(x) centCrossCheck(specCheck,x))

specCheck <- 'dili'

cSizes <- read.delim(paste0(specCheck,'/', specCheck,'.genome'), header=F)

diliScores <- sapply(cSizes$V1, function(x) centCrossCheck(specCheck,x))

specCheck <- 'hoct'

cSizes <- read.delim(paste0(specCheck,'/', specCheck,'.genome'), header=F)

hoctScores <- sapply(cSizes$V1, function(x) centCrossCheck(specCheck,x))

specCheck <- 'nvec200'

cSizes <- read.delim(paste0(specCheck,'/', specCheck,'.genome'), header=F)

nemScores <- sapply(cSizes$V1, function(x) centCrossCheck(specCheck,x))

specCheck <- 'resc'

cSizes <- read.delim(paste0(specCheck,'/', specCheck,'.genome'), header=F)

rescScores <- sapply(cSizes$V1, function(x) centCrossCheck(specCheck,x))

specCheck <- 'amil'

cSizes <- read.delim(paste0(specCheck,'/', specCheck,'.genome'), header=F)

amilScores <- sapply(cSizes$V1, function(x) centCrossCheck(specCheck,x))

plotDF <- data.frame(spec = rep(c('aep','dili','hoct','nem','resc','amil'),c(15,16,9,15,21,14)), 
                     scores = c(aepScores,diliScores,hoctScores,nemScores,rescScores,amilScores))

plotDF$spec <- factor(plotDF$spec, levels=c('aep', 'resc', 'hoct', 'nem', 'dili','amil'))

ggplot(plotDF, aes(x = spec, y = scores, fill=spec)) + geom_violin() + geom_jitter(width = 0.2) + theme_bw()
ggsave('interCentScores.pdf',width = 8, height = 4)

#significance test
cent.lm <- lm(scores ~ spec, data = plotDF)

cent.av <- aov(cent.lm)

tukey.test <- TukeyHSD(cent.av)

tukey.test


#ACA telomere scores

tScores <- data.frame(spec = c('aep','amil','dili','hoct','nem','resc'), 
                      score = c(1.618,1.588,1.283,1.809,1.141,1.183))

tScores$spec <- factor(tScores$spec, levels=c('aep', 'resc', 'hoct', 'nem', 'dili','amil'))

ggplot(tScores, aes(x=spec,y=score, fill=spec)) + geom_col() + theme_bw()
ggsave('interTelScores.pdf',width = 6, height = 4)



