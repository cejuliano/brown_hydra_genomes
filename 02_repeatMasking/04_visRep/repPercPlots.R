library(rstudioapi)
library(ggplot2)

setwd(dirname(getActiveDocumentContext()$path))

#import repeatmasker tables
#(results from both eumetazoa lib and AEP/olig repeatmodeler libs)
fileNames <- c('bothMaskFull.tbl','bothMaskFull105.tbl','olig_genome_combined.fa.tbl')

#function to extract percentage from character string and convert it to numeric data
getPerc <- function(x){
  return(as.numeric(gsub('.*?([^ ]+) \\%$','\\1',x)))
}

#function to extract repeat statistics from repeatmasker tables
percPlot <- function(x){
  
  #import results table text as a character vector
  #one string per line
  rp <- readLines(x)
  
  #extract the percentage of total bases masked
  tmask <- as.numeric(gsub('.*\\((.*)\\%\\).*','\\1',rp[6]))
  
  #extract percentages of masked repeats that were:
  #retroelements, DNA transposons, Unclassified, simple repeats, or low complexity
  percs <- vapply(c(11,27,39,47,48), function(x) getPerc(rp[x]), numeric(1))
  
  #combine simple and low complexity repeats into a single metric
  #also calculate the percent of genome that was non repetitive
  #and the repeats that didn't fit into one of the named categories
  percs <- c(percs[1:3],percs[4] + percs[5],100-tmask,tmask-sum(percs))
  
  fNames <- c('retro','dna','unclass','simple','non-rep','other')
  
  return(data.frame(labs = factor(fNames,fNames[c(5,2,1,4,3,6)]),
                    percs = percs))
}

#combine results from the three different genomes into a single df for plotting
plotDF <- rbind(percPlot(fileNames[1]), percPlot(fileNames[2]), percPlot(fileNames[3]))
#label the species of origin for all the statistics
plotDF$spec <- rep(c('AEP','105','Olig'),c(6,6,6))

#generate barplots of repeat percentages, split by species
ggplot(plotDF,aes(x=labs,y=percs,fill=labs)) + geom_col() + facet_wrap(.~spec) + theme_bw()
ggsave('repPercBar.pdf',width = 10,height = 5)