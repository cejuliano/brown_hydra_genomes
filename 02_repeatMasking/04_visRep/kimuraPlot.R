library(rstudioapi)
library(ggplot2)
library(RColorBrewer)

setwd(dirname(getActiveDocumentContext()$path))

kimP <- function(pref,size){
  
  #get file name that has repeat landscape
  fname <- paste0('kimura',pref,'.txt')
  
  #import data
  kDat <- read.delim(fname,sep = ' ')
  
  #drop a weird last column thats all NAs
  kDat <- kDat[,-ncol(kDat)]
  
  #restructure the data so that all the columns 
  #get combined into a single row, with an 
  #additional column to indicate the repeat 
  #class
  kDat.plot <- lapply(2:ncol(kDat),function(x){
    newDF <- kDat[,c(1,x)]
    colnames(newDF) <- c('perc','cov')
    newDF$ident <- colnames(kDat)[x]
    return(newDF)
  })
  
  kDat.plot <- do.call(rbind,kDat.plot)
  
  #combine results from repeat subfamilies
  kDat.plot$ident <- gsub('[.].*','',kDat.plot$ident)
  
  kDat.plot <- aggregate(kDat.plot$cov,list(kDat.plot$perc,kDat.plot$ident),sum)
  
  colnames(kDat.plot) <- c('perc','ident','cov')
  
  #calculate percent coverage based on genome size
  kDat.plot$cov <- kDat.plot$cov/size
  
  #drop RNA results (things like small RNAs and whatnot)
  kDat.plot <- kDat.plot[!grepl('RNA',kDat.plot$ident),]
  
  #drop artefact results
  kDat.plot <- kDat.plot[!grepl('ARTEFACT',kDat.plot$ident),]
  
  #define colors used in plot
  colourCount <- length(unique(kDat.plot$ident))
  pal <- colorRampPalette(c('#556270','#4ECDC4','#C7F464','#FF6B6B','#C44D58'))(colourCount)
  set.seed(12345)
  pal <- sample(pal)
  
  #generate stacked bar plot of repeat landscape
  ggplot(kDat.plot,aes(x=perc,y=cov,fill=ident)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) + 
    scale_fill_manual(name="Families", values = pal) +
    xlim(-1,50) +
    labs(x="Kimura Substitution Level (%)", y="Genome Proportion") + theme_bw()
  
  #save plot
  pname <- paste0('repFamKimura',pref,'.pdf')
  ggsave(pname,width=8,height=6)
  
  
  #repeat the plot data formating process
  kDat.plot <- lapply(2:ncol(kDat),function(x){
    newDF <- kDat[,c(1,x)]
    colnames(newDF) <- c('perc','cov')
    newDF$ident <- colnames(kDat)[x]
    return(newDF)
  })
  
  kDat.plot <- do.call(rbind,kDat.plot)
  
  #this time we'll keep most of the subfamily information
  #we'll just drop the subfamily number
  kDat.plot$ident <- gsub('(^[^.]+[.][^.]+)[.].*','\\1',kDat.plot$ident)
  
  #fix some formatting issues with extra dots
  kDat.plot$ident <- gsub('[.]+','.',kDat.plot$ident)
  
  kDat.plot$ident <- gsub('[.]$','',kDat.plot$ident)
  
  #dot to slash for clarity/readability
  kDat.plot$ident <- gsub('[.]','/',kDat.plot$ident)
  
  #collapse data based on subfamilies
  kDat.plot <- aggregate(kDat.plot$cov,list(kDat.plot$perc,kDat.plot$ident),sum)
  
  colnames(kDat.plot) <- c('perc','ident','cov')
  
  #normalize values by total genome size
  kDat.plot$cov <- kDat.plot$cov/size
  
  #drop RNA repeats
  kDat.plot <- kDat.plot[!grepl('RNA',kDat.plot$ident),]
  
  #generate color set for plotting
  colourCount <- length(unique(kDat.plot$ident))
  #pal <- colorRampPalette(c('#556270','#4ECDC4','#C7F464','#FF6B6B','#C44D58'))(colourCount)
  pal <- c("#ea6519","#494edc","#43cd2d","#8732d9","#97c61d","#b94ff1",
           "#4fb83c","#d640d2","#31cb6b","#ed2bb1","#88b930","#8b36bb",
           "#cdb623","#846bf2","#b0b034","#bc65e2","#5bb354","#ce3aa8",
           "#85af42","#5d53bd","#e08d24","#577ff0","#c99f30","#877ce1",
           "#5e8b2f","#c66dd7","#47a45a","#df73d1","#4fc58a","#e3371f",
           "#36c8d4","#e03842","#38b99a","#e23568","#34b9e1","#d05827",
           "#468ae0","#c07a30","#5558a4","#b49c40","#92489f","#6f9554",
           "#da3f87","#4c905e","#c056a0","#677021","#9b87da","#b0ad5f",
           "#4375b7","#bd442f","#66a1e5","#a05624","#cb94d6","#846625",
           "#87629f","#d7935a","#a4426e","#e69976","#e07aa5","#ab6f4b",
           "#ce4858","#9d4f35","#c3586b","#e36f5b","#c76e66")
  set.seed(3821)
  pal <- sample(pal)
  
  #generate subfamily repeat landscape stacked barplot
  ggplot(kDat.plot,aes(x=perc,y=cov,fill=ident)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) + 
    scale_fill_manual(name="Families", values = pal) +
    xlim(-1,50) +
    labs(x="Kimura Substitution Level (%)", y="Genome Proportion") + theme_bw()
  
  pname <- paste0('repSubFamKimura',pref,'.pdf')
  ggsave(pname,width=13,height=6)
}

prefs <- list(c('105',786368896),c('AEP',900878955),c('Olig',1274416349))

lapply(prefs, function(x) kimP(x[1],as.numeric(x[2])))
