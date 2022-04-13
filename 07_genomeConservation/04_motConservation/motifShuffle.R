library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path))

#jaspar motif file

jm <- data.frame(V1=readLines('pooledJasparNR.txt'))

jm <- split(jm,rep(1:(nrow(jm)/5),each=5))

jm.names <- vapply(jm,function(x) x[1,1],"")

jm <- lapply(jm,function(x){
  newDF <- x[-1,]
  newDF <- gsub('.*\\[','',newDF)
  newDF <- gsub('\\].*','',newDF)
  newDF <- gsub(' +',' ',newDF)
  newDF <- gsub('^ | $','',newDF)
  newDF <- strsplit(newDF,split=' ')
  newDF <- do.call(rbind,newDF)
  rownames(newDF) <- c('A','C','G','T')
  return(newDF)
})

lapply(1:length(jm), function(x){

  motName <- names(jm)[x]
  
  motTest <- jm[[x]]
  
  motL <- ncol(motTest)
  
  cutoff <- 5
  
  while(T){
    
    for(i in 1:20){
      motTest <- motTest[,sample(ncol(motTest))]
      
      motTest <- apply(motTest,1,paste,collapse=' ')
      
      motTest <- vapply(motTest, function(x) paste0('[ ',x, ' ] '),'')
      
      motTest <- data.frame(base=c('A','C','T','G'),freq=motTest)
      
      motTest <- apply(motTest,1,paste, collapse= ' ')
      
      motTest <- c(motName,motTest)
      
      writeLines(motTest,'shufMot.txt')
      
      system('/Users/Jcazet/meme/libexec/meme-5.4.1/jaspar2meme -bundle shufMot.txt > shufMot.meme.txt')
      
      system('/Users/Jcazet/meme/bin/tomtom -thresh 1 -text shufMot.meme.txt pooledJasparNR.meme.txt > shufMotScore.txt')
      
      matchRes <- read.delim('shufMotScore.txt')
      
      matchRes <- matchRes[complete.cases(matchRes),]
      
      if(min(matchRes$E.value) >= cutoff){
        break
      }
      
      print('Motif too similar. Retrying.')
      
      motTest <- jm[[x]]
      
    }
    
    if(min(matchRes$E.value) >= cutoff){
      break
    }
    
    print('Lowering cutoff to:')
    
    cutoff <- cutoff/2
    
    print(cutoff)
    
  }
  
  system('cat shufMot.txt >> shuffledJasparMotifs.txt')
  
})

