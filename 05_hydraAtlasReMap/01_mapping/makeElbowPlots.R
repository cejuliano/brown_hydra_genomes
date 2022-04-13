#!/usr/bin/env Rscript


#distance functions
distancePointLine <- function(x, y, slope, intercept) {
  ## x, y is the point to test.
  ## slope, intercept is the line to check distance.
  ##
  ## Returns distance from the line.
  ##
  ## Returns 9999 on 0 denominator conditions.
  x1 <- x-10
  x2 <- x+10
  y1 <- x1*slope+intercept
  y2 <- x2*slope+intercept
  distancePointSegment(x,y, x1,y1, x2,y2)
}

distancePointSegment <- function(px, py, x1, y1, x2, y2) {
  ## px,py is the point to test.
  ## x1,y1,x2,y2 is the line to check distance.
  ##
  ## Returns distance from the line, or if the intersecting point on the line nearest
  ## the point tested is outside the endpoints of the line, the distance to the
  ## nearest endpoint.
  ##
  ## Returns 9999 on 0 denominator conditions.
  lineMagnitude <- function(x1, y1, x2, y2) sqrt((x2-x1)^2+(y2-y1)^2)
  ans <- NULL
  ix <- iy <- 0   # intersecting point
  lineMag <- lineMagnitude(x1, y1, x2, y2)
  if( lineMag < 0.00000001) {
    warning("short segment")
    return(9999)
  }
  u <- (((px - x1) * (x2 - x1)) + ((py - y1) * (y2 - y1)))
  u <- u / (lineMag * lineMag)
  if((u < 0.00001) || (u > 1)) {
    ## closest point does not fall within the line segment, take the shorter distance
    ## to an endpoint
    ix <- lineMagnitude(px, py, x1, y1)
    iy <- lineMagnitude(px, py, x2, y2)
    if(ix > iy)  ans <- iy
    else ans <- ix
  } else {
    ## Intersecting point is on the line, use the formula
    ix <- x1 + u * (x2 - x1)
    iy <- y1 + u * (y2 - y1)
    ans <- lineMagnitude(px, py, ix, iy)
  }
  ans
}

distancePointLineTest <- function() {
  if(abs(distancePointSegment(  5,   5,  10, 10, 20, 20) - 7.07106781186548)>.0001)
    stop("error 1")
  if(abs(distancePointSegment( 15,  15,  10, 10, 20, 20) - 0)>.0001)
    stop("error 2")
  if(abs(distancePointSegment( 15,  15,  20, 10, 20, 20) - 5)>.0001)
    stop("error 3")
  if(abs(distancePointSegment(  0,  15,  20, 10, 20, 20) - 20)>.0001)
    stop("error 4")
  if(abs(distancePointSegment(  0,  25,  20, 10, 20, 20) - 20.6155281280883)>.0001)
    stop("error 5")
  if(abs(distancePointSegment(-13, -25, -50, 10, 20, 20) - 39.8808224589213)>.0001)
    stop("error 6")
  if(abs(distancePointSegment(  0,   3,   0, -4,  5,  0) - 5.466082)>.0001)
    stop("error 7")
  if(abs(distancePointSegment(  0,   9,   0, -4,  0, 15) - 0)>.0001)
    stop("error 8")
  if(abs(distancePointSegment(  0,   0,   0, -2,  2,  0)^2 - 2)>.0001)
    stop("error 9")
  return(TRUE)
}

inF <- list.files(pattern='out_cell_readcounts.txt.gz',full.names = T, recursive = T)

inF <- inF[!is.na(inF)]

for(i in 1:length(inF)) {
  a=read.table(inF[i], header=F, stringsAsFactors=F)
  x=cumsum(a$V1)
  x=x/max(x)
  
  #get slope
  checkVal <- min(which(x > 0.85))

  inSlope <- x[checkVal]/checkVal
  
  testDF <- data.frame(x = 1:checkVal, y = x[1:checkVal])
  
  testDF$res <- 0
  
  for(j in 1:nrow(testDF)){
    testDF$res[j] <-  distancePointLine(testDF[j,'x'],testDF[j,'y'],inSlope,0)
  }
  
  cutoff <- which.max(testDF$res)

  dirName=gsub('[.]/(.*)/.*','\\1',inF[i])
  
  newName=paste0(dirName,'/',dirName,'.pdf')
  
  pdf(newName, width = 6, height = 6)
  
  plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads",xlim=c(1,checkVal*2))
  
  abline(v=cutoff)
  
  text(checkVal,0.5,paste0('estimated # Cells = ',cutoff))
  
  dev.off()
  
  cellCountName <- gsub('.pdf','_cellCount.txt',newName)
  write.table(cutoff,file=cellCountName,col.names = F, row.names = F, quote = F)
}
