library(Gviz)
library(rstudioapi)
library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(grDevices)

setwd(dirname(getActiveDocumentContext()$path))

genomeInfo <- read.table("../../genome/aep.final.genome.fa.fai", stringsAsFactors = F)[,1:2]

chrominfo <- data.frame(chrom=genomeInfo$V1, length=genomeInfo$V2, is_circular=rep(FALSE,length(genomeInfo$V1)))
txDb.mi <- makeTxDbFromGFF(file="HVAEP1.GeneModels.pfix.gff3", format="gff", dataSource="Gene Models", organism ="Hydra vulgaris", chrominfo=chrominfo)

options(ucscChromosomeNames=F)

pal <- colorRampPalette(c('#3FB8AF','#7FC7AF','#DAD8A7','#FF9E9D','#FF3D7F'),bias=1)

palCols <- pal(8)

AtacPlot <- function(chr, height, left, right, buffer=0,fpath = "plot.pdf", minCon=0) {
  chrID <- chr
  
  #This specifies the maximum value on the y-axis
  z <- height
  
  #chromosome map
  gtrack <- GenomeAxisTrack(name = chrID,add53=T, add35 = T, fontsize = 13, fontcolor.title = "black", 
                            fontsize.title = 13, showTitle = F, rotation.title = 0, grid = T,
                            cex = 0.6, labelPos = "below")
  #gene models
  grtrack <- GeneRegionTrack(txDb.mi, chromosome=chrID, name="Genes", transcriptAnnotation = "gene", 
                             col = "black", fontsize.group = 13, fontcolor.group = "black", fill = "black", 
                             fontsize=25, rotation.title = 0, background.title = "white", col.line = "black",
                             just.group = "below", collapseTranscripts = "longest")
  
  reptrack <- AnnotationTrack(range = '../../genome/repeats/bothMaskFull.out.sorted.gff', name = 'repeats',
                              chromosome = chrID,shape = "box",col = "black",background.title = "white",
                              fill = "black",cex=0)
  
  magCon <- DataTrack(range = "105.rolling100.cactus.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c(palCols[5],palCols[5]), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(minCon,1), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  oligCon <- DataTrack(range = "olig.rolling100.cactus.bw", 
                      type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c(palCols[6],palCols[6]), 
                      chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(minCon,1), 
                      background.title = "white", fontcolor.title = "black", col.axis = "black", 
                      span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  viridCon <- DataTrack(range = "virid.rolling100.cactus.bw", 
                       type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c(palCols[7],palCols[7]), 
                       chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(minCon,1), 
                       background.title = "white", fontcolor.title = "black", col.axis = "black", 
                       span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  clytiaCon <- DataTrack(range = "clytia.rolling100.cactus.bw", 
                       type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c(palCols[8],palCols[8]), 
                       chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(minCon,1), 
                       background.title = "white", fontcolor.title = "black", col.axis = "black", 
                       span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  ATAC <- DataTrack(range = "../ATAC/AEP_MG_final_shift.bw", 
                       type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c(palCols[1],palCols[1]), 
                       chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z[1]), 
                       background.title = "white", fontcolor.title = "black", col.axis = "black", 
                       span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  H41 <- DataTrack(range = "H41_MG.bw", 
                    type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c(palCols[2],palCols[2]), 
                    chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z[2]), 
                    background.title = "white", fontcolor.title = "black", col.axis = "black", 
                    span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  H43 <- DataTrack(range = "H43_MG.bw", 
                   type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c(palCols[3],palCols[3]), 
                   chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z[3]), 
                   background.title = "white", fontcolor.title = "black", col.axis = "black", 
                   span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  H273 <- DataTrack(range = "H273_MG.bw", 
                   type = "polygon", col = "black", lwd.mountain = 0, fill.mountain = c(palCols[4],palCols[4]), 
                   chromosome = chrID, name = "Control", fontsize = 13, showAxis = T, ylim=c(0,z[4]), 
                   background.title = "white", fontcolor.title = "black", col.axis = "black", 
                   span = 0.5, lwd = 0, showTitle = F, cex.axis = 0.6)
  
  
  pdf(file = fpath, width=7, height=7)
  plotTracks(list(ATAC,H41,H43,H273,
                  grtrack,gtrack,reptrack,magCon,
                  oligCon,viridCon,clytiaCon), 
             from=(left - buffer), to=(right + buffer), title.width = 0.7, margin = 0, sizes = c(rep(5,4),3,3,1,rep(5,4)))
  dev.off()
}

AtacPlot("chr-3",c(3,7,5,2),40525583,40580186,0, "brachyurySample.pdf")
AtacPlot("chr-4",c(5,16,3,4),17489011,17570761,0, "bmp5-8cSample.pdf")

#conservation plots
AtacPlot("chr-6",c(3,7,5,2),9544624,9546772,0, "wnt3ConSample.pdf",minCon = 0.5)
AtacPlot("chr-6",c(5,16,3,4),24222683,24225675,0, "sp5ConSample.pdf",minCon=0.5)
