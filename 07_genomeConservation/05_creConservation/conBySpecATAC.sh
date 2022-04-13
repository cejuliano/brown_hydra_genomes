#!/bin/bash

computeMatrix scale-regions \
	-R ../ATAC/consensusAEP.bed \
	-S 105.cactus.bw olig.cactus.bw virid.cactus.bw clytia.cactus.bw \
	-o geneMatrix.ATAC.specCon.txt \
	--outFileNameMatrix geneMatrix.ATAC.specCon.names.txt \
	--outFileSortedRegions geneMatrix.ATAC.specCon.regions.txt \
	-m 1000 \
	-b 5000 \
	-a 5000 \
	--missingDataAsZero \
	-p 4