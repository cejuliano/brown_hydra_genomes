#!/bin/bash

computeMatrix scale-regions \
	-R consensusH273.bed \
	-S 105.cactus.bw olig.cactus.bw virid.cactus.bw clytia.cactus.bw \
	-o geneMatrix.h273.specCon.txt \
	--outFileNameMatrix geneMatrix.h273.specCon.names.txt \
	--outFileSortedRegions geneMatrix.h273.specCon.regions.txt \
	-m 1000 \
	-b 5000 \
	-a 5000 \
	--missingDataAsZero \
	-p 4