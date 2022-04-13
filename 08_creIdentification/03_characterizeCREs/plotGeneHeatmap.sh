#!/bin/bash

plotHeatmap -m geneMatrix.txt \
	-o geneBodyHeatmap.pdf \
	--colorMap magma \
	-max 1.5 1 0.6 0.6 0.6 0.6 0.6 \
	--yMax 1.5 4.5 1.75 1.75 1.75 1.75 1.75
	
	