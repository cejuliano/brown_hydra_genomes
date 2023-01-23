#!/bin/bash

plotHeatmap -m tadMat_hetCon.gz -o tadHeat_hetCon.pdf \
	--colorList "white,darkblue" \
	--heatmapHeight 5 \
	--yMax 0.1 0.4 0.9 0.07 0.07 0.08 \
	--yMin 0 0.1 0.7 0 0 0.04 \
	--zMax 1 1 1 0.08 0.2 0.2 
	
	