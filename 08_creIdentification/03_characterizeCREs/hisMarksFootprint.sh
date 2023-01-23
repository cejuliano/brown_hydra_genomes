#!/bin/bash

computeMatrix scale-regions \
	-R consensusH41.bed consensusH43.bed consensusH273.bed \
	-S AEP_footprints.bw \
	-o histoneFootprints.txt \
	--outFileNameMatrix histoneFootprints.names.txt \
	--outFileSortedRegions histoneFootprints.regions.txt \
	-m 1000 \
	-b 5000 \
	-a 5000 \
	--missingDataAsZero \
	-p 6
	