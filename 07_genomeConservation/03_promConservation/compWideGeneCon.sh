#!/bin/bash

computeMatrix scale-regions \
	-R HVAEP1.GeneModels.rmod.gtf \
	-S aepCon.bw \
	-o broadConGeneMatrix.txt \
	--outFileNameMatrix broadConGeneMatrix.names.txt \
	--outFileSortedRegions broadConGeneMatrix.regions.txt \
	-m 750 \
	-b 10000 \
	-a 10000 \
	--missingDataAsZero \
	-p 6 \
	--metagene
	