#!/bin/bash

computeMatrix scale-regions \
	-R HVAEP1.GeneModels.gtf \
	-S aepCon.bw NS_RNA.bw ../ATAC/AEP_MG_final_shift.bw IGG_MG.bw H43_MG.bw H41_MG.bw H273_MG.bw \
	-o geneMatrix.txt \
	--outFileNameMatrix geneMatrix.names.txt \
	--outFileSortedRegions geneMatrix.regions.txt \
	-m 500 \
	-b 5000 \
	-a 5000 \
	--missingDataAsZero \
	-p 6 \
	--metagene
