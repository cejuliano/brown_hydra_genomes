#!/bin/bash

computeMatrix scale-regions -S '../../../Cut&Stuff/CnT/H273_MG.bw' \
	../../alignment_conservation/windows/aepCon.bw \
	../../repeats/repDensity.bw \
	'../../../Cut&Stuff/CnT/H41_MG.bw' \
	'../../../Cut&Stuff/CnT/H43_MG.bw' \
	'../../../Cut&Stuff/ATAC/AEP_MG_final_shift.bw' \
	-R /Volumes/Data/genome/hic/aep16k_domains.bed \
	-o tadMat_hetCon.gz \
	-m 100000 -b 100000 -a 100000 \
	--averageTypeBins median \
	-bs 1000 \
	--missingDataAsZero -p 6
