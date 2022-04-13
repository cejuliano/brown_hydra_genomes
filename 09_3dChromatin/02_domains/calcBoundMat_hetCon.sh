#!/bin/bash

computeMatrix scale-regions -S '../../../Cut&Stuff/CnT/H273_MG.bw' \
	../../alignment_conservation/windows/aepCon.bw \
	-R /Volumes/Data/genome/hic/aep16k_boundaries.bed \
	-o boundMat_hetCon.gz \
	-m 30000 -b 100000 -a 100000 \
	--missingDataAsZero -p 6