#!/bin/bash

hicFindTADs -m hic_corrected16k.cool --outPrefix aep16k --correctForMultipleTesting fdr \
	--minDepth 48000 --maxDepth 160000 --step 16000 -p 4 --thresholdComparisons 0.05