#!/bin/bash

for arg in mgBeds/*fa
do
	outName="${arg/_scoredProms.fa/}"
	outName="${outName/mgBeds\//}"
	ame --oc enOutShuf/$outName \
		--bfile clytiaBG.txt \
		$arg shuffledJasparMotifs.meme.txt
done