#!/bin/bash

for arg in mgBeds/*fa
do
	outName="${arg/_scoredProms.fa/}"
	outName="${outName/mgBeds\//}"
	ame --oc enOut/$outName \
		--bfile clytiaBG.txt \
		$arg pooledJasparNR.meme.txt
done