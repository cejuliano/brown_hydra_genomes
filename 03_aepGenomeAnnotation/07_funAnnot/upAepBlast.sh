#!/bin/bash

diamond blastp -p 6 -d up \
	-o uniprotBlast.txt -f 6 \
	-q HVAEP1.prot.longestIso.fa -k 1 \
	-M 12
	
gsed -i 's/sp|//g' uniprotBlast.txt