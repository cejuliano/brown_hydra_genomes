#! /bin/bash

specUse="$1"

specDir="$2"

awk  -v OFMT='%f' -F '\t' '{print $1, int($2/2), int($2/2)+10000}' "$specUse".chroms.genome > "$specUse".pseudochroms.bed

./aidenlab-3d-dna-cb63403/supp/build-aca-hic.sh \
	"$specUse".chroms.genome \
	"$specUse".pseudochroms.bed \
	"$specDir"/aligned/merged_nodups.txt

