#!/bin/bash

gsed 's/\(mRNA\.RE[0-9]\+\)\t/\1_/g' \
	Results_Sep15_1/Phylogenetic_Hierarchical_Orthogroups/"$1".tsv \
	> Results_Sep15_1/Phylogenetic_Hierarchical_Orthogroups/"$1"Mod.tsv