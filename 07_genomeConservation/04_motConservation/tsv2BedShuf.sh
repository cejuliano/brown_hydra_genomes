#!/bin/bash -l
#SBATCH --job-name=makeBed
#SBATCH -p med
#SBATCH -c 1
#SBATCH -t 7-0
#SBATCH --mem=16G
#SBATCH --error=bed.err
#SBATCH --output=bed.out

for arg in *Mots.tsv
do
	newName=${arg/tsv/bed}
	awk 'BEGIN {OFS="\t"}; {if (NR!=1) {print($3,$4-1,$5,$1,$8,$6)}}' $arg > $newName
done

source ~/reference/alignments/venv/bin/activate

halLiftover evolverHydra.hal olig olig.shufMots.bed aep oligShufMotifsAEP.bed
halLiftover evolverHydra.hal 105 105.shufMots.bed aep 105ShufMotifsAEP.bed
halLiftover evolverHydra.hal virid virid.shufMots.bed aep viridShufMotifsAEP.bed
