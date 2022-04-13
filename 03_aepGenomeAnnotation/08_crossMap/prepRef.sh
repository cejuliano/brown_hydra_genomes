#!/bin/bash
#SBATCH -p med
#SBATCH --job-name=prepR
#SBATCH -c 12
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=pRef.err
#SBATCH --output=pRef.out

module load rsem

rsem-prepare-reference --gtf hydra.augustus.gtf \
	--star -p 12 --transcript-to-gene-map dv.g2tMap.txt \
	Hm105_Dovetail_Assembly_1.0.fasta dvStar
