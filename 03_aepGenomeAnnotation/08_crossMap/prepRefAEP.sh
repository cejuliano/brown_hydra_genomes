#!/bin/bash
#SBATCH -p med
#SBATCH --job-name=prepR
#SBATCH -c 12
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=pRef.err
#SBATCH --output=pRef.out

module load rsem

rsem-prepare-reference --gtf HVAEP1.GeneModels.gtf \
	--star -p 12 --transcript-to-gene-map HVAEP1.t2gMap.txt \
	HVAEP1_genome.fa aepGenome
