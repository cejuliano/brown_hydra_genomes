#!/bin/bash
#SBATCH -p med
#SBATCH --job-name=ACAscores
#SBATCH -c 1
#SBATCH -t 30-0
#SBATCH --mem=20G
#SBATCH --error=ACAscore.err
#SBATCH --output=ACAscore.out

for i in *chroms.genome.aca.hic; do
	echo "$i"
	./aidenlab-3d-dna-cb63403/supp/score-aca.sh "$i"
done


