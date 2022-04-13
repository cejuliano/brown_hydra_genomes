#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=purge
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=purge.err
#SBATCH --output=purge.out

purge_haplotigs contigcov -i posSort.10x.bam.gencov -l 17 -m 67 -h 185

purge_haplotigs purge -d -g ../resources/references/canPilTig/tigmint.fa -c coverage_stats.csv -t 24 -b posSort.10x.bam
