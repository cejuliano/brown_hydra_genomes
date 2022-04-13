#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=sub
#SBATCH -c 2
#SBATCH -t 60-0
#SBATCH --mem=60G
#SBATCH --error=sub.err
#SBATCH --output=sub.out

zcat aepNP/aepNP.correctedReads.fasta.gz | split -l 800000 - NP
