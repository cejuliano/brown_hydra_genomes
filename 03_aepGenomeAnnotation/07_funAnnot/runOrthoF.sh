#! /bin/bash
#SBATCH -p med 
#SBATCH --job-name=orthoF
#SBATCH --exclusive
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=orthoF.err
#SBATCH --output=orthoF.out

cd OrthoFinder

./orthofinder -S diamond_ultra_sens -M msa -s ../orthoTree.txt -o ../orthoOut -I 1.3 -t 24 -f ../primary_transcripts
