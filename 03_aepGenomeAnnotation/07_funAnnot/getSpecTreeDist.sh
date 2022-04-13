#! /bin/bash
#SBATCH -p bigmemh 
#SBATCH --job-name=orthoF
#SBATCH --exclusive
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=orthoF.err
#SBATCH --output=orthoF.out

cd OrthoFinder

./orthofinder -S diamond_ultra_sens \
        -M msa \
        -I 1.3 \
        -t 60 \
        -b ../orthoOut/Results_Sep15_1/WorkingDirectory/
