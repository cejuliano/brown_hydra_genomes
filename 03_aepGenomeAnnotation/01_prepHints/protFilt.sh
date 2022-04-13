#!/bin/bash
#SBATCH --job-name=cdhit
#SBATCH -p bigmemh
#SBATCH --exclusive
#SBATCH --nodes=1
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=cdhit.err
#SBATCH --output=cdhit.out

cd-hit-2d -i proteins.fa -i2 cnido_prot_sequence.fa -o cdhit.out -c 0.95 -M 360000 -T 0 -s2 0.9
