#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=mh
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=mh.err
#SBATCH --output=mh.out

purge_haplotigs readhist -b posSort.10x.bam -g ../resources/references/canPilTig/tigmint.fa -t 24
