#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=sort
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=36G
#SBATCH --error=sort.err
#SBATCH --output=sort.out

samtools sort -@ 24 -n possorted_bam.bam -o nameSorted.bam
