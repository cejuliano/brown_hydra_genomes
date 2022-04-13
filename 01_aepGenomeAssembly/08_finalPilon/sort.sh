#! /bin/bash
#SBATCH --job-name=sort
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=60-0
#SBATCH --partition=med
#SBATCH --error=sort.err
#SBATCH --output=sort.out

samtools sort -@ 24 -o pb.sort.bam pb.bam
samtools sort -@ 24 -o on.sort.bam on.bam

samtools index pb.sort.bam
samtools index on.sort.bam
