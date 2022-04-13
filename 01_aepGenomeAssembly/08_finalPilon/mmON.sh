#!/bin/sh
#SBATCH --job-name=mmON
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=60-0
#SBATCH --partition=med
#SBATCH --error=mmON.err
#SBATCH --output=mmON.out

minimap2 -ax map-ont -t 22 \
	../resources/references/hicPbj/aepChr.gapfill.onref.mmi \
	../resources/reads/correctedLR/NPa*.fasta |
	samtools view -b - > on.bam
