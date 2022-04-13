#!/bin/sh
#SBATCH --job-name=mmPB
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=60-0
#SBATCH --partition=med
#SBATCH --error=mmPB.err
#SBATCH --output=mmPB.out

minimap2 -ax map-pb -t 22 \
	../resources/references/hicPbj/aepChr.gapfill.pbref.mmi \
	../resources/reads/correctedLR/aepPB.correctedReads.fasta |
	samtools view -b - > pb.bam
