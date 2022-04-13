#!/bin/bash -l
#SBATCH --job-name=clCount1
#SBATCH -p bigmemm
#SBATCH -c 8
#SBATCH -t 60-0
#SBATCH --mem=64G
#SBATCH --error=clCount_1.err
#SBATCH --output=clCount_1.out

cellranger-6.0.2/cellranger count \
	--id=lane1 \
	--fastqs=raw/lane1 \
	--transcriptome=clytiaTran \
	--sample=FT-SA16888,FT-SA16889,FT-SA16890,FT-SA16891 \
	--chemistry=SC3Pv2 \
	--localcores=8 \
	--localmem=64
