#!/bin/bash -l
#SBATCH --job-name=clCount2
#SBATCH -p bigmemm
#SBATCH -c 8
#SBATCH -t 60-0
#SBATCH --mem=64G
#SBATCH --error=clCount_2.err
#SBATCH --output=clCount_2.out

cellranger-6.0.2/cellranger count \
	--id=lane2 \
	--fastqs=raw/lane2 \
	--transcriptome=clytiaTran \
	--sample=FT-SA16892,FT-SA16893,FT-SA16894,FT-SA16895 \
	--chemistry=SC3Pv2 \
	--localcores=8 \
	--localmem=64
