#!/bin/bash
#SBATCH --job-name=aggr
#SBATCH -p bigmemm 
#SBATCH -c 1
#SBATCH -t 60-0
#SBATCH --mem=16G
#SBATCH --error=aggr.err
#SBATCH --output=aggr.out

conda deactivate
conda deactivate

cellranger-6.0.2/cellranger aggr \
	--id=clAggr \
	--csv=aggr.csv \
	--normalize=none \
	--nosecondary
