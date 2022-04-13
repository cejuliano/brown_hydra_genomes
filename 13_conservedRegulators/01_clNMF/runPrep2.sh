#!/bin/bash
#SBATCH -p med
#SBATCH --job-name=prep
#SBATCH -c 2
#SBATCH -t 60-0
#SBATCH --mem=24G
#SBATCH --error=prep2.err
#SBATCH --output=prep2.out

source ~/.bash_profile

conda activate cnmf_env

python ../cNMF/cnmf.py prepare --name cl_fine \
	-c cl.raw.counts.tsv \
	-k 35 36 37 38 39 40 41 42 43 44 45 \
	--n-iter 200 \
	--total-workers 15 \
	--seed 12345 \
	--tpm cl.norm.counts.tsv \
	--genes-file cl.genes.tsv

