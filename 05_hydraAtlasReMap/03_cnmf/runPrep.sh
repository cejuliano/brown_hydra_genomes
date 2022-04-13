#!/bin/bash
#SBATCH -p med
#SBATCH --job-name=prep
#SBATCH -c 4
#SBATCH -t 60-0
#SBATCH --mem=60G
#SBATCH --error=prep.err
#SBATCH --output=prep.out

source ~/.bash_profile

conda activate cnmf_env

python ../cNMF/cnmf.py prepare --name whole_unfilt_course \
	-c unfilt.whole.raw.counts.tsv \
	-k 15 20 25 30 35 40 45 50 55 60 65 70 75 80 85 90 \
	--n-iter 200 \
	--total-workers 15 \
	--seed 12345 \
	--tpm unfilt.whole.norm.counts.tsv \
	--genes-file unfilt.whole.genes.tsv

