#!/bin/bash
#SBATCH -p med
#SBATCH --job-name=prep
#SBATCH -c 4
#SBATCH -t 60-0
#SBATCH --mem=60G
#SBATCH --error=prep2.err
#SBATCH --output=prep2.out

source ~/.bash_profile

conda activate cnmf_env

python ../cNMF/cnmf.py prepare --name whole_unfilt_fine_broad \
	-c unfilt.whole.raw.counts.tsv \
	-k 20 21 22 23 24 25 26 27 28 29 30 \
	--n-iter 200 \
	--total-workers 15 \
	--seed 12345 \
	--tpm unfilt.whole.norm.counts.tsv \
	--genes-file unfilt.whole.genes.tsv

python ../cNMF/cnmf.py prepare --name whole_unfilt_fine_narrow \
        -c unfilt.whole.raw.counts.tsv \
        -k 50 51 52 53 54 55 56 57 58 59 60 \
        --n-iter 200 \
        --total-workers 15 \
        --seed 12345 \
        --tpm unfilt.whole.norm.counts.tsv \
        --genes-file unfilt.whole.genes.tsv
