#!/bin/bash
#SBATCH --job-name=consensus
#SBATCH -p bigmemm
#SBATCH -c 1
#SBATCH -t 60-0
#SBATCH --mem=24G
#SBATCH --error=consensus.err
#SBATCH --output=consensus.out

source ~/.bash_profile

conda activate cnmf_env

python ../cNMF/cnmf.py consensus --name cl_fine --local-density-threshold 0.13 --components 37 --show-clustering
