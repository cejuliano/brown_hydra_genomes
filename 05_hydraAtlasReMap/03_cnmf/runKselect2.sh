#!/bin/bash
#SBATCH --job-name=combine
#SBATCH -p bigmemm 
#SBATCH -c 1
#SBATCH -t 60-0
#SBATCH --mem=24G
#SBATCH --error=combine.err
#SBATCH --output=combine.out

source ~/.bash_profile

conda activate cnmf_env

python ../cNMF/cnmf.py k_selection_plot --name whole_unfilt_fine_narrow 
