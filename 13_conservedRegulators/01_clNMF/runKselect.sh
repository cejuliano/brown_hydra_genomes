#!/bin/bash
#SBATCH --job-name=kSel
#SBATCH -p bigmemm 
#SBATCH -c 1
#SBATCH -t 60-0
#SBATCH --mem=24G
#SBATCH --error=kSel.err
#SBATCH --output=kSel.out

source ~/.bash_profile

conda activate cnmf_env

python ../cNMF/cnmf.py k_selection_plot --name cl_course 
