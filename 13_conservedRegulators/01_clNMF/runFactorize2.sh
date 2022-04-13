#!/bin/bash
#SBATCH --job-name=factorize
#SBATCH -p med
#SBATCH -c 2
#SBATCH -t 60-0
#SBATCH --array=0-14
#SBATCH --mem=8G
#SBATCH --error=factorizeB_%a.err
#SBATCH --output=factorizeB_%a.out

source ~/.bash_profile

conda activate cnmf_env

python ../cNMF/cnmf.py factorize --name cl_fine --worker-index $SLURM_ARRAY_TASK_ID
