#! /bin/bash -l
#SBATCH --job-name=exo
#SBATCH -p med
#SBATCH -c 4
#SBATCH -t 60-0
#SBATCH --array=0-23
#SBATCH --mem=8G
#SBATCH --error=exo_%a.err
#SBATCH --output=exo_%a.out

source ~/.bash_profile

conda activate agatEnv

array=(query.fa.split/*)

./gbMap.sh ${array[$SLURM_ARRAY_TASK_ID]}
