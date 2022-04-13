#! /bin/bash -l
#SBATCH -p med
#SBATCH --job-name=pbj
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=pbj_%a.err
#SBATCH --output=pbj_%a.out
#SBATCH --array=1-10

conda deactivate
conda deactivate

source /home/jacazet/PBSuite_15.8.24/setup.sh

Jelly.py mapping --debug mCon$SLURM_ARRAY_TASK_ID.xml
