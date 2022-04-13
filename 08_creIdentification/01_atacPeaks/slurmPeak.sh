#! /bin/bash -l
#SBATCH --job-name=ATAC_Peaks
#SBATCH -p med
#SBATCH -c 16
#SBATCH -t 7-0
#SBATCH --mem=0
#SBATCH --error=ATAC_Peaks.err
#SBATCH --output=ATAC_Peaks.out

conda activate deepEnv

resources/ATAC_Peak_Pipeline.sh AEP
