#!/bin/bash
#SBATCH --job-name=CnT_Map
#SBATCH -c 24
#SBATCH -p med
#SBATCH -t 60-0
#SBATCH --array=0-11
#SBATCH --mem=0
#SBATCH --error=CnT_Map_%a.err
#SBATCH --output=CnT_Map_%a.out

module load fastqc/0.11.4
module load java/1.8

array=( H273_1 H273_2 H273_3 \
H41_1 H41_2 H41_3 \
H43_1 H43_2 H43_3 \
IGG_1 IGG_2 IGG_3 )

./resources/mapPipe.sh ${array[$SLURM_ARRAY_TASK_ID]}
