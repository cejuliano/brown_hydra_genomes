#!/bin/bash -l
#SBATCH -p med
#SBATCH --job-name=RNA_Pipeline
#SBATCH -c 12
#SBATCH -t 7-0
#SBATCH --array=0-8
#SBATCH --mem=24G
#SBATCH --error=RNA_Pipeline_%a.err
#SBATCH --output=RNA_Pipeline_%a.out

module load bowtie2
module load fastqc/0.11.4
module load rsem

array=(F1 F2 F3 M1 M2 M3 NS1 NS2 NS3)

resources/RNA_Mapping_Pipeline.sh ${array[$SLURM_ARRAY_TASK_ID]}
