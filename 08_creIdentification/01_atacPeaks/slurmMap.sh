#!/bin/bash -l
#SBATCH --job-name=ATAC_Pipeline
#SBATCH -p med
#SBATCH -c 16
#SBATCH -t 7-0
#SBATCH --array=0-2
#SBATCH --mem=0
#SBATCH --error=ATAC_Pipeline_%a.err
#SBATCH --output=ATAC_Pipeline_%a.out

module load bowtie2
module load samtools
module load fastqc/0.11.4
module load bedtools
module load java/1.8

array=(AEP1 AEP2 AEP3)

./resources/ATAC_Mapping_Pipeline.sh ${array[$SLURM_ARRAY_TASK_ID]}

