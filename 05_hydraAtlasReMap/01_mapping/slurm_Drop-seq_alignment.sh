#!/bin/bash
#SBATCH --job-name=dsalign
#SBATCH -c 16 
#SBATCH -t 4-0
#SBATCH --mem=30G
#SBATCH --array=0-15
#SBATCH --error=logs/ds.%A_%a.err
#SBATCH --output=logs/ds.%A_%a.out

array=(D01-D1_S1 D01-P2_S4 \
D02-CO_S2 D02-P1_S3 D02-PB_S4 \
D03-KI_S1 D03-MA_S2 D03-FM_S3 \
D06-FM_S1 D06-MA_S3 D06-KI_S4 \
D11-PO_S1 D11-BU_S2 D12-N1_S1 D12-N2_S2 \
D12-UN_S0) \

file_prefix=${array[$SLURM_ARRAY_TASK_ID]}
echo "$file_prefix"

module load star
bash Drop-seq_alignment.sh \
	-g STAR \
	-r HVAEP1.tran.final.mito.ds.fasta.gz \
	-d /group/julianolab/analyses/dropseq/Drop-seq_tools-2.4.0 \
	-o out/${file_prefix} \
	-t temp/${file_prefix} \
	-k \
	bam/${array[$SLURM_ARRAY_TASK_ID]}.1.bam

