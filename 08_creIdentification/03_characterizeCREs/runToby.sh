#! /bin/bash -l
#SBATCH --job-name=toby
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=36G
#SBATCH --error=tobias.err
#SBATCH --output=tobias.out
#SBATCH -p med

conda activate tobiasEnv

TOBIAS ATACorrect --bam AEP_MG_final.bam \
                --genome HVAEP1_genome.fa \
                --peaks consensusATAC.bed \
                --outdir corrected \
                --cores 24

TOBIAS FootprintScores --signal corrected/AEP_MG_final_corrected.bw \
                --regions consensusATAC.bed \
                --output AEP_footprints.bw \
                --cores 24
