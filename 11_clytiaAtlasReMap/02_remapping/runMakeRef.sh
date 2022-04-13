#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=cRanger
#SBATCH -c 1
#SBATCH -t 60-0
#SBATCH --mem=16G
#SBATCH --error=cRanger_mkref.err
#SBATCH --output=cRanger_mkref.out

cellranger-6.0.2/cellranger mkref --genome=clytiaTran --fasta=clFinal.longestIso.tran.fa --genes=clFinal.transcriptome.gtf
