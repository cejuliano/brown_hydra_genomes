#!/bin/bash
#SBATCH -p med
#SBATCH --job-name=blastp
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=blastp.err
#SBATCH --output=blastp.out

diamond blastp --query gmCandidates.fa \
		--db allPrimProts --sensitive \
		--outfmt 6 --evalue 1e-5 -p 24 > blastpPrimProt.outfmt6
