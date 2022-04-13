#!/bin/bash
#SBATCH --job-name=ipr
#SBATCH -p med 
#SBATCH -c 8
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=ipr.err
#SBATCH --output=ipr.out

cd interproscan-5.51-85.0 

./interproscan.sh -d final -cpu 8 -dp -f TSV, GFF3 -goterms -i ../HVAEP1.prot.longestIso.fa -iprlookup -pa
