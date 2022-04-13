#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=canuCor
#SBATCH -c 8
#SBATCH -t 60-0
#SBATCH --mem=60G
#SBATCH --error=canuCor_%j.err
#SBATCH --output=canuCor_%j.out

canu -correct -s canuSpecNano.txt -p aepNP -d aepNP genomeSize=1.25g -nanopore-raw nanoReads.fastq.gz
