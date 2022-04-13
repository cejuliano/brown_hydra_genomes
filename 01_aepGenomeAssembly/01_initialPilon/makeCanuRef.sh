#! /bin/bash
#SBATCH -p bigmemm
#SBATCH --job-name=lRanger
#SBATCH -c 4
#SBATCH -t 60-0
#SBATCH --mem=16G
#SBATCH --error=lRanger_mkref_%j.err
#SBATCH --output=lRanger_mkref_%j.out

../../longranger-2.2.2/longranger mkref hydra_aep.canu.contigs.fasta
