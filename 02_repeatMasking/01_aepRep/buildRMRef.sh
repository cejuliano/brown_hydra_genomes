#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=RMDB
#SBATCH -c 8
#SBATCH -t 60-0
#SBATCH --mem=16G
#SBATCH --error=RMDB.err
#SBATCH --output=RMDB.out

../../RepeatModeler-2.0.1/BuildDatabase -name aep aep.final.genome.fa
