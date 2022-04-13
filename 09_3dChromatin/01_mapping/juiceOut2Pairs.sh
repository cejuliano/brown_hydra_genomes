#! /bin/bash -l
#SBATCH -p med
#SBATCH --job-name=con
#SBATCH -c 8
#SBATCH -t 60-0
#SBATCH --mem=16G
#SBATCH --error=con.err
#SBATCH --output=con.out

merged_nodup2pairs.pl -m 29 -s 6 work/aligned/merged_nodups.txt aep.genome aep
