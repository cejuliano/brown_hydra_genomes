#!/bin/bash -l
#SBATCH --job-name=bg2bw
#SBATCH -p bigmemm 
#SBATCH -c 1
#SBATCH -t 7-0
#SBATCH --mem=360G
#SBATCH --error=bg2bw.err
#SBATCH --output=bg2bw.out

cat *rolling100* > aepCon100.bedgraph

bedSort aepCon100.bedgraph aepCon100.sort.bedgraph

bedGraphToBigWig aepCon100.sort.bedgraph aep.genome aepCon100.bw
