#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=canuCor
#SBATCH -c 8
#SBATCH -t 60-0
#SBATCH --mem=60G
#SBATCH --error=canuCor_%j.err
#SBATCH --output=canuCor_%j.out

canu -correct -s canuSpecPB.txt -p aepPB -d aepPB genomeSize=1.25g -pacbio-raw pb.subreads.fq.gz
