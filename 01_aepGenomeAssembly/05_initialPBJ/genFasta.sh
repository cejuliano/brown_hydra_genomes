#! /bin/bash
#SBATCH -p bigmemm
#SBATCH --job-name=bm2fa
#SBATCH -c 2
#SBATCH -t 60-0
#SBATCH --mem=16G
#SBATCH --error=b2f_%j.err
#SBATCH --output=b2f_%j.out

samtools bam2fq pb.subreads.bam | gzip > pb.subreads.fq.gz
