#! /bin/bash
#SBATCH -p med 
#SBATCH --job-name=align
#SBATCH -t 60-0
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --error=align.err
#SBATCH --output=align.out

STAR --version

STAR --runThreadN 20 \
	--genomeDir ./ref \
	--readFilesIn ../reads/combined_R1.fq.gz ../reads/combined_R2.fq.gz \
	--readFilesCommand gunzip -c \
	--outFileNamePrefix ./out/aep \
	--outSAMprimaryFlag AllBestScore \
	--outSAMtype BAM SortedByCoordinate \
	--twopassMode Basic \
	--outFilterScoreMinOverLread 0.3 \
	--outFilterMatchNminOverLread 0.3 \
	--limitBAMsortRAM 12316579964 
