#!/bin/bash
#SBATCH --job-name=readCounts
#SBATCH -c 1 
#SBATCH -t 60-0
#SBATCH --mem=8G
#SBATCH -p med
#SBATCH --error=readCounts.err
#SBATCH --output=readCounts.out

for arg in */final.bam
do
	echo $arg
	
	outDir="${arg/\/final.bam/}"
	
	echo $outDir
	
	/group/julianolab/analyses/dropseq/Drop-seq_tools-2.4.0/BamTagHistogram \
		I=$arg \
		O="$outDir"/out_cell_readcounts.txt.gz \
		TAG=XC
done
