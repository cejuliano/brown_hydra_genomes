#!/bin/bash
#SBATCH --job-name=makeDGE
#SBATCH -c 1 
#SBATCH -t 60-0
#SBATCH --mem=8G
#SBATCH -p med
#SBATCH --error=makeDGE.err
#SBATCH --output=makeDGE.out

for arg in */final.bam
do
	echo $arg
	
	outDir="${arg/\/final.bam/}"
	
	echo $outDir
	
	cellCount=$(cat $outDir/*cellCount.txt)
	
	echo $cellCount
	
	/group/julianolab/analyses/dropseq/Drop-seq_tools-2.4.0/DigitalExpression \
		I=$arg \
		O="$outDir"/"$outDir".dge.txt.gz \
		NUM_CORE_BARCODES="$cellCount"
done
