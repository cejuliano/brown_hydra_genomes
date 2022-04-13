#! /bin/bash
#SBATCH --job-name=makeRef
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=36G
#SBATCH --error=makeRef.err
#SBATCH --output=makeRef.out
#SBATCH -p med

STAR --version

STAR --runThreadN 24 \
	--runMode genomeGenerate \
	--genomeDir ./ref \
	--genomeFastaFiles aep.genome.cplxmask.fa \
	--genomeSAindexNbases 13
