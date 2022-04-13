#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=Pilon
#SBATCH -c 45
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=pilon_%j.err
#SBATCH --output=pilon_%j.out

module load pilon

chrs=( chr-1 \
	chr-2 \
	chr-3 \
	chr-4 \
	chr-5 \
	chr-6 \
	chr-7 \
	chr-8 \
	chr-9 \
	chr-10 \
	chr-11 \
	chr-12 \
	chr-13 \
	chr-14 \
	chr-15 )

for arg in "${chrs[@]}"
do

	echo "$arg"

	samtools faidx ../resources/references/hicPbj/aepChr.gapfill.fa $arg > subGenome.fa

	java -Xmx460G -jar /share/apps/pilon-1.23/pilon-1.23.jar \
		--genome subGenome.fa --bam possorted_bam.bam \
		--nanopore on.sort.bam --pacbio pb.sort.bam \
		--output $arg --outdir pilOut --threads 45

done

rm subGenome.fa
