#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=3D
#SBATCH -c 60
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=3D_%j.err
#SBATCH --output=3D_%j.out

module load bio

cd 3d-dna

./run-asm-pipeline-post-review.sh -r jelly.shrink.rawchrom.review.assembly \
	-g 100 \
	--sort-output \
	../../resources/references/pbj/jelly.shrink.fasta \
	../work/aligned/merged_nodups.txt
