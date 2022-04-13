#!/bin/bash -l
#SBATCH -p med
#SBATCH --job-name=TD
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=60G
#SBATCH --error=TD.err
#SBATCH --output=TD.out

module load TransDecoder/5.2.0

for arg in *cdhit.fa
do
	echo "$arg"

	TransDecoder.LongOrfs -t "$arg"

	diamond blastp --query "$arg".transdecoder_dir/longest_orfs.pep \
		--db ../proteins --max-target-seqs 1 \
		--outfmt 6 --evalue 1e-5 -p 24 --sensitive > blastp.outfmt6

	TransDecoder.Predict -t "$arg" --single_best_only --retain_blastp_hits blastp.outfmt6
done


