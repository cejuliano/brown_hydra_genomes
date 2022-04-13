#!/bin/bash -l
#SBATCH -p med
#SBATCH --job-name=TD
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=60G
#SBATCH --error=TD.err
#SBATCH --output=TD.out

module load TransDecoder/5.2.0

TransDecoder.LongOrfs -t Hech-trinity.fa

echo "blasting"

~/bin/blastp -query Hech-trinity.fa.transdecoder_dir/longest_orfs.pep \
	-db proteins -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 24 > blastpHeck.outfmt6

TransDecoder.Predict -t Hech-trinity.fa --single_best_only --retain_blastp_hits blastpHeck.outfmt6

