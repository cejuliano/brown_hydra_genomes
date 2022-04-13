#!/bin/bash -l
#SBATCH --job-name=TD
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=36G
#SBATCH --error=TD.err
#SBATCH --output=TD.out

module load TransDecoder/5.2.0

TransDecoder.LongOrfs -t aepLRv2.fasta

echo "blasting"

~/bin/blastp -query aepLRv2.fasta.transdecoder_dir/longest_orfs.pep \
	-db ~/blastdb/nr -max_target_seqs 1 \
    -outfmt 6 -evalue 1e-5 -num_threads 24 > blastp.outfmt6

TransDecoder.Predict -t aepLRv2.fasta --single_best_only --retain_blastp_hits blastp.outfmt6

