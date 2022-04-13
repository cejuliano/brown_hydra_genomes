#!/bin/bash -l
#SBATCH --job-name=fimo
#SBATCH -p med
#SBATCH -c 1
#SBATCH -t 7-0
#SBATCH --array=0-3
#SBATCH --mem=16G
#SBATCH --error=fimo_%a.err
#SBATCH --output=fimo_%a.out

conda activate meme

array=(aep 105 olig virid)

specUse=${array[$SLURM_ARRAY_TASK_ID]}

echo $specUse

fimo -bfile genome.markov.txt \
	--max-strand \
	--skip-matched-sequence \
	pooledJasparNR.meme.txt seqs/$specUse.hm.fa > $specUse.mots.tsv
