#! /bin/bash -l
#SBATCH -p med
#SBATCH --job-name=pasaC
#SBATCH -c 4
#SBATCH -t 60-0
#SBATCH --mem=16G
#SBATCH --error=pasaC.err
#SBATCH --output=pasaC.out

module load singularity

singularity exec -B /home/jacazet/reference/makerAnnotations/aepAnnot/pasa \
	~/pasa.sif /usr/local/src/PASApipeline/bin/seqclean Trinity-GG.fasta
