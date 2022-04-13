#! /bin/bash -l
#SBATCH -p med
#SBATCH --job-name=pasaC
#SBATCH -c 1
#SBATCH -t 60-0
#SBATCH --mem=16G
#SBATCH --error=pasaC.err
#SBATCH --output=pasaC.out

module load singularity

singularity exec -B /home/jacazet/reference/makerAnnotations/aepAnnot/pasa/clPasaRedo \
	~/pasa.sif /usr/local/src/PASApipeline/bin/seqclean clytia.fa
