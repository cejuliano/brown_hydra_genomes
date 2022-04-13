#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=pasaP
#SBATCH -c 60
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=pasaP.err
#SBATCH --output=pasaP.out

module load singularity

singularity exec -B /home/jacazet/reference/makerAnnotations/aepAnnot/pasa \
	~/pasa.sif ./runAlignment.sh\
