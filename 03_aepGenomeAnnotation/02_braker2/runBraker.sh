#!/bin/bash
#SBATCH --job-name=braker
#SBATCH -p med
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=braker.err
#SBATCH --output=braker.out

module load singularity

singularity exec -H "/home/jacazet" -B /home/jacazet/reference/makerAnnotations/aepAnnot/maker_braker/braker ~/braker2_2.1.5.sif ./brakerScript.sh
