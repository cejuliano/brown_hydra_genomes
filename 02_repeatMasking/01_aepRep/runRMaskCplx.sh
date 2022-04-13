#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=RMask
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=RMask_%j.err
#SBATCH --output=RMask_%j.out

module load singularity

singularity exec -B /home/jacazet/reference/makerAnnotations/aepAnnot/repMask ~/maker-plus_3.01.03.sif RepeatMasker -nolow -norna -pa 24 -lib aep-families.fa aep.final.genome.fa
