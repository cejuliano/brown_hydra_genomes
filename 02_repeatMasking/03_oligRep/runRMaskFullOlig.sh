#! /bin/bash
#SBATCH -p bigmemm
#SBATCH --job-name=RMask
#SBATCH --exclusive
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=RMask_%j.err
#SBATCH --output=RMask_%j.out

module load singularity

singularity exec -B /home/jacazet/reference/makerAnnotations/aepAnnot/repMask ~/maker-plus_3.01.03.sif RepeatMasker -pa 24 -lib oligConsensi.fa.classified olig_genome.fa
