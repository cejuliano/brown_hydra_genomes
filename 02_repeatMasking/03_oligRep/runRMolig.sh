#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=RM
#SBATCH --exclusive
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=RM_%j.err
#SBATCH --output=RM_%j.out

module load singularity

./dfam-tetools.sh --singularity -- RepeatModeler -database olig -pa 8 -LTRStruct
