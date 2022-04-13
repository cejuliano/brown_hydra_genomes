#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=RM
#SBATCH -c 16
#SBATCH -t 60-0
#SBATCH --mem=360G
#SBATCH --error=RM_%j.err
#SBATCH --output=RM_%j.out

export PATH="/home/jacazet/reference/makerAnnotations/RepeatMasker:$PATH"
export PATH="/home/jacazet/reference/makerAnnotations/RepeatModeler-2.0.1:$PATH"

../../RepeatModeler-2.0.1/RepeatModeler -database aep -pa 4 -LTRStruct
