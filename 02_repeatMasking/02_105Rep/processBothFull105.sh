#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=procRep
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=process.err
#SBATCH --output=process.out

../../RepeatMasker/ProcessRepeats -a -species eumetazoa -gff bothMaskFull105.cat
