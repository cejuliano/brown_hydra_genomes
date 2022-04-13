#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=repland
#SBATCH -c 1
#SBATCH -t 60-0
#SBATCH --mem=16G
#SBATCH --error=repland.err
#SBATCH --output=repland.out

perl ../../RepeatMasker/util/calcDivergenceFromAlign.pl -s 105.divsum bothMaskFull105.align
perl ../../RepeatMasker/util/createRepeatLandscape.pl -g 853782670 -div 105.divsum > 105RepLand.html 
