#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=repland
#SBATCH -c 1
#SBATCH -t 60-0
#SBATCH --mem=16G
#SBATCH --error=repland.err
#SBATCH --output=repland.out

perl ../../RepeatMasker/util/calcDivergenceFromAlign.pl -s olig.divsum olig_genome_combined.fa.align
perl ../../RepeatMasker/util/createRepeatLandscape.pl -g 1274417049 -div olig.divsum > oligRepLand.html 
