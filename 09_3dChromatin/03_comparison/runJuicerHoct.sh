#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=jLaunch
#SBATCH -c 1
#SBATCH -t 30-0
#SBATCH --mem=8G
#SBATCH --error=jLaunch_%j.err
#SBATCH --output=jLaunch_%j.out

./scripts/juicerMod.sh \
        -g hoct -z hoct.genome.fa \
        -p hoct.genome \
	-s Arima \
        -A jacazet -q med -Q 30-0 -l med -L 30-0 -t 8 \
        -D /home/jacazet/reference/revision/hic \
        -d /home/jacazet/reference/revision/hic/work
