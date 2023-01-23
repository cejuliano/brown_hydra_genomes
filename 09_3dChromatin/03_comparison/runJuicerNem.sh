#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=jLaunch
#SBATCH -c 1
#SBATCH -t 30-0
#SBATCH --mem=8G
#SBATCH --error=jLaunch_%j.err
#SBATCH --output=jLaunch_%j.out

./scripts/juicerMod.sh \
        -g nem -z Nvec200.fasta \
        -p nem.genome \
		-s DpnII \
        -A jacazet -q med -Q 30-0 -l med -L 30-0 -t 8 \
        -D /home/jacazet/reference/revision/hic \
        -d /home/jacazet/reference/revision/hic/work
