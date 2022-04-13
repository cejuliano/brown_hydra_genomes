#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=jLaunch
#SBATCH -c 1
#SBATCH -t 60-0
#SBATCH --mem=8G
#SBATCH --error=jLaunch_%j.err
#SBATCH --output=jLaunch_%j.out

./scripts/juicerMod.sh \
        -g pbj -z ../resources/references/pbj/jelly.shrink.fasta \
        -p ../resources/references/pbj/jelly.shrink.genome \
        -q med -Q 60-0 -l med -L 60-0 -t 8 \
        -D /home/jacazet/reference/aepAssembly/06_HiC \
        -d /home/jacazet/reference/aepAssembly/06_HiC/work

