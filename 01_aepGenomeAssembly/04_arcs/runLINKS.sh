#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=arcs
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=36G
#SBATCH --error=arcs_LINKS.err
#SBATCH --output=arcs_LINKS.out

LINKS -f ../resources/references/purgedCanPilTig/curated.fasta -s empty.fof -k 15 -b arcs -l 5 -t 2 -a 0.3
