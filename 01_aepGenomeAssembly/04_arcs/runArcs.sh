#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=arcs
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=arcs.err
#SBATCH --output=arcs.out

arcs -f ../resources/references/purgedCanPilTig/curated.fasta -s 95 -c 5 -l 0 -z 500 \
        -m 20-100000 -d 0 -e 200000 -r 0.05 -v -b arcs nameSorted.bam
