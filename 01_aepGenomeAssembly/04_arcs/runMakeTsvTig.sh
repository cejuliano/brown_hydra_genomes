#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=arcs
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=36G
#SBATCH --error=arcs_py.err
#SBATCH --output=arcs_py.out

python3 ../resources/arcs/Examples/makeTSVfile.py arcs_original.gv arcs.tigpair_checkpoint.tsv ../resources/references/purgedCanPilTig/curated.fasta
