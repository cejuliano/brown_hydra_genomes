#!/bin/bash
#SBATCH --job-name=makeMaf
#SBATCH -p med
#SBATCH -c 1
#SBATCH -t 60-0
#SBATCH --array=1-15
#SBATCH --mem=8G
#SBATCH --error=makeMaf_%a.err
#SBATCH --output=makeMaf_%a.out

source ../venv/bin/activate

arg=$SLURM_ARRAY_TASK_ID

echo "$arg"

hal2maf --onlyOrthologs --noAncestors --noDupes --refSequence chr-"$arg" --refGenome aep evolverHydra.hal cactus.chr"$arg".maf
