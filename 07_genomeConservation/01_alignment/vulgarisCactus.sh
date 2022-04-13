#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=60-0
#SBATCH --partition=bigmemm
#SBATCH --error=cactus.err
#SBATCH --output=cactus.out
#SBATCH --job-name=cactus

source venv/bin/activate

cactus jobstore evolverVulgaris.txt evolverHydra.hal --realTimeLogging
