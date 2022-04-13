#!/bin/bash
#SBATCH --job-name=CnT_Peak
#SBATCH -c 24
#SBATCH -p med
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=CnT_Peak.err
#SBATCH --output=CnT_Peak.out

echo "processing IGG samples"
./resources/iggProcessing.sh

echo "running peak pipeline"
./resources/peakPipe.sh H41

./resources/peakPipe.sh H43

./resources/peakPipe.sh H273
