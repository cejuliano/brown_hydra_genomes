#!/bin/bash -l
#SBATCH --job-name=window
#SBATCH -p med
#SBATCH -c 1
#SBATCH -t 7-0
#SBATCH --mem=16G
#SBATCH --error=window.err
#SBATCH --output=window.out

python mafWindows105.py
python mafWindowsVirid.py
python mafWindowsOlig.py
python mafWindowsClytia.py
