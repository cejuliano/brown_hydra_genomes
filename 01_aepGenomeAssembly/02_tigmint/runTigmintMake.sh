#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=tig
#SBATCH -c 60
#SBATCH -t 60-0
#SBATCH --mem=200G
#SBATCH --error=tigmint_%j.err
#SBATCH --output=tigmint_%j.out

source ../resources/venv/bin/activate

module load bwa

../resources/tigmint/bin/tigmint-make-mod tigmint draft=draft reads=reads
