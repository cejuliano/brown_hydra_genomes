#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=pbj
#SBATCH -c 60
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=pbj.err
#SBATCH --output=pbj.out

source /home/jacazet/PBSuite_15.8.24/setup.sh

Jelly.py setup --debug config.xml
