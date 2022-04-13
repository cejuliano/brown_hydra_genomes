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
Jelly.py mapping --debug config.xml
Jelly.py support --debug config.xml
Jelly.py extraction --debug config.xml
Jelly.py assembly config.xml -x "-p 10000000 -n 60 -w 1000000000 --debug"
Jelly.py output --debug config.xml
