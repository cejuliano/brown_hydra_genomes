#! /bin/bash
#SBATCH -p med
#SBATCH --job-name=mmi
#SBATCH -c 4
#SBATCH -t 60-0
#SBATCH --mem=50G
#SBATCH --error=mmi.err
#SBATCH --output=mmi.out

minimap2 -d aepChr.gapfill.pbref.mmi -x map-pb aepChr.gapfill.fa
minimap2 -d aepChr.gapfill.onref.mmi -x map-ont aepChr.gapfill.fa
