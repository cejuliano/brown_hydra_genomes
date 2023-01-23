#!/bin/bash
#SBATCH -p med
#SBATCH --job-name=ACA
#SBATCH --exclusive
#SBATCH -t 30-0
#SBATCH --mem=0
#SBATCH --error=ACA.err
#SBATCH --output=ACA.out

dirAr=( workAEP workAmil workHoct workResc workDili workNvec200 )
specAr=( aep amil hoct resc dili Nvec200 )


for i in {1..6..1}; do
	dirUse=${dirAr[$i]}
	specUse=${specAr[$i]}
	echo "$dirUse" "$specUse"
	./pseudoACA.sh "$specUse" "$dirUse"
done
