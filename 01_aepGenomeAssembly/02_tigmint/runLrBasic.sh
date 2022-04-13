#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=lRangerB
#SBATCH -c 60
#SBATCH -t 60-0
#SBATCH --mem=350G
#SBATCH --error=lRanger_basic_%j.err
#SBATCH --output=lRanger_basic_%j.out

../resources/longranger-2.2.2/longranger basic \
	--id=aep10x_b \
	--fastqs=../resources/reads/10x/ \
	--sample=Hydra \
	--localcores=60 \
	--localmem=350
