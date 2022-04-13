#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=lRanger
#SBATCH -c 32
#SBATCH -t 60-0
#SBATCH --mem=300G
#SBATCH --error=lRanger_align.err
#SBATCH --output=lRanger_align.out

../resources/longranger-2.2.2/longranger align \
	--id=Canu10xMap \
	--fastqs=../resources/reads/10x/ \
	--sample=Hydra \
	--localcores=32 \
	--localmem=300 \
	--reference=../resources/references/canuDraft/refdata-hydra_aep.canu.contigs


