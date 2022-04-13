#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=lRanger
#SBATCH -c 60
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=lRanger_align.err
#SBATCH --output=lRanger_align.out

../resources/longranger-2.2.2/longranger align \
        --id=purgedMap \
        --fastqs=../resources/reads/10x/ \
        --sample=Hydra \
        --localcores=60 \
        --localmem=500 \
        --reference=../resources/references/purgedCanPilTig/refdata-curated

