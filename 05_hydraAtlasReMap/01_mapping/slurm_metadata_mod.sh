#!/bin/bash
#SBATCH --job-name=metadata
#SBATCH -c 16 
#SBATCH -t 1-0
#SBATCH --mem=24G
#SBATCH --error=logs/metadata.err
#SBATCH --output=logs/metadata.out

	module load star
	module load samtools

        bash create_Drop-seq_reference_metadata_mod.sh \
	-n HVAEP1.tran.final.mito.ds \
        -r HVAEP1.tran.final.mito.fasta \
        -s Hydra \
        -g HVAEP1.transcriptome.mito.gtf \
        -v \
        -d /group/julianolab/analyses/dropseq/Drop-seq_tools-2.4.0 \
        -t /group/julianolab/analyses/dropseq/HVAEP1_transcriptome_final/dropseq/temp
