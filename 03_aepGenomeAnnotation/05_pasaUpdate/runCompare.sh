#! /bin/bash
#SBATCH -p med 
#SBATCH --job-name=pasaCmp
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=60G
#SBATCH --error=pasaCmp.err
#SBATCH --output=pasaCmp.out

module load singularity

singularity exec -B /home/jacazet/reference/makerAnnotations/aepAnnot/pasa/pasaUpdate ~/pasa.sif /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
	-c annotationCompare.config -A \
	-g aep.final.genome.fa \
	-t Trinity-GG.fasta.clean \
	-L \
	--annots HVAEP1.baseline.geneModels.gff3 \
	--CPU 24
