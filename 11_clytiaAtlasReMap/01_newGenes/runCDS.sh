#! /bin/bash -l
#SBATCH -p med
#SBATCH --job-name=cds
#SBATCH -c 1
#SBATCH -t 60-0
#SBATCH --mem=8G
#SBATCH --error=cds.err
#SBATCH --output=cds.out

module load singularity

singularity exec -B /home/jacazet/reference/makerAnnotations/aepAnnot/pasa/clPasaRedo \
	~/pasa.sif /usr/local/src/PASApipeline/scripts/pasa_asmbls_to_training_set.dbi \
	--pasa_transcripts_fasta clPasa.sqlite.assemblies.fasta \
	--pasa_transcripts_gff3 clPasa.sqlite.pasa_assemblies.gff3
