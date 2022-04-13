#! /bin/bash -l
#SBATCH -p med
#SBATCH --job-name=pasaP
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=pasaP.err
#SBATCH --output=pasaP.out

module load singularity

singularity exec -B /home/jacazet/reference/makerAnnotations/aepAnnot/pasa/clPasaRedo \
	~/pasa.sif /usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
        -c alignAssembly.config -C -R --CPU 12 \
        --ALIGNER gmap,blat -g clytiaG.fa -t clytia.fa.clean \
        -T -u clytia.fa --TRANSDECODER \
        -d
