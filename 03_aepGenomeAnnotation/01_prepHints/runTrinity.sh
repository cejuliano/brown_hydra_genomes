#! /bin/bash -l
#SBATCH -p med
#SBATCH --job-name=trinity
#SBATCH -c 24
#SBATCH -t 60-0
#SBATCH --mem=0
#SBATCH --error=trinity.err
#SBATCH --output=trinity.out

source ~/perl5/perlbrew/etc/bashrc

source venv/bin/activate

module load jellyfish

which perl

which bowtie2

$TRINITY_HOME/Trinity --genome_guided_bam ../align/out/aepAligned.sortedByCoord.out.bam \
         --genome_guided_max_intron 20000 \
         --max_memory 60G --CPU 24 \
	 --SS_lib_type RF
