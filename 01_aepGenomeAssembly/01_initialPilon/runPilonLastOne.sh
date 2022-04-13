#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=Pilon
#SBATCH -c 60
#SBATCH -t 60-0
#SBATCH --mem=350G
#SBATCH --error=pilon_%j.err
#SBATCH --output=pilon_%j.out

module load pilon

COUNTER=22

conts=$(head -n 22 ../resources/references/canuDraft/contGroups.txt | tail -n 1)

cut -d " " -f 1 "$conts"


samtools faidx ../resources/references/canuDraft/hydra_aep.canu.contigs.fasta $conts > subGenome.fa

java -Xmx300G -jar /share/apps/pilon-1.23/pilon-1.23.jar \
		--genome subGenome.fa --bam possorted_bam.bam \
		--output $COUNTER --outdir pilOut --threads 60

rm subGenome.fa
