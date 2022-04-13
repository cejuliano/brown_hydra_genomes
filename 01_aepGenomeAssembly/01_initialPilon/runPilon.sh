#! /bin/bash
#SBATCH -p bigmemh
#SBATCH --job-name=Pilon
#SBATCH -c 60
#SBATCH -t 60-0
#SBATCH --mem=350G
#SBATCH --error=pilon_%j.err
#SBATCH --output=pilon_%j.out

module load pilon

COUNTER=0

while read p; do
        echo "$p"

        let COUNTER=COUNTER+1

        echo "$COUNTER"

        samtools faidx ../resources/references/canuDraft/hydra_aep.canu.contigs.fasta $p > subGenome.fa

        java -Xmx300G -jar /share/apps/pilon-1.23/pilon-1.23.jar \
                --genome subGenome.fa --bam possorted_bam.bam \
                --output $COUNTER --outdir pilOut --threads 60

done < ../resources/references/canuDraft/contGroups.txt

rm subGenome.fa
