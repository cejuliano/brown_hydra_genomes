#! /bin/bash -l

shopt -s nullglob

echo "pooling replicates"

samtools merge -f IGG_MG.final.bam IGG_[1-9].final.bam
samtools index IGG_MG.final.bam

for arg in 1 2 3 "MG"
do
	samtools sort -T IGG.sort -n -o IGG_"$arg".final.ns.bam IGG_"$arg".final.bam

	bedtools bamtobed -bedpe -i IGG_"$arg".final.ns.bam > IGG.bed

	cut -f 1,2,6 IGG.bed | \
		sort -k1,1 -k2,2n -k3,3n > IGG.fragments.bed

	bedtools genomecov -bg -i IGG.fragments.bed -g resources/aep.genome > IGG_"$arg".bg

	rm IGG.bed IGG.fragments.bed IGG_"$arg".final.ns.bam
done
