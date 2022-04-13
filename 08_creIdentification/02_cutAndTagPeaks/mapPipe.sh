#!/bin/bash

source ~/.bash_profile

prefix="$1"
echo "$prefix"

[ ! -d "fastqc_raw" ] && mkdir fastqc_raw

echo "initial fastqc report"
fastqc -o fastqc_raw -t 2 resources/reads/raw/"$prefix"_R*.fastq.gz

echo "filtering reads"
java -jar resources/trimmomatic-0.36.jar PE -threads 24 -phred33 \
	resources/reads/raw/"$prefix"_R1.fastq.gz \
	resources/reads/raw/"$prefix"_R2.fastq.gz \
	resources/reads/trim/"$prefix"_R1_trim_P.fastq.gz \
	resources/reads/trim/"$prefix"_R1_trim_UP.fastq.gz \
	resources/reads/trim/"$prefix"_R2_trim_P.fastq.gz \
	resources/reads/trim/"$prefix"_R2_trim_UP.fastq.gz \
	ILLUMINACLIP:resources/NexteraPE-PE.fa:2:30:10:2:keepBothReads \
	TRAILING:3 \
	LEADING:3 \
	SLIDINGWINDOW:4:15 \
	MINLEN:32

rm resources/reads/trim/"$prefix"*trim_UP.fastq.gz

[ ! -d "fastqc_filt" ] && mkdir fastqc_filt

echo "post-filtering fastqc report" 
fastqc -o fastqc_filt -t 2 resources/reads/trim/"$prefix"_R*_trim_P.fastq.gz

echo "mapping to genome reference"
bowtie2 -X 700 \
	-I 10 \
	--very-sensitive-local \
	--threads 22 \
	--no-mixed \
	--no-discordant \
	-x resources/hydra_genome \
	-1 resources/reads/trim/"$prefix"_R1_trim_P.fastq.gz \
	-2 resources/reads/trim/"$prefix"_R2_trim_P.fastq.gz | samtools view -b -o "$prefix".gen.bam -
	
samtools view -Su -F 524 -q 2 "$prefix".gen.bam |
	samtools sort -n -T "$prefix" -o "$prefix"_allMapped.tmp.bam

echo "mapping to spike-in reference"
bowtie2 -X 700 \
	-I 10 \
	--very-sensitive-local \
	--threads 22 \
	--no-mixed \
	--no-discordant \
	--no-overlap \
	--no-dovetail \
	-x resources/ecoli_genome \
	-1 resources/reads/trim/"$prefix"_R1_trim_P.fastq.gz \
	-2 resources/reads/trim/"$prefix"_R2_trim_P.fastq.gz |
	samtools view -Su -F 524 -q 2 - |
	samtools sort -n -T "$prefix" -o "$prefix"_ecoli_mapped.tmp.bam

echo "removing duplicates"
samtools fixmate -@ 24 -r -m "$prefix"_allMapped.tmp.bam "$prefix".fix.tmp.bam

samtools sort -T "$prefix" -@ 24 -o "$prefix".fix.sort.tmp.bam "$prefix".fix.tmp.bam

samtools markdup -@ 24 -T "$prefix" -r -s "$prefix".fix.sort.tmp.bam "$prefix".final.bam > "$prefix".dupStats.txt


samtools fixmate -@ 24 -r -m "$prefix"_ecoli_mapped.tmp.bam "$prefix".ecoli.fix.tmp.bam

samtools sort -T "$prefix" -@ 24 -o "$prefix".ecoli.fix.sort.tmp.bam "$prefix".ecoli.fix.tmp.bam

samtools markdup -@ 24 -T "$prefix" -r -s "$prefix".ecoli.fix.sort.tmp.bam "$prefix".final.ecoli.bam > "$prefix".ecoli.dupStats.txt

samtools index "$prefix".final.bam

samtools index "$prefix".final.ecoli.bam

rm "$prefix"*.tmp.bam

echo "done"

