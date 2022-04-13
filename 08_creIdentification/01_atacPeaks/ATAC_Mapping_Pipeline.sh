#! /bin/bash

# This code takes raw ATAC-seq reads in a fastq.gz format and generates a BAM file
# This BAM file contains only unambiguously mapped, unique, non-mitochondrial reads
# These BAM files are then used for downstream peak calling and differential 
# accessibility analyses

# This code requires:

# Trimmomatic (V. 0.36)
# Bowtie2 (V. 2.2.6)
# Fastqc (V. 0.11.4)
# Picard Tools (V. 2.17.8)
# Samtools (V. 1.2)
# Bedtools (V. 2.25)

# Read files should be formated as such:
# prefix_ATAC_R#.fastq.gz

# R# refers to the read number in the read pair

# Prefix is the unique sample identifier that has the following format:
# Time[Treatment]StructureRep

# As an example, the first biological replicate for head regenerates collected 
# at 12 hours post amputation treated with iCRT14 would have the prefix:
# 12iH1

# The untreated equivalent would have the prefix:
# 12H1


########################


# The only argument needed for this script is the unique sample ID
prefix="$1"

# create an initial fastqc report of both forward and reverse reads
# we'll compare this to the filtered reads to ensure that trimming 
# was effective 
echo "initial fastqc report"
fastqc -o . -t 2 resources/"$prefix"_ATAC_R*.fastq.gz

mv "$prefix"_ATAC_R1_fastqc.html "$prefix"_ATAC_R1_unfiltered_fastqc.html
mv "$prefix"_ATAC_R2_fastqc.html "$prefix"_ATAC_R2_unfiltered_fastqc.html

mv "$prefix"_ATAC_R1_fastqc.zip "$prefix"_ATAC_R1_unfiltered_fastqc.zip
mv "$prefix"_ATAC_R2_fastqc.zip "$prefix"_ATAC_R2_unfiltered_fastqc.zip

# use trimmomatic to remove adapter sequences and
# stretches of low quality base calls.
# It also sets the minimum read size to 32
echo "filtering reads"
java -jar resources/trimmomatic-0.36.jar PE -threads 16 -phred33 \
	resources/"$prefix"_ATAC_R1.fastq.gz resources/"$prefix"_ATAC_R2.fastq.gz \
	"$prefix"_ATAC_R1_trim_P.fastq.gz "$prefix"_ATAC_R1_trim_UP.fastq.gz \
	"$prefix"_ATAC_R2_trim_P.fastq.gz "$prefix"_ATAC_R2_trim_UP.fastq.gz \
	ILLUMINACLIP:resources/NexteraPE-PE.fa:2:30:10:2:true TRAILING:3 \
	SLIDINGWINDOW:4:15 MINLEN:32
	
rm "$prefix"*trim_UP.fastq.gz
	
# create a fastqc report of filtered reads
# this ensures filtering worked well

# you should see: high phred scores across the entire
# read length, a read length distribution that reflects
# expectations for ATAC-seq data, and no overrepresented sequences

# looking at the GC distribution can also give you a
# sense of bacterial contamination (curvibacter GC% is around 60)
echo "post-filtering fastqc report" 
fastqc -o . -t 2 "$prefix"_ATAC_R*_trim_P.fastq.gz

# map the filtered reads to the H. mag 105 genome assembly (dovetail 1.0)
# here I'm opting for quite sensitive, but relatively slow mapping

# This will require that you generate a bowtie2 reference in the
# resources folder using the H. mag 105 reference genome
# That can be found here:
# research.nhgri.nih.gov/hydra/download/?dl=asl
echo "mapping to genome reference"
bowtie2 -X 1000 --very-sensitive-local --mm --threads 16 -S "$prefix"_genome \
	-x resources/hydra_genome \
	-1 "$prefix"_ATAC_R1_trim_P.fastq.gz -2 "$prefix"_ATAC_R2_trim_P.fastq.gz
	
cat "$prefix".log "$prefix"_genome.metfile

rm "$prefix"_genome.metfile

# because the H. mag genome has mitochondrial sequence contamination,
# we need to map reads separately to the mitochondrial reference
# this let's us identify which reads are mitochondrial, so they
# can be removed downstream

# we only need the read IDs, so we just pull the mapped read 
# IDs and disregard all other output
echo "mapping to mitochondrial reference"
bowtie2 -X 1000 --very-sensitive-local --mm --threads 12 -x resources/hydra_mito \
	-1 "$prefix"_ATAC_R1_trim_P.fastq.gz -2 "$prefix"_ATAC_R2_trim_P.fastq.gz | \
	samtools view -S -F 4 - | cut -f 1 | sort | uniq > "$prefix"_mito_IDs.txt
	
cat "$prefix".log "$prefix"_mito.metfile

rm "$prefix"_mito.metfile

# We then filter out mitochondrial reads from the genome mapped reads
# using the read IDs we generated above
echo "removing mitochondrial reads"
java -Xmx16g -jar resources/picard.jar FilterSamReads \
        I="$prefix"_genome \
        O="$prefix"_genome_MF.sam \
        READ_LIST_FILE="$prefix"_mito_IDs.txt \
        FILTER=excludeReadList

# Here we filter reads to only include mapped read pairs with a 
# mapping quality of 3 or greater 

# We then sort them by name, which is required for the fixmate command
# we're performing next
echo "filtering ambiguous mappers and incoherent read pairs"
samtools view -Su -F 524 -q 3 "$prefix"_genome_MF.sam | samtools sort -n -T "$prefix" -o "$prefix".bam

# fixmate fills in information about read pairs, such as the distance 
# between R1 and R2
samtools fixmate -r "$prefix".bam "$prefix".fix.tmp.bam

# We'll again filter, this time removing improperly paired reads
# We then sort by coordinate
samtools view -u -F 524 -f 2 "$prefix".fix.tmp.bam | samtools sort -T "$prefix".fix -o "$prefix".fix.bam

rm "$prefix".fix.tmp.bam "$prefix".bam "$prefix"_genome \
	"$prefix"_mito_IDs.txt "$prefix"_genome_MF.sam

# We'll then use picard to mark perfectly duplicated read pairs
# so that they can be removed later
echo "marking PCR duplicates"
java -Xmx50g -jar resources/picard.jar MarkDuplicates \
        I="$prefix".fix.bam \
        O="$prefix"_markedDup.bam \
        M="$prefix"_markedDup_metrics.txt \
        ASSUME_SORT_ORDER=coordinate
        
# the code below  creates a report giving you a sense of how many PCR duplicates
# were found in the data
# these metrics were established by ENCODE and give a sense of how
# much bottlenecking occured during the library prep

# this code is taken from the ENCODE ATAC pipeline (v1)
# it can be found here:
# encodeproject.org/pipelines/ENCPL792NWO/

# the table is formated as such:
# TotalReadPairs [tab] DistinctReadPairs [tab] OneReadPair [tab] 
# TwoReadPairs [tab] NRF=Distinct/Total [tab] PBC1=OnePair/Distinct [tab] 
# PBC2=OnePair/TwoPair

# for more info on these metrics, see:
# encodeproject.org/atac-seq/#standards
echo "generating PBC report"
bedtools bamtobed -i "$prefix"_markedDup.bam \
        | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' \
        | sort \
        | uniq -c \
        | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > "$prefix"_PBC.txt

# We then create the final bam for downstream analysis by removing PCR duplicates
echo "removing PCR duplicates"
samtools view -F 1804 -b -o "$prefix"_final.bam "$prefix"_markedDup.bam

rm "$prefix"_markedDup.bam "$prefix".fix.bam

echo "done"

