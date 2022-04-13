# Assembling and Annotating a Draft Genome Assembly for *Hydra oligactis*

This document describes our approach for assembling and annotating a draft genome for *Hydra oligactis*. This entailed assembling and polishing the genome using Oxford Nanopore reads, and generating gene models with BRAKER2 using previously published whole-animal RNA-seq data.  

[TOC]

### Downloading SRA from *Hydra oligactis* RNA-seq dataset from NCBI

To prepare a reference transcriptome to provide some intrinsic information about potential proteins for the protein prediction, we accessed publically available RNA-seq datasets from NCBI read archive. 
* SRR9924176,  SRR9924177,  SRR9924178 were male animals 
* SRR11037672, SRR11037671, SRR11037670 were female animals)  
sratoolkit v2.10.5 was used to dump the *R1 and *R2 files seperately. 

```bash
$SOFTWARE/sratoolkit.2.10.5-ubuntu64/bin/fastq-dump SRR9924176 --split-e & 
$SOFTWARE/sratoolkit.2.10.5-ubuntu64/bin/fastq-dump SRR9924177 --split-e & 
$SOFTWARE/sratoolkit.2.10.5-ubuntu64/bin/fastq-dump SRR9924178 --split-e & 
$SOFTWARE/sratoolkit.2.10.5-ubuntu64/bin/fastq-dump SRR11037672 --split-e & 
$SOFTWARE/sratoolkit.2.10.5-ubuntu64/bin/fastq-dump SRR11037671 --split-e & 
$SOFTWARE/sratoolkit.2.10.5-ubuntu64/bin/fastq-dump SRR11037670 --split-e  
```


*R1 and *R2 accross the three sequencing runs from each dataset were concatenated. 

```bash
cat SRR9924176_1.fastq SRR9924177_1.fastq SRR9924178_1.fastq > male_1_1.fastq  
cat SRR9924176_2.fastq SRR9924177_2.fastq SRR9924178_2.fastq > male_1_2.fastq  
cat SRR11037672_1.fastq SRR11037671_1.fastq SRR11037670_1.fastq > female_2_1.fastq  
cat SRR11037672_2.fastq SRR11037671_2.fastq SRR11037670_2.fastq > female_2_2.fastq
```

All four files were compressed using pigz v2.4

```bash
$SOFTWARE/pigz *male*
```

Next, the script "Trinity_and_Trinotate_final_version_PHIL" was invoked to 1) quality assess the reads, 2) correct them using rcorrector, 3) trim them, 4) curate them using TranscriptomeAssemblyTools, and 5) quality assess them again prior to being assembled using Trinity. The finished transcriptome was then quality assessed using BUSCO in transcriptome mode. 


```bash
#!/bin/bash


####################################
TR_SPECIMENNAME="Hydra_male+female"
TR_DATE="2021-04-19_phil"
####################################


#load envs
export PATH=$PATH:$SOFTWARE/rcorrector/
export PATH=$PATH:$SOFTWARE/salmon-latest_linux_x86_64/bin
export TRINITY_HOME=$SOFTWARE/trinityrnaseq
export PATH=$PATH:$SOFTWARE/FastQC
export PATH=$PATH:$SOFTWARE/TrimGalore-0.6.6
export PATH="$SOFTWARE/augustus-3.3.3/bin:$PATH"
export PATH="$SOFTWARE/augustus-3.3.3/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="$SOFTWARE/augustus-3.3.3/config/"
export PATH=$PATH:$SOFTWARE/Prodigal
export PATH=$PATH:$SOFTWARE/sepp-master
export PATH=$PATH:$SOFTWARE/hmmer-3.3/src
export BUSCO_CONFIG_FILE="$SOFTWARE/busco/config/myconfig.ini"
export PATH="$SOFTWARE/ncbi-blast-2.10.0+/bin:$PATH"
export PATH=$PATH:$SOFTWARE/pigz-2.4/

 ##Define variables
 TR_LOG="log/general.log"
 TR_RESULTS_DIR="${TR_SPECIMENNAME}_trinity_${TR_DATE}_running"
 TR_FINAL_DIR="${TR_SPECIMENNAME}_trinity_${TR_DATE}"

timestamp() {
  date +"%Y-%m-%d_%H-%M-%S"
}


## Create folder structures 
if [ ! -d log ]; then
    mkdir log
    mkdir toBeDeleted
    mkdir ${TR_RESULTS_DIR}
    mkdir ${TR_RESULTS_DIR}/fastqc
    mkdir ${TR_RESULTS_DIR}/busco
    mkdir ${TR_RESULTS_DIR}/rcorrector
    mkdir ${TR_FINAL_DIR}
    mkdir ${TR_FINAL_DIR}/fastqc
    mkdir ${TR_FINAL_DIR}/busco
    mkdir ${TR_FINAL_DIR}/trinotate
	echo "All directories were created." 2>&1 | tee -a $TR_LOG
fi


##versions of the software:
echo "******************************************************" 2>&1 | tee -a $TR_LOG
echo "              Version of this script                  " 2>&1 | tee -a $TR_LOG
echo "         [12] 12. March  2021, Philip Bertemes        " 2>&1 | tee -a $TR_LOG
echo " parts of this script originate from T. Ostermann     " 2>&1 | tee -a $TR_LOG
echo "  *Trinity v.2.12 installed, changed version test     " 2>&1 | tee -a $TR_LOG
echo "  *Busco 5.0.0 installed                              " 2>&1 | tee -a $TR_LOG
echo "******************************************************" 2>&1 | tee -a $TR_LOG
echo "******************************************************" 2>&1 | tee -a $TR_LOG
echo "   Versions of the software used in this pipeline     " 2>&1 | tee -a $TR_LOG
echo "******************************************************" 2>&1 | tee -a $TR_LOG
echo "Rcorrector: commit ce5d06b" 2>&1 | tee -a $TR_LOG
echo "******************************************************" 2>&1 | tee -a $TR_LOG
salmon -v 2>&1 | tee -a $TR_LOG
echo "******************************************************" 2>&1 | tee -a $TR_LOG
$TRINITY_HOME/Trinity --version 2>&1 | tee -a $TR_LOG 
echo "******************************************************" 2>&1 | tee -a $TR_LOG
fastqc -v 2>&1 | tee -a $TR_LOG
echo "******************************************************" 2>&1 | tee -a $TR_LOG
trim_galore -v 2>&1 | tee -a $TR_LOG
echo "******************************************************" 2>&1 | tee -a $TR_LOG
augustus --version 2>&1 | tee -a $TR_LOG
echo "******************************************************" 2>&1 | tee -a $TR_LOG
prodigal -v 2>&1 | tee -a $TR_LOG 
echo "******************************************************" 2>&1 | tee -a $TR_LOG
echo "sepp: commit bd26318" 2>&1 | tee -a $TR_LOG 
echo "******************************************************" 2>&1 | tee -a $TR_LOG
echo "hmmer: HMMER 3.3 November 2019)" 2>&1 | tee -a $TR_LOG 
echo "******************************************************" 2>&1 | tee -a $TR_LOG
blastn -version 2>&1 | tee -a $TR_LOG 
echo "******************************************************" 2>&1 | tee -a $TR_LOG
busco -version 2>&1 | tee -a $TR_LOG 
echo "******************************************************" 2>&1 | tee -a $TR_LOG
echo "TranscriptomeAssemblyTools: commit e2df226" 2>&1 | tee -a $TR_LOG 
echo "******************************************************" 2>&1 | tee -a $TR_LOG
pigz --version 2>&1 | tee -a $TR_LOG 
echo "******************************************************" 2>&1 | tee -a $TR_LOG


#THIS TEST WAS ADDED IN v8 OF THIS SCRIPT
if [ ! -s ${TR_RESULTS_DIR}/fastqc/unfixrm* ]
then

##fastqc prior to corrections
timestamp 2>&1 | tee -a $TR_LOG
echo "Start with fastqc with raw data"  2>&1 | tee -a $TR_LOG
$SOFTWARE/FastQC/fastqc *fastq.gz --threads 63 --outdir ${TR_RESULTS_DIR}/fastqc 2>&1 | tee -a $TR_LOG

##counter to find pairs of PE and add pass the right variable to "-s" in FilterUncorrectabledPEfastq.py
counter=1
while [ $counter -le 5 ]
do

actualfile="_""$counter""_"
echo $actualfile  2>&1 | tee -a $TR_LOG

if ls -la | grep $actualfile 
then
    ##R-corrector
    timestamp  2>&1 | tee -a $TR_LOG
    echo "A file with ID " $actualfile " has been found in the directory. Start with rcorrector..."  2>&1 | tee -a $TR_LOG
    perl $SOFTWARE/rcorrector/run_rcorrector.pl -1 *$actualfile*R1* -2 *$actualfile*R2* -t 63 -od ${TR_RESULTS_DIR}/rcorrector 2>&1 | tee -a log/output_rcorrector.dat
    timestamp  2>&1 | tee -a $TR_LOG
    echo "Done with rcorrector on files with ID " $counter "."  2>&1 | tee -a $TR_LOG

    #TranscriptomeAssemblyTools https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html 
    echo "Start with TranscriptomeAssemblyTools on rcorrected-files with ID " $actualfile " in directory " ${TR_RESULTS_DIR}"/rcorrector/"  2>&1 | tee -a $TR_LOG
    python $SOFTWARE/TranscriptomeAssemblyTools/FilterUncorrectabledPEfastq.py -1 ${TR_RESULTS_DIR}/rcorrector/*$actualfile*R1* -2 ${TR_RESULTS_DIR}/rcorrector/*$actualfile*R2* -s ${TR_SPECIMENNAME}_$counter  2>&1 | tee -a log/output_TranscriptomeAssemblyTools.dat
    timestamp  2>&1 | tee -a $TR_LOG
    echo "Done with TranscriptomeAssemblyTools on files with ID " $actualfile "."  2>&1 | tee -a $TR_LOG
    
    ##Trim-TrimGalore https://github.com/FelixKrueger/TrimGalore/blob/master/Docs/Trim_Galore_User_Guide.md
    echo "Start with Trimgalore on rcorrected-files and TAT-fixed files with ID " $actualfile " in the main directory, named unfixrm_......"$actualfile"...cor.fq"  2>&1 | tee -a $TR_LOG
    trim_galore --paired --retain_unpaired --phred33 --output_dir trimmed_reads --length 36 -q 5 --stringency 1 -e 0.1 unfixrm*$actualfile*R1* unfixrm*$actualfile*R2* --cores 64 --gzip 2>&1 | tee -a log/output_trimgalore.dat
    timestamp  | tee -a log/output_trimgalore.dat 2>&1
    echo "Done with TrimGalore on files with ID " $actualfile "."  2>&1 | tee -a $TR_LOG
else
    echo "no pair with ID $actualfile found"   2>&1 | tee -a $TR_LOG ; 
fi
((counter++))
done

timestamp 2>&1 | tee -a $TR_LOG
echo "Concatenating all corrected and validated R1 files to left.fq.gz" 2>&1 | tee -a $TR_LOG
zcat trimmed_reads/*R1*cor_val_1.fq.gz | pigz > trimmed_reads/left.fq.gz
timestamp 2>&1 | tee -a $TR_LOG
echo "Concatenating all corrected and validated R2 files to right.fq.gz" 2>&1 | tee -a $TR_LOG
zcat trimmed_reads/*R2*cor_val_2.fq.gz | pigz > trimmed_reads/right.fq.gz

# ##fastqc after corrections
timestamp 2>&1 | tee -a $TR_LOG
echo "Start with fastqc after rcorrector, TranscriptomeAssemblyTools, and TrimGalore!"   2>&1 | tee -a $TR_LOG
$SOFTWARE/FastQC/fastqc trimmed_reads/*cor_val_*.fq.gz --threads 63 --outdir ${TR_RESULTS_DIR}/fastqc 2>&1 | tee -a log/output_fastqc.dat | tee -a $TR_LOG


# TEST ENDS HERE
fi




##Trinity
timestamp  2>&1 | tee -a $TR_LOG
echo "Start with Trinity..."   2>&1 | tee -a $TR_LOG
# ATTENTION we can use --monitor to start collectl. 
$TRINITY_HOME/Trinity \
  --seqType fq \
  --max_memory 250G \
  --CPU 63 \
  --no_salmon \
  --SS_lib_type RF \
  --left trimmed_reads/left.fq.gz \
  --right trimmed_reads/right.fq.gz \
  --output ${TR_RESULTS_DIR} 2>&1 | tee -a log/output_trinity.dat

#THIS TEST IF CONSTRUCT WAS ADDED IN V8 OF THIS SCRIPT
  if tail -100 log/output_trinity.dat | grep "Trinity assemblies are written to"
then
    timestamp 2>&1 | tee -a $TR_LOG
    echo "TRINITY ASSEMBLIES WERE FOUND. CONTINUE WITH THE PIPELINE!" 2>&1 | tee -a $TR_LOG
    
##Trinity STATS generieren
# ATTENTION we can use --monitor to start collectl. If --monitor was used in Trinity command, uncomment the examine_resource_usage_profiling.pl script in the next line
$TRINITY_HOME/util/TrinityStats.pl ${TR_RESULTS_DIR}/Trinity.fasta >> ${TR_RESULTS_DIR}/TrinityStats.txt
#$TRINITY_HOME/trinity-plugins/COLLECTL/examine_resource_usage_profiling.pl collectl


## Path-Informationen aus Fasta-File entfernen (von Trinity):
cp ${TR_RESULTS_DIR}/Trinity.fasta ./log/${TR_SPECIMENNAME}_${TR_DATE}_Trinity-ORIGINAL.fasta
sed 's/ path.*//' ${TR_RESULTS_DIR}/Trinity.fasta >> ${TR_RESULTS_DIR}/Trinity_ohnePath.fasta
sed "s/>TRINITY/>${TR_SPECIMENNAME}/" ${TR_RESULTS_DIR}/Trinity_ohnePath.fasta >> ${TR_RESULTS_DIR}/${TR_SPECIMENNAME}_${TR_DATE}.fasta

##Run the new busco v.4. ATTENTION 20200115phil_: BLAST v2.6 has an issue with multi-thread in Busco. Busco constraints tblastn to one single core!
	python3 $SOFTWARE/busco/src/busco/run_BUSCO.py \
	-o busco \
	-i ${TR_RESULTS_DIR}/${TR_SPECIMENNAME}_${TR_DATE}.fasta \
	-l metazoa_odb10 \
	-m transcriptome \
	--cpu 63 2>&1 | tee -a $TR_LOG

##this was added in verison 12 of my script
##rename the busco output
mv busco/short_summary.* busco/short_summary.specific.metazoa_odb10.${TR_SPECIMENNAME}_${TR_DATE}.txt
	
## generate Figure...
python3 $SOFTWARE/busco/scripts/generate_plot.py -wd busco 2>&1 | tee -a log/output_busco.log


if [ -s ${TR_RESULTS_DIR}/${TR_SPECIMENNAME}_${TR_DATE}.fasta ]
    then
        mv *.log log/
        mv busco* ${TR_FINAL_DIR}/busco
        mv ${TR_RESULTS_DIR}/fastqc ${TR_FINAL_DIR}/fastqc
        mv ${TR_RESULTS_DIR}/${TR_SPECIMENNAME}_${TR_DATE}.fasta ${TR_FINAL_DIR}
        mv ${TR_RESULTS_DIR}/TrinityStats.txt ${TR_FINAL_DIR}
        mv log/ ${TR_FINAL_DIR}/log
        mv unfix* toBeDeleted/
        mv *fastq* toBeDeleted/
        mv ${TR_RESULTS_DIR} toBeDeleted/
        mv trimmed_reads toBeDeleted/
        #mv *collectl* log/
    else
        echo "Trinity has not finished correctly. Please check what went wrong, and restart!" 2>&1 | tee -a $TR_LOG
fi
```


### Basecalling the Nanopore reads
Two libraries were prepared to generate long reads for *Hydra oligactis*. The first library was loaded 2 times on a Nanopore flow cell, resulting in two read files (lib1_1 and lib1_2). The second library was loaded 6 times. All 8 files were basecalled using the guppy basecaller v4.5.2 using the high accuracy moodel (HAC) over a NVidia 2070Super GPU. 

```bash
$SOFTWARE/ont-guppy/bin/guppy_basecaller --min_qscore 7 -i ../../nanopore_raw/Hydra_oligactis/Hydra_oli_SRE_XS_lib1_1/*/fast5/ -s Hydra_oli_SRE_XS_lib1_1_HAC_guppy452 -c $SOFTWARE/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg -x 'cuda:0' --num_callers 4 --gpu_runners_per_device 8 
$SOFTWARE/ont-guppy/bin/guppy_basecaller --min_qscore 7 -i ../../nanopore_raw/Hydra_oligactis/Hydra_oli_SRE_XS_lib1_2/*/fast5/ -s Hydra_oli_SRE_XS_lib1_2_HAC_guppy452 -c $SOFTWARE/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg -x 'cuda:0' --num_callers 4 --gpu_runners_per_device 8 
$SOFTWARE/ont-guppy/bin/guppy_basecaller --min_qscore 7 -i ../../nanopore_raw/Hydra_oligactis/Hydra_oli_SRE_XS_lib2_1/*/fast5/ -s Hydra_oli_SRE_XS_lib2_1_HAC_guppy452 -c $SOFTWARE/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg -x 'cuda:0' --num_callers 4 --gpu_runners_per_device 8 
$SOFTWARE/ont-guppy/bin/guppy_basecaller --min_qscore 7 -i ../../nanopore_raw/Hydra_oligactis/Hydra_oli_SRE_XS_lib2_2/*/fast5/ -s Hydra_oli_SRE_XS_lib2_2_HAC_guppy452 -c $SOFTWARE/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg -x 'cuda:0' --num_callers 4 --gpu_runners_per_device 8
$SOFTWARE/ont-guppy/bin/guppy_basecaller --min_qscore 7 -i ../../nanopore_raw/Hydra_oligactis/Hydra_oli_SRE_XS_lib2_3/*/fast5/ -s Hydra_oli_SRE_XS_lib2_3_HAC_guppy452 -c $SOFTWARE/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg -x 'cuda:0' --num_callers 4 --gpu_runners_per_device 8 
$SOFTWARE/ont-guppy/bin/guppy_basecaller --min_qscore 7 -i ../../nanopore_raw/Hydra_oligactis/Hydra_oli_SRE_XS_lib2_4/*/fast5/ -s Hydra_oli_SRE_XS_lib2_4_HAC_guppy452 -c $SOFTWARE/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg -x 'cuda:0' --num_callers 4 --gpu_runners_per_device 8 
$SOFTWARE/ont-guppy/bin/guppy_basecaller --min_qscore 7 -i ../../nanopore_raw/Hydra_oligactis/Hydra_oli_SRE_XS_lib2_5/*/fast5/ -s Hydra_oli_SRE_XS_lib2_5_HAC_guppy452 -c $SOFTWARE/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg -x 'cuda:0' --num_callers 4 --gpu_runners_per_device 8 
$SOFTWARE/ont-guppy/bin/guppy_basecaller --min_qscore 7 -i ../../nanopore_raw/Hydra_oligactis/Hydra_oli_SRE_XS_lib2_6/*/fast5/ -s Hydra_oli_SRE_XS_lib2_6_HAC_guppy452 -c $SOFTWARE/ont-guppy/data/dna_r9.4.1_450bps_hac.cfg -x 'cuda:0' --num_callers 4 --gpu_runners_per_device 8 
```

### Quality assessement of the Nanopore reads
Next, we ran NanoPlot v1.30.1 to quality assess each of the "runs". 

```bash
NanoPlot --summary Hydra_oli_SRE_XS_lib1_1_HAC_guppy452/sequencing_summary.txt --N50 --o NANOPLOT_Hydra_oli_SRE_XS_lib1_1_HAC_guppy452 -t 60
NanoPlot --summary Hydra_oli_SRE_XS_lib1_2_HAC_guppy452/sequencing_summary.txt --N50 --o NANOPLOT_Hydra_oli_SRE_XS_lib1_2_HAC_guppy452 -t 60
NanoPlot --summary Hydra_oli_SRE_XS_lib2_1_HAC_guppy452/sequencing_summary.txt --N50 --o NANOPLOT_Hydra_oli_SRE_XS_lib2_1_HAC_guppy452 -t 60
NanoPlot --summary Hydra_oli_SRE_XS_lib2_2_HAC_guppy452/sequencing_summary.txt --N50 --o NANOPLOT_Hydra_oli_SRE_XS_lib2_2_HAC_guppy452 -t 60
NanoPlot --summary Hydra_oli_SRE_XS_lib2_3_HAC_guppy452/sequencing_summary.txt --N50 --o NANOPLOT_Hydra_oli_SRE_XS_lib2_3_HAC_guppy452 -t 60
NanoPlot --summary Hydra_oli_SRE_XS_lib2_4_HAC_guppy452/sequencing_summary.txt --N50 --o NANOPLOT_Hydra_oli_SRE_XS_lib2_4_HAC_guppy452 -t 60
NanoPlot --summary Hydra_oli_SRE_XS_lib2_5_HAC_guppy452/sequencing_summary.txt --N50 --o NANOPLOT_Hydra_oli_SRE_XS_lib2_5_HAC_guppy452 -t 60
NanoPlot --summary Hydra_oli_SRE_XS_lib2_6_HAC_guppy452/sequencing_summary.txt --N50 --o NANOPLOT_Hydra_oli_SRE_XS_lib2_6_HAC_guppy452 -t 60
NanoPlot --summary Hydra_oli_SRE_XS_*/sequencing_summary.txt --N50 --o NANOPLOT_Hydra_oli_SRE_XS_lib1+2complete_HAC_guppy452 -t 60
```

### Concatenating all the QC-passed reads for assembly
All the reads from the 8 individual "runs" passing the guppy basecaller QC criteria (Q-Score >7) were now concatenated into a single fastq file. 
```bash
cat Hydra_oli_SRE_XS_lib*_HAC_guppy452/pass/Hydra*.fastq > Hydra_lib1+2_complete_q7_HACguppy452.fastq
```

### assembling of the genome
Then, the Flye v2.8.3 assembler was run to generate the genome
```bash
$SOFTWARE/Flye/bin/flye --nano-raw Hydra_lib1+2_complete_q7_HACguppy452.fastq --threads 62 -o Hydra_lib1+2_complete_q7_guppy452_flye283 
mv Hydra_lib1+2_complete_q7_guppy452_flye283/assembly.fasta Hydra_lib1+2_complete_q7_guppy452_flye283.fasta
```

Assembly on flye 283 finished with the following metrics:
  Total length:    1283640785 
     Fragments:    18685 
     Fragments N50:    272153 
     Largest frg:    6425239 
     Scaffolds:    114 
     Mean coverage:    17 

### QC with BUSCO on the newly created genome file
We checked the genome file for completeness using BUSCO v5.0.0 in genome mode

```bash
busco -i Hydra_lib1+2_complete_q7_guppy452_flye283.fasta -o busco -l metazoa_odb10 -m genome -c 62
```
        |Results from dataset metazoa_odb10                
        -------------------------------------------------- 
        |C:89.2%[S:76.1%,D:13.1%],F:5.7%,M:5.1%,n:954      
        |851    Complete BUSCOs (C)                        
        |726    Complete and single-copy BUSCOs (S)        
        |125    Complete and duplicated BUSCOs (D)         
        |54     Fragmented BUSCOs (F)                      
        |49     Missing BUSCOs (M)                         
        |954    Total BUSCO groups searched                
        -------------------------------------------------- 

### The genome file was then polished using medaka
"Raw" draft assemblies using uncorrected Nanopore reads are known to contain some errors, especially in homopolymer-regions. A software (medaka) was developped to find and correct these errors. We used the best medaka model, however, it was trained on guppy v3.0.3 basecalled data, which was not as accurate as the guppy v4.5.2 used to basecall the data from the present genome. 

. $SOFTWARE/medaka/bin/activate 
PATH=$SOFTWARE/htslib-1.11/:$PATH 
PATH=$SOFTWARE/samtools-1.11/:$PATH 
PATH=$SOFTWARE/minimap2/:$PATH 
PATH=$SOFTWARE/bcftools/:$PATH 

```bash
medaka_consensus -i Hydra_lib1+2_complete_q7_HACguppy452.fastq -d Hydra_lib1+2_complete_q7_HACguppy452_flye283.fasta -o medaka_1_on_flye_original -t 60 m r941_min_high_g303
### This script breaks before the stitch, so it has to be invoked after medaka has finished
medaka stitch medaka_1_on_flye_original/consensus_probs.hdf Hydra_lib1+2_complete_q7_HACguppy452_flye283.fasta Hydra_lib1+2_complete_q7_HACguppy452_flye283_medaka1.fasta --threads 60
```

### QC with BUSCO on the polished genome file
The polished genome was rechecked with BUSCO v5.0.0 to check for improvements

```bash
busco -i Hydra_lib1+2_complete_q7_guppy452_flye283_medaka1.fasta -o busco -l metazoa_odb10 -m genome -c 62
```

    |Results from dataset metazoa_odb10                
    -------------------------------------------------- 
    |C:90.2%[S:77.1%,D:13.1%],F:5.3%,M:4.5%,n:954
    |861    Complete BUSCOs (C)
    |736    Complete and single-copy BUSCOs (S) 
    |125    Complete and duplicated BUSCOs (D)
    |51     Fragmented BUSCOs (F)   
    |42     Missing BUSCOs (M)   
    |954 Total BUSCO groups searched
    --------------------------------------------------



### Soft-masking repetitive elements in the genome
Before generating gene models for the *oligactis* genome, we first soft-masked all repeats in the draft assembly. The process for generating the masked version of the genome (`olig_genome.sm.fa`) is described in the `02_repeatMasking.md` document.

### Mapping the de-novo transcriptome to the genome
The next part is the protein prediction on the new genome file. First, the newly generated transcriptome gets mapped to the genome using minimap2 

```bash
$SOFTWARE/minimap2/minimap2 -ax splice olig_genome.sm.fa Hydra_male+female_2021-04-19_phil.fasta -t 62  > transcriptome_malefemale.sam; $SOFTWARE/samtools-1.11/samtools sort -@ 60 transcriptome_malefemale.sam > transcriptome_malefemale.bam
```

### Protein prediction with braker2 pipeline
Finally, we started the braker2 pipeline (in a conda environment) to predict proteins

```bash
conda activate braker
braker.pl --species=Hydraoli_20211201 --genome=olig_genome.sm.fa --bam=transcriptome_malefemale.bam --cores=62 -gff3 --softmasking
```

### QC of the predicted proteins
In addition, we ran BUSCO v5.2.2 (now in a conda environment) in protein mode to check for completeness of the predicted proteins. 

```bash
conda activate busco
busco -i Hydra_lib1+2_complete_q7_guppy452_flye283_medaka1.fasta -o busco -l metazoa_odb10 -m genome -c 62
```

       |Results from dataset metazoa_odb10                
        -------------------------------------------------- 
        |C:86.4%[S:71.0%,D:15.4%],F:7.8%,M:5.8%,n:954
        |824     Complete BUSCOs (C)
        |677     Complete and single-copy BUSCOs (S)
        |147     Complete and duplicated BUSCOs (D)
        |74      Fragmented BUSCOs (F)
        |56      Missing BUSCOs (M)
        |954     Total BUSCO groups searched
        -------------------------------------------------- 


​        
### Final information
The computational needs for this complete project did not exceed the power provided by a workstation containing 
    * AMD Ryzen Threadripper 3970X 32-Core Processor
    * 256G memory
    * NVidia 2070Super

### Files Associated with This Document



```
04_oligactisDraftGenome/
├── assembly.fasta
		Initial draft assembly of the H. oligactis genome generated by Flye using Nanopore reads.
├── assemblyPolish.fasta
		Updated version of the H. oligactis draft genome following error correction with Medaka
		using Nanopore data.
├── augustusSpeciesHyOli.tar.gz
		Folder containing the optimized Augustus gene prediction parameters for the H. oligactis 
		assembly generated by BRAKER2.
├── braker.gff3
		GFF3-formatted file containing the Initial gene models for the H. oligactis assembly 
		generated by Braker2 using previously published whole-animal RNA-seq data as hints.
├── braker.gtf
		GTF-formatted file containing the Initial gene models for the H. oligactis assembly 
		generated by Braker2 using previously published whole-animal RNA-seq data as hints.
└── Hydra_male+female_2021-04-19_phil.fasta
		De-novo H. oligactis transcriptome generated with Trinity using previously published
		whole-animal RNA-seq data. Used as hints for the BRAKER2 annotation pipeline.
```

