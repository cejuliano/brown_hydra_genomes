#!/bin/bash

braker.pl --genome=aep.genome.fullsoft.fa \
        --prot_seq=allPrimProts.fa \
        --bam=rna.bam \
        --etpmode \
        --softmasking \
        --cores 24 \
        --species=HyVul \
        --AUGUSTUS_CONFIG_PATH=/home/jacazet/reference/makerAnnotations/aepAnnot/maker_braker/braker/config \
        --AUGUSTUS_BIN_PATH=/home/Augustus/bin \
        --AUGUSTUS_SCRIPTS_PATH=/home/Augustus/scripts \
        --gff3

