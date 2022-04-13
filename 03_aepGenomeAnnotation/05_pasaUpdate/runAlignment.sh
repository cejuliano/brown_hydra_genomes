#! /bin/bash

/usr/local/src/PASApipeline/Launch_PASA_pipeline.pl \
        -c alignAssembly.config -C -R --CPU 30 \
        --ALIGNER gmap,blat -g aep.final.genome.fa -t Trinity-GG.fasta.clean \
        -T -u Trinity-GG.fasta --TRANSDECODER \
        --stringent_alignment_overlap 30.0 -d
