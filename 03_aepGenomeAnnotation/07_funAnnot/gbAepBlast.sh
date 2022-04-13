#!/bin/bash

blastn -query HVAEP1.tran.longestIso.fa -db genbank -evalue 1e-80 -outfmt 6 -max_target_seqs 3 -num_threads 6 > gb2AEP.txt