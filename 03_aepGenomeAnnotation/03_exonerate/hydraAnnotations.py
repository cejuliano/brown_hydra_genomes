#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 14 10:04:21 2020

@author: Jcazet
"""
from Bio import SeqIO
import re

gb = SeqIO.parse('/Users/Jcazet/Google_Drive/Juliano_lab/References/genbank/"Hydra vulgaris"[porgn] AND (biomol_mrna[PROP] AND ddbj_embl_genbank[filter]).gb', "genbank")

# c = 0

gbKeep = []

while True:
    try:
        f = next(gb)
    except StopIteration:
        break
    else:
        tests = []
        tests.append(re.search('est project', f.annotations['references'][0].title, re.IGNORECASE))
        tests.append(re.search('rna interference', f.annotations['references'][0].title, re.IGNORECASE))
        tests.append(re.search('Hydra vulgaris cDNAs', f.annotations['references'][0].title, re.IGNORECASE))
        tests.append(re.search('Comparative analysis of septic', f.annotations['references'][0].title, re.IGNORECASE))
        
        if all(v is None for v in tests):
            f.description = re.sub(r'\, .*cds$',r'',f.description, re.IGNORECASE)
           # print(f.description)
            f.description = f.description + '|' + f.annotations['references'][0].title
            f.id = f.id + '|'
            gbKeep.append(f)


SeqIO.write(gbKeep, "/Users/Jcazet/Google_Drive/Juliano_lab/References/genbank/hydraAnnotations.gb","gb")
SeqIO.write(gbKeep, "/Users/Jcazet/Google_Drive/Juliano_lab/References/genbank/hydraAnnotations.fasta","fasta")

gbKeep[60].annotations
