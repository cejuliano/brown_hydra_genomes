#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 11 18:31:26 2021

@author: Jcazet
"""
from Bio import AlignIO
import numpy as np
import pandas as pd
import glob

files = glob.glob('../*.maf')

for file in files:
    #creates a generator
    alignment = AlignIO.parse(file, "maf")

    #initialize empty results table
    conDF = []

    i = 0

    #go through each chunk of the maf file
    for subA in alignment:
        
        i += 1
        
        if i % 1000 == 0:
            print(i)
        
        #get the list the species for each row
        species = [sq.id.split(sep='.')[0] for sq in subA._records]
        
        if i == 1:
          chrName = [sq.id.split(sep='.')[1] for sq in subA._records][0]
        
        #if AEP is the only sequence in the chunk, just write all zeros and move on
        #or if olig is not in species list
        if len(species) == 1 or 'olig' not in species:
            outL = subA.get_alignment_length()
            align_array = [0] * outL
            #print(align_array)
            conDF.extend(align_array)
            continue

        #convert alignment to character array
        align_array = np.array([list(rec) for rec in subA], str)
            
        #make all characters uppercase (so softmasking doesn't mess with checking for a match)
        align_array = np.char.upper(align_array)

        #check which bases match the AEP sequence
        align_array = pd.DataFrame(align_array[0,:] == align_array[1:,:])

        #get the species names for the match matrix
        align_array['species'] = species[1:]
        
        #drop non-olig hits
        align_array = align_array.loc[lambda df: df['species'] == 'olig', :]

        #in cases where multiple sequences from a species are in a chunk,
        #just look to see if there is any match at a given spot
        #This allows us to collapse species into single rows
        align_array = align_array.groupby('species',0).any()
        
        #count the number of species that matched at each position
        align_array = align_array.sum(0).tolist()
        
        #print(align_array)
        
        #add to results dataframe
        conDF.extend(align_array)
        
    conBG = pd.DataFrame(data={'chrom': chrName,
                               'start': range(len(conDF)),
                               'end': range(1,len(conDF) + 1),
                               'score' : conDF
                               }
                         )

    conBG.to_csv(chrName + '.olig.cactus.bedgraph',sep='\t',header=False,index=False)


    conBG_roll = conBG.copy()

    conBG_roll['score'] = conBG_roll['score'].rolling(window=100,min_periods=1).mean()

    conBG_roll.to_csv(chrName + '.olig.rolling10.cactus.bedgraph',sep='\t',header=False,index=False)



