#!/bin/bash

bedtools closest -D ref -a HVAEP1.genes.sorted.bed -b aep16k_boundaries.sorted.bed > genesCloseTads.bed