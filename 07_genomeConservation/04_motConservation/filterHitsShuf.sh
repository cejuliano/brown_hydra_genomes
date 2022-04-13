#!/bin/bash -l
#SBATCH --job-name=filt
#SBATCH -p med
#SBATCH -c 1
#SBATCH -t 7-0
#SBATCH --mem=0
#SBATCH --error=filt.err
#SBATCH --output=filt.out

bedtools intersect -f 1 -u -a aep.shufMots.bed -b consensusAEP.bed > aepShufMotifsATAC.bed
bedtools intersect -v -a aepShufMotifsATAC.bed -b HVAEP1.cds.bed > aepShufMotifsATACnCDS.bed

bedtools intersect -f 1 -u -a viridShufMotifsAEP.bed -b consensusAEP.bed > viridShufMotifsAepAtac.bed
bedtools intersect -v -a viridShufMotifsAepAtac.bed -b HVAEP1.cds.bed > viridShufMotifsAepAtacNcds.bed

bedtools intersect -f 1 -u -a 105ShufMotifsAEP.bed -b consensusAEP.bed > 105ShufMotifsAepAtac.bed
bedtools intersect -v -a 105ShufMotifsAepAtac.bed -b HVAEP1.cds.bed > 105ShufMotifsAepAtacNcds.bed

bedtools intersect -f 1 -u -a oligShufMotifsAEP.bed -b consensusAEP.bed > oligShufMotifsAepAtac.bed
bedtools intersect -v -a oligShufMotifsAepAtac.bed -b HVAEP1.cds.bed > oligShufMotifsAepAtacNcds.bed

sort -k1,1 -k2,2n viridShufMotifsAEP.bed > viridShufMotifsAEP.sort.bed
sort -k1,1 -k2,2n 105ShufMotifsAEP.bed > 105ShufMotifsAEP.sort.bed
sort -k1,1 -k2,2n oligShufMotifsAEP.bed > oligShufMotifsAEP.sort.bed
sort -k1,1 -k2,2n aepShufMotifsATACnCDS.bed > aepShufMotifsATACnCDS.sort.bed


bedtools intersect -r -f 0.8 -s -wa -wb -sorted -a aepShufMotifsATACnCDS.sort.bed -b viridShufMotifsAEP.sort.bed > viridShufAepOlap.bed
bedtools intersect -r -f 0.8 -s -wa -wb -sorted -a aepShufMotifsATACnCDS.sort.bed -b oligShufMotifsAEP.sort.bed > oligShufAepOlap.bed
bedtools intersect -r -f 0.8 -s -wa -wb -sorted -a aepShufMotifsATACnCDS.sort.bed -b 105ShufMotifsAEP.sort.bed > 105ShufAepOlap.bed

awk 'BEGIN {OFS="\t"}; {if ($4 == $10) {print}}' 105ShufAepOlap.bed > 105ShufAepOlapCon.bed
awk 'BEGIN {OFS="\t"}; {if ($4 == $10) {print}}' oligShufAepOlap.bed > oligShufAepOlapCon.bed
awk 'BEGIN {OFS="\t"}; {if ($4 == $10) {print}}' viridShufAepOlap.bed > viridShufAepOlapCon.bed
