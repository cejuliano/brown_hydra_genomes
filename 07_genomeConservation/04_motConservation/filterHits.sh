#!/bin/bash -l
#SBATCH --job-name=filt
#SBATCH -p med
#SBATCH -c 1
#SBATCH -t 7-0
#SBATCH --mem=0
#SBATCH --error=filt.err
#SBATCH --output=filt.out

bedtools intersect -f 1 -u -a aep.mots.bed -b consensusAEP.bed > aepMotifsATAC.bed
bedtools intersect -v -a aepMotifsATAC.bed -b HVAEP1.cds.bed > aepMotifsATACnCDS.bed

bedtools intersect -f 1 -u -a viridMotifsAEP.bed -b consensusAEP.bed > viridMotifsAepAtac.bed
bedtools intersect -v -a viridMotifsAepAtac.bed -b HVAEP1.cds.bed > viridMotifsAepAtacNcds.bed

bedtools intersect -f 1 -u -a 105MotifsAEP.bed -b consensusAEP.bed > 105MotifsAepAtac.bed
bedtools intersect -v -a 105MotifsAepAtac.bed -b HVAEP1.cds.bed > 105MotifsAepAtacNcds.bed

bedtools intersect -f 1 -u -a oligMotifsAEP.bed -b consensusAEP.bed > oligMotifsAepAtac.bed
bedtools intersect -v -a oligMotifsAepAtac.bed -b HVAEP1.cds.bed > oligMotifsAepAtacNcds.bed

sort -k1,1 -k2,2n viridMotifsAEP.bed > viridMotifsAEP.sort.bed
sort -k1,1 -k2,2n 105MotifsAEP.bed > 105MotifsAEP.sort.bed
sort -k1,1 -k2,2n oligMotifsAEP.bed > oligMotifsAEP.sort.bed
sort -k1,1 -k2,2n aepMotifsATACnCDS.bed > aepMotifsATACnCDS.sort.bed


bedtools intersect -r -f 0.8 -s -wa -wb -sorted -a aepMotifsATACnCDS.sort.bed -b viridMotifsAEP.sort.bed > viridAepOlap.bed
bedtools intersect -r -f 0.8 -s -wa -wb -sorted -a aepMotifsATACnCDS.sort.bed -b oligMotifsAEP.sort.bed > oligAepOlap.bed
bedtools intersect -r -f 0.8 -s -wa -wb -sorted -a aepMotifsATACnCDS.sort.bed -b 105MotifsAEP.sort.bed > 105AepOlap.bed

awk 'BEGIN {OFS="\t"}; {if ($4 == $10) {print}}' 105AepOlap.bed > 105AepOlapCon.bed
awk 'BEGIN {OFS="\t"}; {if ($4 == $10) {print}}' oligAepOlap.bed > oligAepOlapCon.bed
awk 'BEGIN {OFS="\t"}; {if ($4 == $10) {print}}' viridAepOlap.bed > viridAepOlapCon.bed
