#!/bin/bash -l
#SBATCH --job-name=bw
#SBATCH -p bigmemm 
#SBATCH -c 1
#SBATCH -t 7-0
#SBATCH --mem=0
#SBATCH --error=bw.err
#SBATCH --output=bw.out

cat *olig.rolling* > olig.rolling10.cactus.bedgraph
bedSort olig.rolling10.cactus.bedgraph olig.rolling10.cactus.sort.bedgraph
bedGraphToBigWig olig.rolling10.cactus.sort.bedgraph ../../aep.genome olig.rolling10.cactus.bw

cat *virid.rolling* > virid.rolling10.cactus.bedgraph
bedSort virid.rolling10.cactus.bedgraph virid.rolling10.cactus.sort.bedgraph
bedGraphToBigWig virid.rolling10.cactus.sort.bedgraph ../../aep.genome virid.rolling10.cactus.bw

cat *105.rolling* > 105.rolling10.cactus.bedgraph
bedSort 105.rolling10.cactus.bedgraph 105.rolling10.cactus.sort.bedgraph
bedGraphToBigWig 105.rolling10.cactus.sort.bedgraph ../../aep.genome 105.rolling10.cactus.bw

cat *clytia.rolling* > clytia.rolling10.cactus.bedgraph
bedSort clytia.rolling10.cactus.bedgraph clytia.rolling10.cactus.sort.bedgraph
bedGraphToBigWig clytia.rolling10.cactus.sort.bedgraph ../../aep.genome clytia.rolling10.cactus.bw

