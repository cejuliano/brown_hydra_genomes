#!/bin/bash

bedtools genomecov -i bothMaskCplx.out.gff -bga -g ../aep.genome > repDensity.bg

bedGraphToBigWig repDensity.bg ../aep.genome repDensity.bw