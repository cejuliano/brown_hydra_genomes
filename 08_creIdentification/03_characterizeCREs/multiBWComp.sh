#!/bin/bash

multiBigwigSummary bins -b H41_1.bw H41_2.bw H41_3.bw \
	H43_1.bw H43_2.bw H43_3.bw \
	H273_1.bw H273_2.bw H273_3.bw \
	IGG_1.bw IGG_2.bw IGG_3.bw \
	../ATAC/AEP1_final_shift.bw \
	../ATAC/AEP2_final_shift.bw \
	../ATAC/AEP3_final_shift.bw \
	-o corplot.npz \
	-l H41_1 H41_2 H41_3 H43_1 H43_2 H43_3 \
	H273_1 H273_2 H273_3 IGG_1 IGG_2 IGG_3 \
	ATAC_1 ATAC_2 ATAC_3 \
	-p 6