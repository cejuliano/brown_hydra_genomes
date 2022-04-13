#!/bin/bash

plotCorrelation -in corplot.npz \
	-c spearman \
	-p heatmap \
	-o corPlot.pdf \
	-min -0.05 \
	--plotNumbers \
	--colorMap magma
