## Script to plot bedgraph plots from bismark output
## Sophia Horigan
## Oct 26 2022
library(methrix)
install.packages(methrix)
library(Sushi)
install.packages(Sushi)

## Load in bismark files
homewd= "~/Desktop/Epigenetic_Clock/methylationAging/"

setwd(paste0(homewd, "/output/bismark-bedgraph/"))

data = consensus_probe_sequences_MT_CpG

chrom = 
chromstart = 
chromend = 

plotBedgraph(data)
