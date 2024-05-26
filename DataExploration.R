

library(dplyr)
library(ggplot2)

# add ID 
mammalianarray_default_batpredictions$BrookLabID <- sub1_IDlist$sampleid
colnames(mammalianarray_default_batpredictions)[29] <- "sampleid"

# add age 
tmp <- inner_join(mammalianarray_default_batpredictions, aged.since.2018, by = "sampleid")
