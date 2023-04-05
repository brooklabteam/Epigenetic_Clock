#This R script is for general analysis of known-aged bats
#these are samples that we will use for building our clocks, or our clock samples


#There are multiple types of aged bats relevant to these analyses
#1) tooth-aged bats
#2) young-of-year-aged bats
#3) bats sampled before wing punches started being taken (DATE???)
#4) bats for which Luminex was run (DATE???)


## Each section will have code to total subsections of this data in order to compile lists of samples
## to be used for various analyses.

homewd = "/Users/sophiahorigan/Documents/Github/Epigenetic_Clock/"
setwd(paste0(homewd,"/Aged-Samples/"))

aged_bats <- read.csv(file="test.csv", header = T, stringsAsFactors = F)

######################################
### Section 1: Samples for VirScan ###
######################################
# Goal: Bats sampled from 2018 onwards, for which Luminex was previously run.

aged_bats$collection_date <- as.Date(aged_bats$collection_date, format="%Y-%m-%d")
# something is wrong, there are NAs


# subset bats sampled in 2018 or later
aged.since.2018 <- aged_bats[aged_bats$collection_date > "2017-12-31",]


df <- df[rowSums(is.na(df)) != ncol(df), ]
aged.since.2018 <- aged.since.2018[rowSums(is.na(aged.since.2018)) != ncol(aged.since.2018),]

write.csv(aged.since.2018, file = "aged-since-2018.csv", row.names = F)






