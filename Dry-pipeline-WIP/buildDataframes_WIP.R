####################################################################################
################# BUILDING THE DATAFRAME(S) FOR CLOCK TRAINING #####################
####################################################################################

#Code for combining the bismark.cov files for each individual bat
#Direct questions to Theresa Laverty, Sophia Horigan, or Cara Brook
#tlaverty@nmsu.edu; shorigan@uchicago.edu; cbrook@uchicago.edu

# The dataframe (df) should be organized with rows = samples, columns = CpGs,... 
# ...and a few metadata columns that need to be removed before clock training.

#Note from Manny on the bismark.cov files:
# The non-“sub” samples are whole genome bisulfite seq samples, 
# and serve as a “does this sequence, and does it map to the genome 
# and/or to the probe FASTA” control. The “sub” samples represent 
# the same library as the non-sub sample, but after probe capture 
# enrichment for the Wilkinson et al CpG sites; these are the 
# representative samples for what things look like after the full pipeline.

#Using the "sub" files for now. 

# bismark file outputs ==
#   <chromosome> <start position> <end position> <methylation percentage> <count methylated> <count unmethylated>


####################################################################################
#Clear work environment
rm(list=ls())

#Set working directory - add yours here
#homewd= "/Users/theresalaverty/Documents/R/R_repositories/Epigenetic_Clock"
homewd= "/Users/shorigan/Documents/GitHub/Epigenetic_Clock"
#Folder where the bismark.cov files of interest are stored
setwd(paste0(homewd,"/methylation-aging/output/bismark_bedgraph/"))

#Load packages
#install.packages("data.table")
#install.packages("dplyr")
#install.packages("janitor")
library(data.table)
library(dplyr)
library(tidyr)
library(tidyverse)
library(janitor)
library(cowplot)

####################################################################################
#PART 1: TRANSFORM METHYLATION DATA
### 1. CpG sites ###
#find all of the *-sub_CpG.gz.bismark.cov files in the working directory
(bismark_files_cpg <- list.files(pattern="*-sub_CpG.gz.bismark.cov", recursive=TRUE))

out_cpg<-data.frame()
for(i in 1:length(bismark_files_cpg)) {
  eg<-read.table(bismark_files_cpg[i]) #read in each file in list as table
  tmp<-as.data.frame(t(eg)) #transpose table
  cols <- paste(tmp[1,], tmp[2,], sep="_") #create vector with probe site _ nucleotide #
  names(tmp) <- cols #make the vector the column names
  tmp <- tmp[-c(1:3,5:6),] #remove all rows except percent methylation data
  tmp <-  cbind("id" = bismark_files_cpg[i], tmp) #add a column with the file name
  out_cpg <- rbindlist(list(out_cpg, tmp), use.names= TRUE, fill=TRUE) #merge all tables together
}
#transform CpG percentages to numeric data
out_cpg <- out_cpg %>%
  mutate(across(-c(id), as.numeric))
out_cpg <- as.data.frame(out_cpg)
#view dataframe
View(out_cpg)
#number of probe site combinations
ncol(out_cpg)-1 #6172

#write dataframe to a file
write.csv(out_cpg, file= paste0(homewd,"/Dry-pipeline-WIP/Dataframes/sub_cpg_merged.csv"), row.names = F)


### 2. CHH sites ###
#find all of the *-sub_CHH.gz.bismark.cov files in the working directory
(bismark_files_chh <- list.files(pattern="*-sub_CHH.gz.bismark.cov", recursive=TRUE))

out_chh<-data.frame()
for(i in 1:length(bismark_files_chh)) {
  eg<-read.table(bismark_files_chh[i]) #read in each file in list as table
  tmp<-as.data.frame(t(eg)) #transpose table
  cols <- paste(tmp[1,], tmp[2,], sep="_") #create vector with probe site _ nucleotide #
  names(tmp) <- cols #make the vector the column names
  tmp <- tmp[-c(1:3,5:6),] #remove all rows except percent methylation data
  tmp <-  cbind("id" = bismark_files_chh[i], tmp) #add a column with the file name
  out_chh <- rbindlist(list(out_chh, tmp), use.names= TRUE, fill=TRUE) #merge all tables together
}
#transform CHH percentages to numeric data
out_chh <- out_chh %>%
  mutate(across(-c(id), as.numeric))
out_chh <- as.data.frame(out_chh)
#view dataframe
View(out_chh)
#number of probe site combinations
ncol(out_chh)-1 #28982

#write dataframe to a file
write.csv(out_chh, file= paste0(homewd,"/Dry-pipeline-WIP/Dataframes/sub_chh_merged.csv"), row.names = F)


### 3. CHG sites ###
#find all of the *-sub_CHG.gz.bismark.cov files in the working directory
(bismark_files_chg <- list.files(pattern="*-sub_CHG.gz.bismark.cov", recursive=TRUE))

out_chg<-data.frame()
for(i in 1:length(bismark_files_chg)) {
  eg<-read.table(bismark_files_chg[i]) #read in each file in list as table
  tmp<-as.data.frame(t(eg)) #transpose table
  cols <- paste(tmp[1,], tmp[2,], sep="_") #create vector with probe site _ nucleotide #
  names(tmp) <- cols #make the vector the column names
  tmp <- tmp[-c(1:3,5:6),] #remove all rows except percent methylation data
  tmp <-  cbind("id" = bismark_files_chg[i], tmp) #add a column with the file name
  out_chg <- rbindlist(list(out_chg, tmp), use.names= TRUE, fill=TRUE) #merge all tables together
}
#transform CHG percentages to numeric data
out_chg <- out_chg %>%
  mutate(across(-c(id), as.numeric))
out_chg <- as.data.frame(out_chg)
#view dataframe
View(out_chg)
#number of probe site combinations
ncol(out_chg)-1 #9524

#write dataframe to a file
write.csv(out_chg, file= paste0(homewd,"/Dry-pipeline-WIP/Dataframes/sub_chg_merged.csv"), row.names = F)


####################################################################################
# PART 2: GENERATE BREADTH COVERAGE STATISTICS
# Suggested next steps from Manny; "Unlike in Wilkinson, we have nuance that makes 
# things trickier. My first stab would be to do a correlation graph between the three 
# methylation sites per probe to see how or if they differ. Then get the weighed 
# average of methylated sites per probe sequence to figure out "how often is this 
# sequence methylated" as a percentage, which is the ultimate readout of microarrays."

#Correlation plot ???

#CALCULATE:
#1. CpG sites
#Proportion of samples with methylation data (i.e., NOT NA) per probe site
out_cpg_cov <- out_cpg[-(1:8),]
#Proportion of samples with methylation data (i.e., NOT NA) per probe site
prop_cpg <- as.numeric(colSums(!is.na(out_cpg[,-1]))/nrow(out_cpg))
#Average percent methylation per probe site (removes NAs)
mean_cpg <- colMeans(out_cpg[,-1], na.rm=T)
#Append data to the dataframe
out_cpg_cov[nrow(out_cpg_cov)+1, ] <- c("proportion_coverage", prop_cpg)
out_cpg_cov[nrow(out_cpg_cov)+1, ] <- c("mean_percent_methylation", mean_cpg)
View(out_cpg_cov)

write.csv(out_cpg_cov, file= paste0(homewd,"/Dry-pipeline-WIP/Dataframes/sub_cpg_merged_cov.csv"), row.names = F)

#2. CHH sites
#Proportion of samples with methylation data (i.e., NOT NA) per probe site
out_chh_cov <- out_chh[-(1:8),]
#Proportion of samples with methylation data (i.e., NOT NA) per probe site
prop_chh <- as.numeric(colSums(!is.na(out_chh[,-1]))/nrow(out_chh))
#Average percent methylation per probe site (removes NAs)
mean_chh <- colMeans(out_chh[,-1], na.rm=T)
#Append data to the dataframe
out_chh_cov[nrow(out_chh_cov)+1, ] <- c("proportion_coverage", prop_chh)
out_chh_cov[nrow(out_chh_cov)+1, ] <- c("mean_percent_methylation", mean_chh)
View(out_chh_cov)

write.csv(out_chh_cov, file= paste0(homewd,"/Dry-pipeline-WIP/Dataframes/sub_chh_merged_cov.csv"), row.names = F)


#3. CHG sites
out_chg_cov <- out_chg[-(1:8),]
#Proportion of samples with methylation data (i.e., NOT NA) per probe site
prop_chg <- as.numeric(colSums(!is.na(out_chg[,-1]))/nrow(out_chg))
#Average percent methylation per probe site (removes NAs)
mean_chg <- colMeans(out_chg[,-1], na.rm=T)
#Append data to the dataframe
out_chg_cov[nrow(out_chg_cov)+1, ] <- c("proportion_coverage", prop_chg)
out_chg_cov[nrow(out_chg_cov)+1, ] <- c("mean_percent_methylation", mean_chg)
View(out_chg_cov)

write.csv(out_chg_cov, file= paste0(homewd,"/Dry-pipeline-WIP/Dataframes/sub_chg_merged_cov.csv"), row.names = F)


# PLOTS OF BREADTH COVERAGE STATISTICS
# CPG
out_cpg_cov_plot <- as.data.frame(t(out_cpg_cov)) %>%
  row_to_names(row_number=1) %>%
  mutate(across(everything(), as.numeric))

p1 <- ggplot(out_cpg_cov_plot, aes(proportion_coverage)) + geom_histogram(binwidth = 0.125) + 
  labs(title ="CPG coverage - 8 samples") + scale_x_continuous(breaks = c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) #number of bats
p1
p2 <- ggplot(out_cpg_cov_plot, aes(mean_percent_methylation)) + geom_histogram() + 
  labs(title ="CPG percent methylation - 8 samples") + scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))
p2

# CHG
out_chg_cov_plot <- as.data.frame(t(out_chg_cov)) %>%
  row_to_names(row_number=1) %>%
  mutate(across(everything(), as.numeric))

p3 <- ggplot(out_chg_cov_plot, aes(proportion_coverage)) + geom_histogram(binwidth = 0.125) + 
  labs(title ="CHG coverage - 8 samples") + scale_x_continuous(breaks = c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) #number of bats
p3
p4 <- ggplot(out_chg_cov_plot, aes(mean_percent_methylation)) + geom_histogram() + 
  labs(title ="CHG percent methylation - 8 samples") + scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))
p4

# CHH
out_chh_cov_plot <- as.data.frame(t(out_chh_cov)) %>%
  row_to_names(row_number=1) %>%
  mutate(across(everything(), as.numeric))

p5 <- ggplot(out_chh_cov_plot, aes(proportion_coverage)) + geom_histogram(binwidth = 0.125) + 
  labs(title ="CHH coverage - 8 samples") + scale_x_continuous(breaks = c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) #number of bats
p5
p6 <- ggplot(out_chh_cov_plot, aes(mean_percent_methylation)) + geom_histogram() + 
  labs(title ="CHH percent methylation - 8 samples") + scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))
p6

plot_grid(p1, p2, p3, p4, p5, p6, ncol = 2)

####################################################################################
# PART 3: CALCULATE DEPTH COVERAGE STATISTICS

# Goal: use fastq files to figure out the coverage depth
# for all sites, to be used for filtering





####################################################################################
# PART 4: FILTERING BASED ON BREADTH AND DEPTH COVERAGE STATISTICS

# Goal: only keep sites that map to at least 1 genome (this criteria may already be in place in bismark) 
# this might need to be 2? based on how training works # look into wilkinson training
# and that are above a read depth threshhold
# 

# Goal: after determining coverage filtering based on 
# statistics generated above, generate a list of 
# sites to keep. Then use this list as a reference
# to trim the original dataframes to only contain
# those sites and their methylation data.
# This will create the final list of sites to be used
# for clock training.




##################################################################################
# PART 5: Add in metadata

# Goal: to make this dataframe ready for clock building, we need to
# add back in associated metadata for each sample.
# This includes the age, method of aging, and tissue type

# Because we are troubleshooting this on the myotis samples which
# do not have known ages, I am just going to make up ages
# in order to practice making the clocks

cpg_clock <- out_cpg
chg_clock <- out_chg
chh_clock <- out_chh

cpg_clock$known_age <- c(1,3,4,5,5,2,1,6)
cpg_clock$id <- c("M_californicus", "M_evotis", "M_evotis", "M_lucifugus", "M_lucifugus", "M_lucifugus", "M_thysanodes", "M_volans") #How do we want these to be written? just spp name? (to be able to build spp clocks, might need to add 'indl' column)
cpg_clock$tissue <- c("wingpunch", "wingpunch", "wingpunch", "cellline", "cellline", "cellline", "cellline", "cellline") #these we do know
cpg_clock$age_method <- c("DNAm", "DNAm", "DMAm", "markrecap", "dentition", "dentition", "dentition", "dentition") #these are all DNAm technically, but adding in some dentition and mark recap to test
write.csv(cpg_clock, file= paste0(homewd,"/Dry-pipeline-WIP/Dataframes/cpg_clock.csv"), row.names = F)

chg_clock$known_age <- c(3,5,3,1,1,2,2,3)
chg_clock$id <- c("M_californicus", "M_evotis", "M_evotis", "M_lucifugus", "M_lucifugus", "M_lucifugus", "M_thysanodes", "M_volans") #How do we want these to be written? just spp name? (to be able to build spp clocks, might need to add 'indl' column)
chg_clock$tissue <- c("wingpunch", "wingpunch", "wingpunch", "cellline", "cellline", "cellline", "cellline", "cellline") #these we do know
chg_clock$age_method <- c("DNAm", "DNAm", "DMAm", "markrecap", "dentition", "dentition", "dentition", "dentition") #these are all DNAm technically, but adding in some dentition and mark recap to test
write.csv(chg_clock, file= paste0(homewd,"/Dry-pipeline-WIP/Dataframes/chg_clock.csv"), row.names = F)

chh_clock$known_age <- c(5,6,7,8,7,6,7,8)
chg_clock$id <- c("M_californicus", "M_evotis", "M_evotis", "M_lucifugus", "M_lucifugus", "M_lucifugus", "M_thysanodes", "M_volans") #How do we want these to be written? just spp name? (to be able to build spp clocks, might need to add 'indl' column)
chh_clock$tissue <- c("wingpunch", "wingpunch", "wingpunch", "cellline", "cellline", "cellline", "cellline", "cellline") #these we do know
chh_clock$age_method <- c("markrecap", "DNAm", "DMAm", "markrecap", "markrecap", "dentition", "dentition", "dentition") #these are all DNAm technically, but adding in some dentition and mark recap to test
write.csv(chh_clock, file= paste0(homewd,"/Dry-pipeline-WIP/Dataframes/chh_clock.csv"), row.names = F)

# need to figure out how to bind together into a giant dataframe

############
# BY NOW, WE SHOULD HAVE DATAFRAMES THAT ARE READY TO BE USED IN CLOCK TRAINING

