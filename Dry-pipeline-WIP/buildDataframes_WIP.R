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
homewd= "/Users/theresalaverty/Documents/R/R_repositories/Epigenetic_Clock"
#Folder where the bismark.cov files of interest are stored
setwd(paste0(homewd,"/methylation-aging/output/bismark_bedgraph/"))

#Load packages
#install.packages("data.table")
#install.packages("dplyr")
library(data.table)
library(dplyr)
####################################################################################

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

# Suggested next steps from Manny; "Unlike in Wilkinson, we have nuance that makes 
# things trickier. My first stab would be to do a correlation graph between the three 
# methylation sites per probe to see how or if they differ. Then get the weighed 
# average of methylated sites per probe sequence to figure out "how often is this 
# sequence methylated" as a percentage, which is the ultimate readout of microarrays."

#Correlation plot ???

#CALCULATE:
#1. CpG sites
#Proportion of samples with methylation data (i.e., NOT NA) per probe site
prop_cpg <- colSums(!is.na(out_cpg[,-1]))/nrow(out_cpg)
#Average percent methylation per probe site (removes NAs)
mean_cpg <- colMeans(out_cpg[,-1], na.rm=T)
#Append data to the dataframe
out_cpg[nrow(out_cpg) + 1,] = c("proportion_coverage",prop_cpg)
out_cpg[nrow(out_cpg) + 1,] = c("mean_percent_methylation",mean_cpg)
View(out_cpg)

#2. CHH sites
#Proportion of samples with methylation data (i.e., NOT NA) per probe site
prop_chh <- colSums(!is.na(out_chh[,-1]))/nrow(out_chh)
#Average percent methylation per probe site (removes NAs)
mean_chh <- colMeans(out_chh[,-1], na.rm=T)
#Append data to the dataframe
out_chh[nrow(out_chh) + 1,] = c("proportion_coverage",prop_chh)
out_chh[nrow(out_chh) + 1,] = c("mean_percent_methylation",mean_chh)
View(out_chh)

#3. CHG sites
#Proportion of samples with methylation data (i.e., NOT NA) per probe site
prop_chg <- colSums(!is.na(out_chg[,-1]))/nrow(out_chg)
#Average percent methylation per probe site (removes NAs)
mean_chg <- colMeans(out_chg[,-1], na.rm=T)
#Append data to the dataframe
out_chg[nrow(out_chg) + 1,] = c("proportion_coverage",prop_chg)
out_chg[nrow(out_chg) + 1,] = c("mean_percent_methylation",mean_chg)
View(out_chg)
