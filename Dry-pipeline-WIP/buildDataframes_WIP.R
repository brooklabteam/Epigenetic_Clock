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
homewd= "/Users/sophiahorigan/Documents/GitHub/Epigenetic_Clock"
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
# PART 3: UNDERSTAND HOW METHYLATION VARIES BY cpXXXXXXXX

# Goal: look into how sites of the same cpXXXXXXXX relate to one another
# May be used for filtering

### 1. CpG sites ###
#transpose data frame so cpXXXXXXXX are rows and samples are columns
cp_cpg_all <- t(data.frame(c("cp_cpg",substr(colnames(out_cpg[,-1]), 1, 10)))) 
colnames(cp_cpg_all)<-colnames(out_cpg)
cp_cpg_all <- rbind(out_cpg, cp_cpg_all)
rownames(cp_cpg_all)<-NULL
cp_cpg_all <-t(cp_cpg_all)
colnames(cp_cpg_all) <- cp_cpg_all[1,]
cp_cpg_all <- cp_cpg_all[-1,]
cp_cpg_all <- as.data.frame(cp_cpg_all)
cp_cpg_all[,1:8] <- as.numeric(cp_cpg_all[,1:8]) 
View(cp_cpg_all)

#Calculate cell counts with data (excluding NAs) by cpXXXXXXXX
count_cp_cpg <- cp_cpg_all %>%
  group_by(cp_cpg) %>%
  select_if(function(x) any(is.na(x))) %>%
  summarise_all(funs(sum(!is.na(.)))) %>%
  mutate(num=rowSums(.[2:9]))
  
#Calculate overall cell counts (includes NAs) by cpXXXXXXXX
num_cells_cp_cpg <- cp_cpg_all %>% 
  count(cp_cpg, name = "no_rows", .drop = F) %>%
  mutate(total_cells = no_rows*(ncol(cp_cpg_all)-1))

#Sum all methylation percentages by cpXXXXXXXX
# make datafram without cp_cpg column
cp_cpg_all <- cp_cpg_all %>%
  mutate(across(-c(cp_cpg), as.numeric))
cp_cpg_all <- as.data.frame(cp_cpg_all)
#sum by cpXXXXXXXX
sum_cp_cpg <- aggregate(. ~ cp_cpg, cp_cpg_all,              
          function(x) sum(x, na.rm=TRUE), na.action = na.pass)
#bind sums with cp_cpg
sum_cp_cpg <- cbind(sum_cp_cpg$cp_cpg, as.numeric(rowSums(sum_cp_cpg[,2:9])))
sum_cp_cpg <- as.data.frame(sum_cp_cpg)
colnames(sum_cp_cpg) <- c("cp_cpg", "sum") #name columns
sum_cp_cpg$sum <- as.numeric(sum_cp_cpg$sum) #make numeric

#Merge cell counts with data, overall cell counts, and sum
together_cp_cpg <- merge(count_cp_cpg, num_cells_cp_cpg, by="cp_cpg")
together_cp_cpg <- merge(together_cp_cpg, sum_cp_cpg, by="cp_cpg")
together_cp_cpg <- together_cp_cpg[,-c(2:9,11)]

#Calculate proportion coverage and mean percent methylation by cpXXXXXXXX
prop_cpg_by_cp <- together_cp_cpg %>% 
  group_by(cp_cpg) %>%
  mutate(prop=num/total_cells, mean=sum/num) 
#prop = number of samples with data at that cpXXXXXXXX site divided by (no_rows per cpXXXXXXXX * no_samples)
#mean = the methylation percentage divided by the number of samples with data at that cpXXXXXXXX site

#Create a table with each cpXXXXXXXX and prop and mean
out_cpg_cov_by_cp <- prop_cpg_by_cp[,c(1, 5:6)]
colnames(out_cpg_cov_by_cp)[2:3] <- c("proportion_coverage", "mean_percent_methylation")
View(out_cpg_cov_by_cp)


### 2. CHH sites ###
#transpose data frame so cpXXXXXXXX are rows and samples are columns
cp_chh_all <- t(data.frame(c("cp_chh",substr(colnames(out_chh[,-1]), 1, 10)))) 
colnames(cp_chh_all)<-colnames(out_chh)
cp_chh_all <- rbind(out_chh, cp_chh_all)
rownames(cp_chh_all)<-NULL
cp_chh_all <-t(cp_chh_all)
colnames(cp_chh_all) <- cp_chh_all[1,]
cp_chh_all <- cp_chh_all[-1,]
cp_chh_all <- as.data.frame(cp_chh_all)
cp_chh_all[,1:8] <- as.numeric(unlist(cp_chh_all[,1:8])) 
View(cp_chh_all)

#Calculate cell counts with data (excluding NAs) by cpXXXXXXXX
count_cp_chh <- cp_chh_all %>%
  group_by(cp_chh) %>%
  select_if(function(x) any(is.na(x))) %>%
  summarise_all(funs(sum(!is.na(.)))) %>%
  mutate(num=rowSums(.[2:9]))

#Calculate overall cell counts (includes NAs) by cpXXXXXXXX
num_cells_cp_chh <- cp_chh_all %>% 
  count(cp_chh, name = "no_rows", .drop = F) %>%
  mutate(total_cells = no_rows*(ncol(cp_chh_all)-1))

#Sum all methylation percentages by cpXXXXXXXX
# make datafram without cp_chh column
cp_chh_all <- cp_chh_all %>%
  mutate(across(-c(cp_chh), as.numeric))
cp_chh_all <- as.data.frame(cp_chh_all)
#sum by cpXXXXXXXX
sum_cp_chh <- aggregate(. ~ cp_chh, cp_chh_all,              
                        function(x) sum(x, na.rm=TRUE), na.action = na.pass)
#bind sums with cp_chh
sum_cp_chh <- cbind(sum_cp_chh$cp_chh, as.numeric(rowSums(sum_cp_chh[,2:9])))
sum_cp_chh <- as.data.frame(sum_cp_chh)
colnames(sum_cp_chh) <- c("cp_chh", "sum") #name columns
sum_cp_chh$sum <- as.numeric(sum_cp_chh$sum) #make numeric

#Merge cell counts with data, overall cell counts, and sum
together_cp_chh <- merge(count_cp_chh, num_cells_cp_chh, by="cp_chh")
together_cp_chh <- merge(together_cp_chh, sum_cp_chh, by="cp_chh")
together_cp_chh <- together_cp_chh[,-c(2:9,11)]

#Calculate proportion coverage and mean percent methylation by cpXXXXXXXX
prop_chh_by_cp <- together_cp_chh %>% 
  group_by(cp_chh) %>%
  mutate(prop=num/total_cells, mean=sum/num) 
#prop = number of samples with data at that cpXXXXXXXX site divided by (no_rows per cpXXXXXXXX * no_samples)
#mean = the methylation percentage divided by the number of samples with data at that cpXXXXXXXX site

#Create a table with each cpXXXXXXXX and prop and mean
out_chh_cov_by_cp <- prop_chh_by_cp[,c(1, 5:6)]
colnames(out_chh_cov_by_cp)[2:3] <- c("proportion_coverage", "mean_percent_methylation")
View(out_chh_cov_by_cp)


### 3. CHG sites ###
#transpose data frame so cpXXXXXXXX are rows and samples are columns
cp_chg_all <- t(data.frame(c("cp_chg",substr(colnames(out_chg[,-1]), 1, 10)))) 
colnames(cp_chg_all)<-colnames(out_chg)
cp_chg_all <- rbind(out_chg, cp_chg_all)
rownames(cp_chg_all)<-NULL
cp_chg_all <-t(cp_chg_all)
colnames(cp_chg_all) <- cp_chg_all[1,]
cp_chg_all <- cp_chg_all[-1,]
cp_chg_all <- as.data.frame(cp_chg_all)
cp_chg_all[,1:8] <- as.numeric(unlist(cp_chg_all[,1:8])) 
View(cp_chg_all)

#Calculate cell counts with data (excluding NAs) by cpXXXXXXXX
count_cp_chg <- cp_chg_all %>%
  group_by(cp_chg) %>%
  select_if(function(x) any(is.na(x))) %>%
  summarise_all(funs(sum(!is.na(.)))) %>%
  mutate(num=rowSums(.[2:9]))

#Calculate overall cell counts (includes NAs) by cpXXXXXXXX
num_cells_cp_chg <- cp_chg_all %>% 
  count(cp_chg, name = "no_rows", .drop = F) %>%
  mutate(total_cells = no_rows*(ncol(cp_chg_all)-1))

#Sum all methylation percentages by cpXXXXXXXX
# make datafram without cp_chg column
cp_chg_all <- cp_chg_all %>%
  mutate(across(-c(cp_chg), as.numeric))
cp_chg_all <- as.data.frame(cp_chg_all)
#sum by cpXXXXXXXX
sum_cp_chg <- aggregate(. ~ cp_chg, cp_chg_all,              
                        function(x) sum(x, na.rm=TRUE), na.action = na.pass)
#bind sums with cp_chg
sum_cp_chg <- cbind(sum_cp_chg$cp_chg, as.numeric(rowSums(sum_cp_chg[,2:9])))
sum_cp_chg <- as.data.frame(sum_cp_chg)
colnames(sum_cp_chg) <- c("cp_chg", "sum") #name columns
sum_cp_chg$sum <- as.numeric(sum_cp_chg$sum) #make numeric

#Merge cell counts with data, overall cell counts, and sum
together_cp_chg <- merge(count_cp_chg, num_cells_cp_chg, by="cp_chg")
together_cp_chg <- merge(together_cp_chg, sum_cp_chg, by="cp_chg")
together_cp_chg <- together_cp_chg[,-c(2:9,11)]

#Calculate proportion coverage and mean percent methylation by cpXXXXXXXX
prop_chg_by_cp <- together_cp_chg %>% 
  group_by(cp_chg) %>%
  mutate(prop=num/total_cells, mean=sum/num) 
#prop = number of samples with data at that cpXXXXXXXX site divided by (no_rows per cpXXXXXXXX * no_samples)
#mean = the methylation percentage divided by the number of samples with data at that cpXXXXXXXX site

#Create a table with each cpXXXXXXXX and prop and mean
out_chg_cov_by_cp <- prop_chg_by_cp[,c(1, 5:6)]
colnames(out_chg_cov_by_cp)[2:3] <- c("proportion_coverage", "mean_percent_methylation")
View(out_chg_cov_by_cp)


# PLOTS OF BREADTH COVERAGE STATISTICS by cpXXXXXXXX
# CPG
p1_by_cp <- ggplot(out_cpg_cov_by_cp, aes(proportion_coverage)) + geom_histogram(binwidth = 0.125) + 
  labs(title ="CPG coverage by cpXXXXXXXX - 8 samples") + scale_x_continuous(breaks = c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) #number of bats
#p1_by_cp
p2_by_cp <- ggplot(out_cpg_cov_by_cp, aes(mean_percent_methylation)) + geom_histogram() + 
  labs(title ="CPG percent methylation by cpXXXXXXXX - 8 samples") + scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))
#p2_by_cp

# CHG
p3_by_cp <- ggplot(out_chh_cov_by_cp, aes(proportion_coverage)) + geom_histogram(binwidth = 0.125) + 
  labs(title ="CHG coverage by cpXXXXXXXX - 8 samples") + scale_x_continuous(breaks = c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) #number of bats
#p3_by_cp
p4_by_cp <- ggplot(out_chh_cov_by_cp, aes(mean_percent_methylation)) + geom_histogram() + 
  labs(title ="CHG percent methylation by cpXXXXXXXX - 8 samples") + scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))
#p4_by_cp

# CHH
p5_by_cp <- ggplot(out_chg_cov_by_cp, aes(proportion_coverage)) + geom_histogram(binwidth = 0.125) + 
  labs(title ="CHH coverage by cpXXXXXXXX - 8 samples") + scale_x_continuous(breaks = c(0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1)) #number of bats
#p5_by_cp
p6_by_cp <- ggplot(out_chg_cov_by_cp, aes(mean_percent_methylation)) + geom_histogram() + 
  labs(title ="CHH percent methylation by cpXXXXXXXX - 8 samples") + scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))
#p6_by_cp

plot_grid(p1_by_cp, p2_by_cp, p3_by_cp, p4_by_cp, p5_by_cp, p6_by_cp, ncol = 2)

#Look at overlap of cpXXXXXXXX sites across CpG, CHH, CHG sites
colnames(out_cpg_cov_by_cp)[1:3] <- c("cp", "cpg_proportion_coverage", "cpg_mean_percent_methylation")
colnames(out_chh_cov_by_cp)[1:3] <- c("cp", "chh_proportion_coverage", "chh_mean_percent_methylation")
colnames(out_chg_cov_by_cp)[1:3] <- c("cp", "chg_proportion_coverage", "chg_mean_percent_methylation")
#determine how long merged df should be
length(out_cpg_cov_by_cp$cp)#1664
length(out_chh_cov_by_cp$cp)#1873
length(out_chg_cov_by_cp$cp)#1757
#put all data frames into list
df_list <- list(out_cpg_cov_by_cp, out_chh_cov_by_cp, out_chg_cov_by_cp)
#merge all data frames in list
(merged_cov_by_cp <- df_list %>% reduce(full_join, by='cp'))
length(merged_cov_by_cp$cp) #1880
sum(complete.cases(merged_cov_by_cp)) #1588 rows have data for CpG, CHH, and CHG
output <- sapply(merged_cov_by_cp,is.na) 
table(rowSums(output)/2) 
#0 = no. rows in which CpG, CHH, and CHG have data, 1 = no. rows in which two have data, 2 = no. rows in which only one has data
# 0    1    2 
# 1588  238   54 

####################################################################################
# PART 4: CALCULATE DEPTH COVERAGE STATISTICS

# Goal: summarize depth statistics from bismark files
(bismark_files_cpg <- list.files(pattern="*_CpG.gz.bismark.cov", recursive=TRUE))

out_cpg_depth<-data.frame()
for(i in 1:length(bismark_files_cpg)) {
  eg<-read.table(bismark_files_cpg[i]) #read in each file in list as table
  tmp<-as.data.frame(eg)
  tmp$cov<-tmp[,5]+tmp[,6]
  #tmp <- tmp[-c(1:3,5:6),] #remove all rows except percent methylation data
  tmp <-  cbind("id" = bismark_files_cpg[i], tmp) #add a column with the file name
  out_cpg_depth <- rbindlist(list(out_cpg_depth, tmp), use.names= TRUE, fill=TRUE) #merge all tables together
}
write.csv(out_cpg_depth, file= paste0(homewd,"/Dry-pipeline-WIP/Dataframes/cpg_depth.csv"), row.names = F)



####################################################################################
# PART 5: FILTERING BASED ON BREADTH AND DEPTH COVERAGE STATISTICS

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
# PART 6: Add in metadata

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

