##R Script for Epigenetic Clocks (including log and sqrt)
#Main author: Joseph A. Zoller (from Steve Horvath's lab) - Marie-Laurence Cossette tweaked it a bit based on the code from below
#https://figshare.com/articles/online_resource/Epigenetic_clock_construction_and_cross-validation/13522913?file=25958888
#some parts where I'm not sure why we are coding for things might be useful when using a data set includes multiple species?

#some of these packages might not be necessary
library(tidyverse)
library(ggplot2)
library(ggfortify)
library(glmnet)
library(survival)
library(tab) #tabcoxph
library(grid) #gtable
library(gridExtra) #grid.arrange
library(rmcorr)

#set working directory where all files are in 
setwd("~/Desktop/SCHOOL/THESIS/Chapter 2 - Methyl/PAPER/Gitlab uploads")

#run 2.1_Functions.R before running the rest of this script

###datasets needed###
#1. species age info sent by Joseph Zoller
#anAgeUpdatedCaesarVersion38.csv
#2. shrew info from Steve Horvath
#SampleSheetAgeN90final.csv
#3. in my personal dropbox folder in normalized data
datAllSamp_tp.rdata=c('all_probes_sesame_normalized.Rdata')
#4. probes I aligned for Sorex cinereus
#Sorex_cinereus_aln_probes.xlsx

###############################################################################
### LOADING ALL DATA ###
###############################################################################
#load anAge data
anAgeUpdated=read.csv('anAgeUpdatedCaesarVersion38.csv', as.is=T)
anAgeUpdated <- anAgeUpdated %>% dplyr::filter(profiled)
anAgeUpdated$gestationYears <- anAgeUpdated$Gestation.Incubation..days./365

#load sample sheet data
infoN15=read.csv('SampleSheetAgeN90final.csv')
infoAllSamp=read.csv('SampleSheetAgeN90final.csv', as.is=T)
infoAllSamp <- infoAllSamp %>%
  dplyr::select(Basename,SpeciesLatinName,OriginalOrderInBatch,Age,ConfidenceInAgeEstimate,
                CanBeUsedForAgingStudies,Tissue,Female,SpeciesCommonName,ExternalSampleID,Folder)
if ("Female" %in% colnames(infoAllSamp)) {
  infoAllSamp$Female[which(is.na(infoAllSamp$Female))] <- "NA"
  infoAllSamp$Female <- factor(infoAllSamp$Female, levels=c(0,1,"NA"))
  levels(infoAllSamp$Female) <- c("Male","Female","NA")
}

#load reformatted DNA methylation data, link methylation data to the specific shrew (first col is shrew iD, rest are CPG sites)
load(datAllSamp_tp.rdata)
datAllSamp_tp=normalized_betas_sesame
rm(datAllSamp_tp.rdata, normalized_betas_sesame)
datAllSamp <- datAllSamp_tp %>% dplyr::select(-1) %>% t() %>% as.data.frame()
colnames(datAllSamp)=datAllSamp_tp$CGid
datAllSamp <- cbind(Basename=rownames(datAllSamp), datAllSamp)
rownames(datAllSamp)=seq_len(nrow(datAllSamp))
datAllSamp$Basename <- as.character(datAllSamp$Basename)
rm(datAllSamp_tp)

#load my probe alignment (here is the specific species aligned probe file)
probe_marie_cinereusshrewn90 <- read.csv("aln_2.csv")

###############################################################################
### APPENDING anAge DATA TO SAMPLE SHEET ###
#basically adding info on reproductive maturity, gestation time... etc to my shrew N90 table
###############################################################################
anAgeUpdated.temp <- anAgeUpdated %>% dplyr::filter(profiled) %>%
  dplyr::select(SpeciesLatinName,Female.maturity..days.,Male.maturity..days.,
                averagedMaturity.yrs,maxAgeCaesar,gestationYears)
infoAllSamp$idx <- 1:nrow(infoAllSamp)
infoAllSamp <- base::merge(infoAllSamp, anAgeUpdated.temp, by="SpeciesLatinName", all.x=T, sort=F)
infoAllSamp <- infoAllSamp[order(infoAllSamp$idx),]
infoAllSamp <- select(infoAllSamp, -idx)
rm(anAgeUpdated.temp)

###############################################################################
### REFINING DATA ###
###############################################################################
#only keep good data, remove one outlier fetus basically (CBU = can be used)
infoCBUAllSamp <- infoAllSamp %>% dplyr::filter(CanBeUsedForAgingStudies=="yes") %>%
  dplyr::filter(ConfidenceInAgeEstimate>=90) %>%
  dplyr::filter(!is.na(Age)) %>%
  dplyr::filter(Basename %in% datAllSamp$Basename)

###############################################################################
### Sample filtering and partitioning
# anAgeUpdated <- dplyr::select(anAgeUpdated,SpeciesLatinName,Female.maturity..days.,
#                               Male.maturity..days.,averagedMaturity.yrs,maxAgeCaesar)
#select shrews from the species age info sent by Joseph Zoller (only keep shrew row)
anAgeUpdatedcinereusshrewn90 <- anAgeUpdated %>% dplyr::filter(SpeciesLatinName %in% c("Sorex cinereus"))
#not sure what this is? I think has to do with looking in our N90.2021-9098CinereusShrewAaronShafer folder for Folder??
infoAllcinereusshrewn90 <- infoAllSamp %>%
  dplyr::filter(Folder %in% c("N90.2021-9098CinereusShrewAaronShafer"))
#basically rename the whole shrew N90 dataset, only keeping the one that are good to use (CBU = can be used, removed outliers)
infoCBUcinereusshrewn90 <- infoCBUAllSamp %>%
  dplyr::filter(Folder %in% c("N90.2021-9098CinereusShrewAaronShafer"))
#making latin, common name and tissue columns factors
infoCBUcinereusshrewn90$SpeciesLatinName <- factor(infoCBUcinereusshrewn90$SpeciesLatinName)
infoCBUcinereusshrewn90$SpeciesCommonName <- factor(infoCBUcinereusshrewn90$SpeciesCommonName)
infoCBUcinereusshrewn90$Tissue <- factor(infoCBUcinereusshrewn90$Tissue)

# results in an empty dataset because basically selecting nothing with the ! operator
#from N90 dataset select ones that can not be used (NBU)
infoNBUcinereusshrewn90 <- infoAllcinereusshrewn90 %>% dplyr::filter(!Basename %in% infoCBUcinereusshrewn90$Basename) %>%
  dplyr::filter(CanBeUsedForAgingStudies == "yes")
#change columns to factor
infoNBUcinereusshrewn90$SpeciesLatinName <- factor(infoNBUcinereusshrewn90$SpeciesLatinName)
infoNBUcinereusshrewn90$SpeciesCommonName <- factor(infoNBUcinereusshrewn90$SpeciesCommonName)
infoNBUcinereusshrewn90$Tissue <- factor(infoNBUcinereusshrewn90$Tissue)


### Data refinement of CpG Sites, based on shared probe mappings
#remove all columns that don't have matching cgID as column name (removes about 8 thousand cols)
datAllSamp_subCPGcinereusshrewn90 <- datAllSamp[,c(1,which(colnames(datAllSamp) %in% probe_marie_cinereusshrewn90$qname))]

#latin name of shrew
table(as.character(infoCBUcinereusshrewn90$Tissue))
latin2common_cinereusshrewn90 <- unique(dplyr::select(infoCBUcinereusshrewn90,SpeciesLatinName,SpeciesCommonName))

#####USING ALIGNED PROBE DATA####
###Clock pre-analysis, using subsetted CpGs without outlier samples, Training from All shrews 
#setting up training groups fro analysis
set.seed(1236)
yxs.train.list <- alignInfoDat(infoCBUcinereusshrewn90,datAllSamp_subCPGcinereusshrewn90,"Basename","Basename")
yxs.test.list <- alignInfoDat(infoNBUcinereusshrewn90,datAllSamp_subCPGcinereusshrewn90,"Basename","Basename")
ys.train <- yxs.train.list[[1]]
xs.train <- yxs.train.list[[2]]
ys.test <- yxs.test.list[[1]]
xs.test <- yxs.test.list[[2]]
rm(yxs.train.list,yxs.test.list)

###Clock analysis, using subset CpGs without outlier samples , Training from All
OUTVAR="Age"
out.csv='Clock/SCsub_ClockBasedOnAll_EpigeneticAge.csv'
output.csv='Clock/SCsub_ClockBasedOnAll_EpigeneticAge_PredictedValues.csv'
out.png='Clock/SCsub_ClockBasedOnAll_EpigeneticAge.png'
out.png.title='Clock/SCsub_ClockBasedOnAll_EpigeneticAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinTestVAR="AgeAccelinTestbasedOnAll"
ALPHA=0.5
NFOLD=10
#run analysis
ys.output <- saveNetTrained(xs.train,ys.train,xs.test,ys.test,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinTestVAR,ALPHA,NFOLD)
#add output data to infoAllcinereusshrewn90 dataset 
ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinTestVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinTestVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

## Square-Root Transformed
OUTVAR="Age"
out.csv='Clock/SCsub_ClockBasedOnAll_EpigeneticSqrtAge.csv'
output.csv='Clock/SCsub_ClockBasedOnAll_EpigeneticSqrtAge_PredictedValues.csv'
out.png='Clock/SCsub_ClockBasedOnAll_EpigeneticSqrtAge.png'
out.png.title='Clock/SCsub_ClockBasedOnAll_EpigeneticSqrtAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinTestVAR="AgeAccelinTestbasedOnAll"
ALPHA=0.5
NFOLD=10
ys.output <- saveNetTrained(xs.train,ys.train,xs.test,ys.test,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinTestVAR,ALPHA,NFOLD,fun_trans=fun_sqrt_trans,fun_inv=fun_sqrt_inv)
ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinTestVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinTestVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

## Logarithm Transformed
OUTVAR="Age"
out.csv='Clock/SCsub_ClockBasedOnAll_EpigeneticLogAge.csv'
output.csv='Clock/SCsub_ClockBasedOnAll_EpigeneticLogAge_PredictedValues.csv'
out.png='Clock/SCsub_ClockBasedOnAll_EpigeneticLogAge.png'
out.png.title='Clock/SCsub_ClockBasedOnAll_EpigeneticLogAge'
PREDVAR="DNAmAgebasedOnAll"
RESVAR="AgeAccelbasedOnAll"
RESinTestVAR="AgeAccelinTestbasedOnAll"
ALPHA=0.5
NFOLD=10
ys.output <- saveNetTrained(xs.train,ys.train,xs.test,ys.test,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinTestVAR,ALPHA,NFOLD,fun_trans=fun_log_trans,fun_inv=fun_log_inv)
ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinTestVAR)],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinTestVAR)],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

#######other type of clocks that weren't used
#code doesn't work for some reason?

## Log+Linear Transformed
# OUTVAR="Age"
# out.csv='SpeciesSubsetAnalyses/Subset_CinereusShrewN90_Clocks_Final/Subset_CinereusShrewN90_Clock_basedOnAll_EpigeneticLLinAge.csv'
# output.csv='SpeciesSubsetAnalyses/Subset_CinereusShrewN90_Clocks_Final/Subset_CinereusShrewN90_Clock_basedOnAll_EpigeneticLLinAge_PredictedValues.csv'
# out.png='SpeciesSubsetAnalyses/Subset_CinereusShrewN90_Clocks_Final/Subset_CinereusShrewN90_Clock_basedOnAll_EpigeneticLLinAge.png'
# out.png.title='Subset_CinereusShrewN90_Clock_basedOnAll_EpigeneticLLinAge'
# PREDVAR="DNAmAgebasedOnAll"
# RESVAR="AgeAccelbasedOnAll"
# RESinTestVAR="AgeAccelinTestbasedOnAll"
# fun_VAR1="averagedMaturity.yrs"
# fun_VAR2="maxAgeCaesar"
# ALPHA=0.5
# NFOLD=10
# ys.output <- saveNetTrained(xs.train,ys.train,xs.test,ys.test,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinTestVAR,ALPHA,NFOLD,fun_trans=fun_llin_trans,fun_inv=fun_llin_inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2)
# ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinTestVAR)],
#                          "Basename", all=T, sort=F)
# ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinTestVAR)],
#                          "Basename", all=T, sort=F)
# write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)


## Relative Age
##need to tweak to make it work but only necessary for comparing different species
#OUTVAR="RelAge"
#out.csv='Subset_CinereusShrewN90_Clock_subCPGcinereusshrewn90_basedOnAll_EpigeneticRelAge.csv'
#output.csv='Subset_CinereusShrewN90_Clock_subCPGcinereusshrewn90_basedOnAll_EpigeneticRelAge_PredictedValues.csv'
#out.png='Subset_CinereusShrewN90_Clock_subCPGcinereusshrewn90_basedOnAll_EpigeneticRelAge.png'
#out.png.title='Subset_CinereusShrewN90_Clock_subCPGcinereusshrewn90_basedOnAll_EpigeneticRelAge'
#PREDVAR="DNAmRelAgebasedOnAll"
#RESVAR="RelAgeAccelbasedOnAll"
#RESinTestVAR="RelAgeAccelinTestbasedOnAll"
#COLVAR="SpeciesLatinName"
#ALPHA=0.5
#NFOLD=10
#ys.output <- saveNetTrained(xs.train,ys.train,xs.test,ys.test,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinTestVAR,ALPHA,NFOLD)
#ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinTestVAR)],
#                         "Basename", all=T, sort=F)
#ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", "train", PREDVAR, RESVAR, RESinTestVAR)],
#                         "Basename", all=T, sort=F)
#write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
