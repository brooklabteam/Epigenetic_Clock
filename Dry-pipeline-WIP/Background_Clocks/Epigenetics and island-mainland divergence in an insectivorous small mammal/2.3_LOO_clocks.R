##R Script for Epigenetic Clock validation (leave-one-out)
#Main author: Joseph A. Zoller (from Steve Horvath's lab) 

#some of these packages might not be necessary
options(stringAsFactors=F)
library(tidyverse)
library(glmnet)
library(RColorBrewer)
library(devtools)
library(dplyr)

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

#load probe mapping data
#probe_mappability_table <- read_csv(probe_mappability_table.csv, col_types = cols(.default="c"))
#rm(probe_mappability_table.csv)

#load my probe alignment (here is your specific species aligned probe file)
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

#selecting samples that can not be used (NBU)
#we get an empty dataset by selecting nothing with the ! operator because all samples can be used (CBU)
#from N90 dataset select ones that cannot be used 
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

###############################################################################

#####USING ALIGNED PROBE DATA####
###Clock pre-analysis, using subsetted CpGs without outlier samples 
#setting up training groups for analysis
### Pre-Analysis, using subsetted clocks
set.seed(12345)
yxs.list <- alignInfoDat(infoCBUcinereusshrewn90,datAllSamp_subCPGcinereusshrewn90,"Basename","Basename")
ys <- yxs.list[[1]]
xs <- yxs.list[[2]]
rm(yxs.list)

### Analysis, using subsetted clocks
#no transformation for age
OUTVAR="Age"
out.rdata='Clock/SCsubs_LOOclock_EpigeneticAge.RData'
output.csv='Clock/SCsubs_LOOclock_EpigeneticAge_PredictedValues.csv'
out.png='Clock/SCsubs_LOOclock_EpigeneticAge.png'
out.png.title='Clock/SCsubs_LOOclock_EpigeneticAge'
PREDVAR="DNAmAgeLOO"
RESVAR="AgeAccelLOO"
fun_VAR1="averagedMaturity.yrs"
fun_VAR2="maxAgeCaesar"
COLVAR="Tissue"
ALPHA=0.5
NFOLD=10
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,COLVAR=COLVAR)
ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", PREDVAR, RESVAR, "cpg_count_used")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", PREDVAR, RESVAR, "cpg_count_used")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

#square root of age
OUTVAR="Age"
out.rdata='Clock/SCsubs_LOOclock_EpigeneticSqrtAge.RData'
output.csv='Clock/SCsubs_LOOclock_EpigeneticSqrtAge_PredictedValues.csv'
out.png='Clock/SCsubs_LOOclock_EpigeneticSqrtAge.png'
out.png.title='Clock/SCsubs_LOOclock_EpigeneticSqrtAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_sqrt_trans,fun_inv=fun_sqrt_inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", PREDVAR, RESVAR, "cpg_count_used")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", PREDVAR, RESVAR, "cpg_count_used")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

#log of age
OUTVAR="Age"
out.rdata='Clock/SCsubs_LOOclock_EpigeneticLogAge.RData'
output.csv='Clock/SCsubs_LOOclock_EpigeneticLogAge_PredictedValues.csv'
out.png='Clock/SCsubs_LOOclock_EpigeneticLogAge.png'
out.png.title='Clock/SCsubs_LOOclock_EpigeneticLogAge'
ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_log_trans,fun_inv=fun_log_inv,COLVAR=COLVAR)
ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", PREDVAR, RESVAR, "cpg_count_used")],
                         "Basename", all=T, sort=F)
ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", PREDVAR, RESVAR, "cpg_count_used")],
                         "Basename", all=T, sort=F)
write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)

#######other type of clocks that weren't used

## Log+Linear Transformed
# OUTVAR="Age"
# out.rdata='SpeciesSubsetAnalyses/Subset_CinereusShrewN90_AgeLOO_Final_Analysis/Subset_CinereusShrewN90_LOO_Final_subCPGcinereusshrewn90_EpigeneticLLinAge.RData'
# output.csv='SpeciesSubsetAnalyses/Subset_CinereusShrewN90_AgeLOO_Final_Analysis/Subset_CinereusShrewN90_LOO_Final_subCPGcinereusshrewn90_EpigeneticLLinAge_PredictedValues.csv'
# out.png='SpeciesSubsetAnalyses/Subset_CinereusShrewN90_AgeLOO_Final_Analysis/Subset_CinereusShrewN90_LOO_Final_subCPGcinereusshrewn90_EpigeneticLLinAge.png'
# out.png.title='Subset_CinereusShrewN90_LOO_Final_subCPGcinereusshrewn90_EpigeneticLLinAge'
# ys.output <- saveLOOEstimation(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_llin_trans,fun_inv=fun_llin_inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,COLVAR=COLVAR)
# ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", PREDVAR, RESVAR, "cpg_count_used")],
#                          "Basename", all=T, sort=F)
# ys.output <- base::merge(infoAllcinereusshrewn90, ys.output[, c("Basename", PREDVAR, RESVAR, "cpg_count_used")],
#                          "Basename", all=T, sort=F)
# write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)


### Plotting separately by tissue
OUTVAR="Age"
PREDVAR="DNAmAgeLOO"
PANELVAR="Tissue"
input.info_pred <- dplyr::filter(read.csv('Clock/SCsubs_LOOclock_EpigeneticAge_PredictedValues.csv', as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::select(Basename, Age, DNAmAgeLOO, Tissue) %>%
  dplyr::mutate(Tissue = factor(Tissue))
out.png='Clock/SCsubs_LOOclock_EpigeneticAge_TISSUE_PANEL.png'
out.png.title='Leave-One-Out Analysis by Tissue'
saveValidationPanelPlot(input.info_pred,OUTVAR,PREDVAR,PANELVAR,out.png,paste0(out.png.title,'\n'),mfrow=c(2,2),width=9,height=10)
output.csv='Clock/SCsubs_LOOclock_EpigeneticAge_TISSUE_PredictedValues.csv'
input.info_pred <- dplyr::filter(read.csv('Clock/SCsubs_LOOclock_EpigeneticAge_PredictedValues.csv', as.is=T), !is.na(DNAmAgeLOO)) %>%
  dplyr::select(Basename, Age, DNAmAgeLOO, Tissue) %>%
  dplyr::mutate(Tissue = factor(Tissue))
out.png='Clock/SCsubs_LOOclock_EpigeneticAge_TISSUE_PANEL.png'
out.png.title='Leave-One-Out Analysis by Tissue'
saveValidationPanelPlot(input.info_pred,OUTVAR,PREDVAR,PANELVAR,out.png,paste0(out.png.title,'\n'),mfrow=c(2,2),width=9,height=10)



