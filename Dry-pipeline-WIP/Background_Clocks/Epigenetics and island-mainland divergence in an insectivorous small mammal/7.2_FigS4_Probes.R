##R Script for supplemental Figure S4 (probe location and distance)
#Main author: Marie-Laurence Cossette

#set working directory
setwd("~/Desktop/SCHOOL/THESIS/Chapter 2 - Methyl/PAPER/Gitlab uploads/")

#load packages
library(dplyr)
library(ggVennDiagram)

#load annotation data from 1.1_Probe_alignment script
#this dataset has information on the CpGs and where they are located
Anno_SC <- read.csv("Anno_sorexcinereus.csv")
#remove tarseq_ to only have scaffold number in column
Anno_SC<- Anno_SC%>% mutate(seqnames = as.numeric(gsub("tarseq_", "", seqnames)))

#load shrew functional annotation data (this comes from annotating the genome using GenSAS)
#this dataset contains information on where genes are located on the genome
functional_ann <- read.table("functional annotation manh.txt")
#separating one column into multiple to get final column with only gene name to put on plot
functional_ann <- within(functional_ann, V6<-data.frame(do.call('rbind', strsplit(as.character(V6), '|', fixed=TRUE))))
functional_ann <- within(functional_ann, V6$X3<-data.frame(do.call('rbind', strsplit(as.character(V6$X3), '_', fixed=TRUE))))
functional_ann$gene <- functional_ann$V6$X3$X1
functional_ann$transcriptId <- functional_ann$V1

#join CpG annotation data set with functional annotation (gene) data set
Anno_SC  <- Anno_SC %>% full_join(functional_ann, by = c("transcriptId"))


#load EWAS data from the 4.1_EWAS_Island_effect and 4.2_BPI_effect scripts
#load readr to be able to import csv files
library(readr)
tab_BPI <- read_csv("EWAS/EbFit_BPI_effect.csv")
tab_Island <- read_csv("EWAS/EbFit_Island_effect.csv")
#rename first column
tab_BPI$CGid <- tab_BPI$...1
tab_Island$CGid <- tab_Island$...1

#merge each EWAS dataset with the annotations based on CGid column
allData_BPI <- merge(tab_BPI,Anno_SC, by = "CGid")
allData_Island <- merge(tab_Island,Anno_SC, by = "CGid")
#for background probes you can use any EWAS dataset as long as you don't filter out non-significant CpGs in the next step
allData_bg <- merge(tab_Island,Anno_SC, by = "CGid")

#make subset data frames for top hypo and hyper methylated sites for Island and BPI EWAS
#only keep sites with p value < 2e-6 
significant_BPI <- allData_BPI %>% filter(p.value.Location_1yes < 2e-6  )
significant_Island <- allData_Island %>% filter(p.value.Location_1yes <  2e-6   )

#add column with pos or neg methylation based on t.value for each data frame
#if < 0 then negative if > 0 then positive
significant_BPI <- significant_BPI %>%
  mutate(direction = case_when(significant_BPI$t.Location_1yes < 0 ~ "neg", significant_BPI$t.Location_1yes > 0 ~ "pos"))
BPI_neg <- significant_BPI %>% filter(direction == "neg")
BPI_pos <- significant_BPI %>% filter(direction == "pos")
#same as above for Island data frame
significant_Island <- significant_Island %>%
  mutate(direction = case_when(significant_Island$t.Location_1yes < 0 ~ "neg", significant_Island$t.Location_1yes > 0 ~ "pos"))
Island_neg <- significant_Island %>% filter(direction == "neg")
Island_pos <- significant_Island %>% filter(direction == "pos")

#load packages
library(tidyr)
library(stringr)
#here we will create a new column named locations to say near what feature each probe is located (e.g., intron, exon, etc.)
#this is already present in the annotation column but with extra detail (e.g., Intron (Sc.00g106270.m01...))
#using str_extract and "(\\w+)" we basically only copy the first word into the new location column
#do this for each data frame
significant_BPI$locations <- str_extract(significant_BPI$annotation, "(\\w+)")
significant_Island$locations <- str_extract(significant_Island$annotation, "(\\w+)")
allData_bg$locations <- str_extract(allData_bg$annotation, "(\\w+)")

###Here we want to create barplots for the the location of the different probes (e.g., near introns, exons, promoter, etc.)

#load ggplot
library(ggplot2)
#plot significant Bon Portage Island probe locations
#make x axis as location and choose to have color fill based on if the probe is hyper or hypomethylated
BPI_location <- ggplot(significant_BPI, aes(locations, fill = direction)) +
  #use dodge position to avoid hyper and hypomethylated bars to be stacked, puts them side by side instead
  #make width a bit smaller than default
  geom_bar(position = "dodge", width = 0.85) + 
  #make y axis scale jump by 50 
  scale_y_continuous(breaks=c(0,50,100,150,200,250)) +
  #flip axes to have x axis on the y axis and vice versa
  #make y axis (now the new x axis) stop at 225 because we want same range between BPI and Islands location plots
  coord_flip(ylim = c(0, 225)) +
  #select colors (red and blue)
  scale_fill_manual(values=c("#68B3FA","#FF6567")) + 
  #update title
  ggtitle("Bon Portage Island") +
  #update axes labels
  labs(y="CpG count", x = "CpG location") +
  #choose classic theme
  theme_classic() +
  #angle the CpG count value text
  theme(axis.text.x = element_text(angle=45)) +
  #remove legend
  theme(legend.position = "none")
#view
BPI_location

#plot significant Island probe locations
#make x axis as location and choose to have color fill based on if the probe is hyper or hypomethylated
Island_location <- ggplot(significant_Island, aes(locations, fill = direction)) +
  #use dodge position to avoid hyper and hypomethylated bars to be stacked, puts them side by side instead
  #make width a bit smaller than default
  geom_bar(position = "dodge", width = 0.85) + 
  #make y axis scale jump by 50 
  scale_y_continuous(breaks=c(0,50,100,150,200,250)) +
  #flip axes to have x axis on the y axis and vice versa
  #make y axis (now the new x axis) stop at 225 because we want same range between BPI and Islands location plots
  coord_flip(ylim = c(0, 225)) +
  #select colors (red and blue)
  scale_fill_manual(values=c("#68B3FA","#FF6567")) +
  #update title
  ggtitle("Islands") +
  #update axes labels
  labs(y="CpG count", x = "") +
  #choose classic theme
  theme_classic() +
  #angle the CpG count value text
  theme(axis.text.x = element_text(angle=45)) +
  #remove legend
  theme(legend.position = "none") 
#view
Island_location

#plot significant background probe locations
#make x axis as location
Background_location <- ggplot(allData_bg, aes(locations))+
  #use dodge position to avoid hyper and hypomethylated bars to be stacked, puts them side by side instead
  #make width a bit smaller than default
  geom_bar(position = "dodge", width = 0.85) +
  #make y axis scale jump by 5000 
  scale_y_continuous(breaks=c(0,5000,10000)) +
  #make y axis (now the new x axis) stop at 12200
  coord_flip(ylim = c(0,12200)) +
  #update title
  ggtitle("Background") +
  #update axes labels
  labs(y="CpG count", x = "")+
  #choose classic theme
  theme_classic() +
  #angle the CpG count value text and change size to fit the longer number values
  theme(axis.text.x = element_text(angle=45, size = 7)) +
  #remove legend
  theme(legend.position = "none")
#view
Background_location

###Here we want to create histograms for the distribution of probe distances from transcriptional start sites (TSS)

#plot significant Bon Portage Island probe distances from TSS
#make x axis as distance to TSS and choose to have color fill based on where probe is located (e.g., near Intron, Exon, etc.)
BPI_distance <- ggplot(significant_BPI, aes(x=distanceToTSS, fill=locations)) + 
  #make histogram with 500 different columns (bins), black outline color for bin and outline line thickness of 0.3
  geom_histogram(bins = 500, color = "black", size=0.3) +
  #limit x axis from -200,000 to 200,000 bp
  coord_cartesian(xlim=c(-200000,200000)) +
  #select specific color values for each CpG category (gray scale)
  scale_fill_manual(limits = c("Promoter", "Intron", "Exon", "Downstream", "Distal"), values=c("black", "gray35", "gray50", "gray80", "gray95")) + 
  #set theme as minimal and remove legend
  theme_minimal() + theme(legend.position = "none") +
  #to avoid x axis values being in scientific notation (because large numbers), specify the accuracy to 0 decimal place)  
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) 
#view
BPI_distance

#plot significant Island probe distances from TSS
#make x axis as distance to TSS and choose to have color fill based on where probe is located (e.g., near Intron, Exon, etc.)
Island_distance <- ggplot(significant_Island, aes(x=distanceToTSS, fill=locations)) + 
  #make histogram with 500 different columns (bins), black outline color for bin and outline line thickness of 0.3
  geom_histogram(bins = 500, color = "black", size = 0.3) +
  #limit x axis from -200,000 to 200,000 bp
  coord_cartesian(xlim=c(-200000,200000)) +
  #select specific color values for each CpG category (gray scale)
  scale_fill_manual(limits = c("Promoter", "Intron", "Exon", "Downstream", "Distal"), values=c("black", "gray35", "gray50", "gray80", "gray95")) + 
  #set theme as minimal and remove legend
  theme_minimal() + theme(legend.position = "none") +
  #to avoid x axis values being in scientific notation (because large numbers), specify the accuracy to 0 decimal place)  
  scale_x_continuous(labels = scales::number_format(accuracy = 1)) 
#view
Island_distance 

#plot background probe distances from TSS
#make x axis as distance to TSS and choose to have color fill based on where probe is located (e.g., near Intron, Exon, etc.)
Background_distance <- ggplot(allData_bg, aes(x=distanceToTSS, fill=locations)) + 
  #make histogram with 500 different columns (bins), black outline color for bin and outline line thickness of 0.3
  geom_histogram(bins = 500, color = "black", size = 0.3) +
  #limit x axis from -200,000 to 200,000 bp
  coord_cartesian(xlim=c(-200000,200000)) +
  #select specific color values for each CpG category (gray scale)
  scale_fill_manual(limits = c("Promoter", "Intron", "Exon", "Downstream", "Distal"), values=c("black", "gray35", "gray50", "gray80", "gray95")) + 
  #set theme as minimal
  theme_minimal() +
  #to avoid x axis values being in scientific notation (because large numbers), specify the accuracy to 0 decimal place)  
  scale_x_continuous(labels = scales::number_format(accuracy = 1))
#view
Background_distance 

#load cowplot to be able to merge all plots together
library(cowplot)
#merge all 6 plots together and add labels to each, specify that there are three columns 
Final_plot <- plot_grid(BPI_location, Island_location, Background_location, BPI_distance, Island_distance, Background_distance, labels = c('A', 'B', 'C', 'D', 'E', 'F'), ncol=3)
Final_plot



