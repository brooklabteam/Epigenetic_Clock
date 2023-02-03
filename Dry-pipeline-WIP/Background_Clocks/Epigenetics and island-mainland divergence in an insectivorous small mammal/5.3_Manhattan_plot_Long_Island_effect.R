## R Script for EWAS Manhattan plot (Long Island vs other populations)
#Main author: Marie-Laurence Cossette (inspired by plots from Steve Horvath's lab and below links)
#https://danielroelfs.com/blog/how-i-create-manhattan-plots-using-ggplot/
#https://www.r-graph-gallery.com/101_Manhattan_plot.html

#packages to load
library(dplyr)
library(ggplot2)
library(gghighlight)
library(ggrepel)
library(readr)

#set working directory desktop
setwd("~/Desktop/SCHOOL/THESIS/Chapter 2 - Methyl/PAPER/Gitlab uploads/")

#load annotation data from 1.1_Probe_alignment script
Anno_SC <- read.csv("Anno_sorexcinereus.csv")
#remove tarseq_ to only have scaffold number in column
Anno_SC<- Anno_SC%>% mutate(seqnames = as.numeric(gsub("tarseq_", "", seqnames)))
 
#load EWAS data from 4.3_EWAS_Long_Island_effect script
tab <- read_csv("EWAS/EbFit_LongIsland_effect.csv")

#merge EWAS and annotations based on CGid
allData <- merge(tab,Anno_SC, by = "CGid", all=T)

#make numeric variables for base pair (where probe starts) and chromosome (scaffold name)
allData$bp <- as.numeric(allData$probeStart)
allData$chr<- as.numeric(allData$seqnames)

#make column with cumulative base pairs
#selects the largest position for each chromosome, and then calculates the cumulative sum of those
data_cum <- allData %>% 
  group_by(chr) %>% 
  summarise(max_bp = max(bp)) %>% 
  mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% 
  select(chr, bp_add)

#add to main dataset the new column
allData_manhattan <- allData %>% 
  inner_join(data_cum, by = "chr") %>% 
  mutate(bp_cum = bp + bp_add)

#rename p value column
allData_manhattan$p <- allData_manhattan$p.value.Location_1yes
#based on t value ( - or +) make new column indicating hyper or hypomethylated
allData_manhattan <- allData_manhattan %>%
  mutate(direction = case_when(allData_manhattan$t.Location_1yes < 0 ~ "neg", allData_manhattan$t.Location_1yes > 0 ~ "pos"))
#save to computer
write.csv(allData_manhattan,file="allData_manhattan_sorexcinereus_LongIsland.csv")

#to set axis ticks to be in middle of each chromosome
axis_set <- allData_manhattan %>% 
  group_by(chr) %>% 
  summarize(center = mean(bp_cum))

#to set y-axis based on your p values
#transform the largest exponent to positive and add 2, to give some extra space on the top edge of the plot
ylim <- allData_manhattan %>% 
  filter(allData_manhattan$p == min(allData_manhattan$p)) %>% 
  mutate(ylim = abs(floor(log10(p))) + 2) %>% 
  pull(ylim)

#indicate the significance threshold
#bonferroni (0.05/24460 total obs)
sig <- 2e-6

#make subsetted data frames for top hypo and hypermethylated sites (for plotting different colors)
significant<- allData_manhattan %>% filter(p < 2e-6)
significant_pos <- subset(significant, direction =="pos")
significant_neg <- subset(significant, direction =="neg")

#actual gene names are in another table (from GenSAS), so load it up - functional annotations
functional_ann <- read.table("functional annotation manh.txt")
#separating one column into multiple to get final column with only gene name to put on plot
functional_ann <- within(functional_ann, V6<-data.frame(do.call('rbind', strsplit(as.character(V6), '|', fixed=TRUE))))
functional_ann <- within(functional_ann, V6$X3<-data.frame(do.call('rbind', strsplit(as.character(V6$X3), '_', fixed=TRUE))))
functional_ann$gene <- functional_ann$V6$X3$X1
functional_ann$transcriptId <- functional_ann$V1

#join data frame with new one containing functional annotations
allData_manhattan  <- allData_manhattan %>% full_join(functional_ann, by = c("transcriptId"))

#ggplot for Manhattan plot
#x axis is location of CpGs on genome, y is the log of p value, want each scaffold to have alternating color
manhplot <- ggplot(allData_manhattan, aes(x = bp_cum, y = -log10(p))) +
  #line for significant p value threshold (2e-6), make it red and dashed
  geom_hline(yintercept = -log10(sig), color = "#f32712", linetype = "dashed") +
  #color each point based on scaffold and change size
  geom_point(aes(color=as.factor(chr)), size=0.8) +
  #point color for significant + (red) and - (blue) methylated genes
  #make bigger than non significant points
  geom_point(data = significant_pos, colour = "#fb6869" ,size=1.2) +
  geom_point(data = significant_neg, colour = "#66b1ff",size=1.2) +
  #annotation info for top n (want 15 but many don't have gene annotations so select 27 to get 15)
  #add label for gene name, make font size 5 and black
  geom_text_repel(data=allData_manhattan %>% top_n(-27, p), aes(label=gene), size=5, color ="black") +
  #place ticks in middle of each chromosome/scaffold on x axis
  scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
  #y axis limits set and breaks based on p values
  #start at 0 till 24 but labels start at 4 till 20 with jumps of 4
  scale_y_continuous(expand = c(0,0), limits = c(0, ylim), breaks = seq(4, 24, by = 4)) +
  #select colors for non significant plot points
  #make the colors repeat over whole genome
  scale_color_manual(values = rep(c("#bebebe", "#f8e9d5"), unique(length(axis_set$chr)))) +
  #axis labels
  labs(x = "Chromosome",y = "-log10(P)") + 
  #title 
  ggtitle("EWAS of Long Island Effect") +
  #setting up background lines and white background color 
  theme(panel.background = element_rect(fill = "white", colour = "white",size = 0.5, linetype = "solid"),
  panel.grid.major = element_line(size = 0.25, linetype = 'solid',colour = "#ececec"), 
  panel.grid.minor = element_line(size = 0.15, linetype = 'solid',colour = "#ececec"),
  #remove legend
  legend.position = "none",
  #remove vertical lines and chromosome numbers, basically just make the plot cleaner
  panel.border = element_blank(),
  panel.grid.major.x = element_blank(),
  panel.grid.minor.x = element_blank(),
  axis.text.x = element_blank(),
  #position title in the middle and select size and bold
  plot.title = element_text(hjust = 0.5, face="bold")
  )
#plot
plot(manhplot)

#view data for top 30 significant (p value) CpGs
#make new column with what species the gene that match came from
allData_manhattan$species <- allData_manhattan$V6$X3$X2
#select top 30 based on p value
top <- allData_manhattan %>% top_n(-30, p)
#only subset columns of interest
top_clean = subset(top, select = c("CGid", "geneId","gene", "direction", "p", "species", "V5", "annotation") )

