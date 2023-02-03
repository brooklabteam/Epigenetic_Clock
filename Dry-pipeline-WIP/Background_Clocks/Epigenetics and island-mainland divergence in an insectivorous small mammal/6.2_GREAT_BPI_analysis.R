## R Script for GREAT analysis (to find gene pathways) for BPI vs other populations
#Main author: Marie-Laurence Cossette (inspired by plots and code from Steve Horvath's lab, Amin Haghani's work and other mammalian methylation papers)

#set working directory desktop
setwd("~/Desktop/SCHOOL/THESIS/Chapter 2 - Methyl/PAPER/Gitlab uploads/")

#load package
library(dplyr)
#load annotation data from 1.1_Probe_alignment script
Anno_SC <- read.csv("Anno_sorexcinereus.csv")

#actual gene names are in another table (from GenSAS), so load it up - functional annotations
functional_ann <- read.table("functional annotation manh.txt")
#separating one column into multiple to get final column with only gene name to put on plot
functional_ann <- within(functional_ann, V6<-data.frame(do.call('rbind', strsplit(as.character(V6), '|', fixed=TRUE))))
functional_ann <- within(functional_ann, V6$X3<-data.frame(do.call('rbind', strsplit(as.character(V6$X3), '_', fixed=TRUE))))
functional_ann$gene <- functional_ann$V6$X3$X1
functional_ann$transcriptId <- functional_ann$V1

#join anotation data with functional annotations
Anno_SC <- Anno_SC %>% left_join(functional_ann, by = c("transcriptId"))
#remove CpGs that don't have gene information and make new dataset
Anno <- Anno_SC[!is.na(Anno_SC$gene), ]  

#load EWAS data from 4.2_EWAS_BPI_effect script
tab <- read.csv('EWAS/EbFit_BPI_effect.csv')

#set names for output files
options(stringsAsFactors = F)
out.txt= 'HUM_background_location.txt'
out1.txt='HUM_pos_top500_location.txt'
out2.txt='HUM_neg_top500_location.txt'
out3.txt='HUM_pos_permu_top500_location.txt'
out4.txt='HUM_neg_permu_top500_location.txt'

#use human annotation dataset by Amin Haghani
#used V7 for my initial analysis but now only V9 accessible
anno.RDS=c('Human.Homo_sapiens.hg19.Amin.V9.RDS')
myanno_HUM=readRDS(anno.RDS)

#make data set with only one column containing CpG names
Anno <- Anno %>% select("CGid")
#merge with human annotations, basically make dataset with only our relevant CpGs for shrews but 
Anno_final <- Anno %>% inner_join(myanno_HUM, by = c("CGid"))

#make input = tab dataset (which is the EWAS data)
input=tab
#merge EWAS dataset with human annotations
input.all=merge(by='CGid',input,Anno_final)
#look at how many times each chromosome appears in columns
#in other words, which ones have more CpGs included in this analysis
table(input.all$geneChr,useNA='ifany')
table(input.all$seqnames,useNA='ifany')
#remove any rows that have NAs
input.all=subset(input.all,!is.na(geneChr))
#make new column named CHR with seqnames (basically just chromomsome numbers)
input.all$CHR=as.character(input.all$seqnames)

#order p values
input.all=input.all[order(input.all$p.value.Location_1yes),]
#subset hypo and hypermethylated CpGs
pos=subset(input.all,t.Location_1yes>0)
neg=subset(input.all,t.Location_1yes<0)
#export tables with all the CpG's info and for top hypo and hypermethylated CpGs too
#we are interested in where the CpG starts and ends on the genome
write.table(subset(input.all, select=c(CHR,CGstart,CGend,CGid)),out.txt,sep=' ',row.names = F,quote=F,col.names = F)
write.table(subset(pos[1:500,],select=c(CHR,CGstart,CGend,CGid)),out1.txt,sep=' ',row.names = F,quote=F,col.names = F)
write.table(subset(neg[1:500,],select=c(CHR,CGstart,CGend,CGid)),out2.txt,sep=' ',row.names = F,quote=F,col.names = F)

#set seed so that it's consistent every time you run script
set.seed(12345)
#add row to pos dataset (hypermethylated CpGs) with random order for each CpG
#'sample' randomly selects data from the dataset
#'nrow' function returns the number of rows in the dataset
pos$permu=sample(nrow(pos))
#order the permu column
pos=pos[order(pos$permu),]
#make table with 500 permutations for hypermethylated CpGs
write.table(subset(pos[1:500,],select=c(CHR,CGstart,CGend,CGid)),out3.txt,sep=' ',row.names = F,quote=F,col.names = F)
#do same as above but for hypomethylated CpGs
neg$permu=sample(nrow(neg))
neg=neg[order(neg$permu),]
write.table(subset(neg[1:500,],select=c(CHR,CGstart,CGend,CGid)),out4.txt,sep=' ',row.names = F,quote=F,col.names = F)

#use GREAT analysis version 3.0 with default setting and use a background based on Mammalian array
#link for more details on how to use GREAT
#https://great-help.atlassian.net/wiki/spaces/GREAT/overview?homepageId=655450

#load package
library(rGREAT)

#set input file names
background.txt='HUM_background_location.txt'
input.pos.txt='HUM_pos_top500_location.txt'
input.neg.txt='HUM_neg_top500_location.txt'

#set output file name
output.rds='HUM_pos_top500_extend1MB_location.RDS'

#import input
input.pos=read.table(input.pos.txt,header=F)
input.neg=read.table(input.neg.txt,header=F)

#make empty list of 2 vectors
input=vector(len=2,mode='list')
#give names to each vector
names(input)=c('pos','neg')
#make first vector have info on hypermethylated (positive) CpGs
input[[1]]=input.pos
#make second vector have info on hypomethylated (negative) CpGs
input[[2]]=input.neg
#view top of input list
head(input)

#import background dataset
background=read.table(background.txt,header=F)
#view
head(background)

##run the actual GREAT analysis

#set parameter for GREAT job
extend=c(50,1000)
#demo with the default extension of gene regulatory region(1 Mb)
job = submitGreatJob(input[[1]], bg = background,
                     species               = "hg19",
                     includeCuratedRegDoms = TRUE,
                     rule                  = c("basalPlusExt"),
                     adv_upstream          = 5.0,
                     adv_downstream        = 1.0,
                     adv_span              = extend[2],
                     request_interval = 300,
                     version="3.0.0",
                     max_tries = 10)

#run for hypermethylated (positive) CpGs
job.pos = readRDS(system.file("extdata", "job.rds", package = "rGREAT"))
#choose all availble ontologies
ontology.all=availableOntologies(job.pos)

#run loop
output.all={}
for(k in 1:length(ontology.all)){
  print(ontology.all[k])
  out0.list = tryCatch(getEnrichmentTables(job,ontology=ontology.all[k]),error=function(e){NULL})
  if(!is.null(out0.list)){
    db0.list=as.list(names(out0.list))
    output<-Map(cbind,Database=db0.list,out0.list)
    output<-do.call('rbind',output)
    output.all=rbind(output.all,output)
  }
}

#see all the databases that came out
table(output.all$Database)
#make data frame from GREAT output
output.all=data.frame(class='pos',extend='1Mb',output.all)
#view
head(output.all)
#save
saveRDS(output.all,file=output.rds)
#save as csv on computer
write.csv(output.all,"output.pos_BPI.csv")

#look at all the possible ontologies that can be checked for enrichment
#next steps are all for hypermethylated CpGs until we run the GREAT analysis with negative (hypomethylated) dataset
availableOntologies(job.pos)

#look specifically at the mouse phenotypes
mouse_pos = getEnrichmentTables(job, download_by = 'tsv', ontology = c("Mouse Phenotype"))
#the mouse_pos will be a list, so extract dataframe that we need
mouse_pos_df <- mouse_pos[["Mouse Phenotype"]]
#add column to specify "positive"
mouse_pos_df$direction <- 'pos'

#do same as above for any ontology dataset of interest for pos (hypermethylated) CpGs
Go_cell_pos = getEnrichmentTables(job, download_by = 'tsv', ontology = c("GO Cellular Component"))
Go_cell_pos_df <- Go_cell_pos[["GO Cellular Component"]]
Go_cell_pos_df$direction <- 'pos'

Go_mol_pos = getEnrichmentTables(job, download_by = 'tsv', ontology = c("GO Molecular Function"))
Go_mol_pos_df <- Go_mol_pos[["GO Molecular Function"]]
Go_mol_pos_df$direction <- 'pos'

InterPro_pos = getEnrichmentTables(job, download_by = 'tsv', ontology = c("InterPro"))
InterPro_pos_df <- InterPro_pos[["InterPro"]]
InterPro_pos_df$direction <- 'pos'

MGI_Exp_Det_pos = getEnrichmentTables(job, download_by = 'tsv', ontology = c("MGI Expression: Detected"))
MGI_Exp_Det_pos_df <- MGI_Exp_Det_pos[["MGI Expression: Detected"]]
MGI_Exp_Det_pos_df$direction <- 'pos'

GO_bio_pos = getEnrichmentTables(job, download_by = 'tsv', ontology = c("GO Biological Process"))
GO_bio_pos_df <- GO_bio_pos[["GO Biological Process"]]
GO_bio_pos_df$direction <- 'pos'

panther_pos = getEnrichmentTables(job, download_by = 'tsv', ontology = c("PANTHER Pathway"))
panther_pos_df <- panther_pos[["PANTHER Pathway"]]
panther_pos_df$direction <- 'pos'

human_pos = getEnrichmentTables(job, download_by = 'tsv', ontology = c("Human Phenotype"))
human_pos_df <- human_pos[["Human Phenotype"]]
human_pos_df$direction <- 'pos'

MSig_Pert_pos = getEnrichmentTables(job, download_by = 'tsv', ontology = c("MSigDB Perturbation"))
MSig_Pert_pos_df <-MSig_Pert_pos[["MSigDB Perturbation"]]
MSig_Pert_pos_df$direction <- 'pos'

Disease_pos = getEnrichmentTables(job, download_by = 'tsv', ontology = c("Disease Ontology"))
Disease_pos_df <-Disease_pos[["Disease Ontology"]]
Disease_pos_df$direction <- 'pos'

#look at multiple at once and make dataset including them all
big_df_pos <- getEnrichmentTables(job, download_by = 'tsv', ontology = c("GO Molecular Function",            "GO Biological Process",           
                                                                         "GO Cellular Component",            "Mouse Phenotype" ,                
                                                                         "Human Phenotype" ,                 "Disease Ontology", 
                                                                         "MGI Expression: Detected", "MSigDB Perturbation",             
                                                                         "MSigDB Predicted Promoter Motifs", "MSigDB miRNA Motifs"))
big_df_pos <- do.call(rbind.data.frame, big_df_pos) 
big_df_pos$direction <- 'pos'


#run the actual GREAT analysis for hypomethylated (negative) CpGs now
#re-run same thing as before basically
output.rds='HUM_neg_top500_extend1MB_location.RDS'
job = submitGreatJob(input[[2]], bg = background,
                     species               = "hg19",
                     includeCuratedRegDoms = TRUE,
                     rule                  = c("basalPlusExt"),
                     adv_upstream          = 5.0,
                     adv_downstream        = 1.0,
                     adv_span              = extend[2],
                     request_interval = 30,
                     version="3.0.0",
                     max_tries = 10)
#run for hypomethylated (negative) CpGs
job.neg = readRDS(system.file("extdata", "job.rds", package = "rGREAT"))
#choose all available ontologies
ontology.all=availableOntologies(job.neg)

#run loop
output.all={}
for(k in 1:length(ontology.all)){
  print(ontology.all[k])
  out0.list = tryCatch(getEnrichmentTables(job,ontology=ontology.all[k]),error=function(e){NULL})
  if(!is.null(out0.list)){
    db0.list=as.list(names(out0.list))
    output<-Map(cbind,Database=db0.list,out0.list)
    output<-do.call('rbind',output)
    output.all=rbind(output.all,output)
  }
}

#see all the databases that came out
table(output.all$Database)
#make data frame from GREAT output
output.all=data.frame(class='neg',extend='1Mb',output.all)
#view
head(output.all)
#save
saveRDS(output.all,file=output.rds)
#save as csv on computer
write.csv(output.all,"output.neg_BPI.csv")

#look at all the possible ontologies that can be checked for enrichment
#next steps are all for hypomethylated CpGs
availableOntologies(job.neg)

#look specifically at the mouse phenotypes
mouse_neg = getEnrichmentTables(job, download_by = 'tsv', ontology = c("Mouse Phenotype"))
#the mouse_neg will be a list, so extract dataframe that we need
mouse_neg_df <- mouse_neg[["Mouse Phenotype"]]
#add column to specify "negative"
mouse_neg_df$direction <- 'neg'

#do same as above for any ontology dataset of interest for pos (hypermethylated) CpGs
Go_cell_neg = getEnrichmentTables(job, download_by = 'tsv', ontology = c("GO Cellular Component"))
Go_cell_neg_df <- Go_cell_neg[["GO Cellular Component"]]
Go_cell_neg_df$direction <- 'neg'

Go_mol_neg = getEnrichmentTables(job, download_by = 'tsv', ontology = c("GO Molecular Function"))
Go_mol_neg_df <- Go_mol_neg[["GO Molecular Function"]]
Go_mol_neg_df$direction <- 'neg'

InterPro_neg = getEnrichmentTables(job, download_by = 'tsv', ontology = c("InterPro"))
InterPro_neg_df <- InterPro_neg[["InterPro"]]
InterPro_neg_df$direction <- 'neg'

MGI_Exp_Det_neg = getEnrichmentTables(job, download_by = 'tsv', ontology = c("MGI Expression: Detected"))
MGI_Exp_Det_neg_df <- MGI_Exp_Det_neg[["MGI Expression: Detected"]]
MGI_Exp_Det_neg_df$direction <- 'neg'

GO_bio_neg = getEnrichmentTables(job, download_by = 'tsv', ontology = c("GO Biological Process"))
GO_bio_neg_df <- GO_bio_neg[["GO Biological Process"]]
GO_bio_neg_df$direction <- 'neg'

panther_neg = getEnrichmentTables(job, download_by = 'tsv', ontology = c("PANTHER Pathway"))
panther_neg_df <- panther_neg[["PANTHER Pathway"]]
panther_neg_df$direction <- 'neg'

human_neg = getEnrichmentTables(job, download_by = 'tsv', ontology = c("Human Phenotype"))
human_neg_df <- human_neg[["Human Phenotype"]]
human_neg_df$direction <- 'neg'

MSig_Pert_neg = getEnrichmentTables(job, download_by = 'tsv', ontology = c("MSigDB Perturbation"))
MSig_Pert_neg_df <- MSig_Pert_neg[["MSigDB Perturbation"]]
MSig_Pert_neg_df$direction <- 'neg'

Disease_neg = getEnrichmentTables(job, download_by = 'tsv', ontology = c("Disease Ontology"))
Disease_neg_df <-Disease_neg[["Disease Ontology"]]
Disease_neg_df$direction <- 'neg'

#look at multiple at once and make dataset including them all
big_df_neg <- getEnrichmentTables(job, download_by = 'tsv', ontology = c("GO Molecular Function",            "GO Biological Process",           
                                                                         "GO Cellular Component",            "Mouse Phenotype" ,                
                                                                         "Human Phenotype" ,                 "Disease Ontology", 
                                                                         "MGI Expression: Detected", "MSigDB Perturbation",             
                                                                         "MSigDB Predicted Promoter Motifs", "MSigDB miRNA Motifs"))
big_df_neg <- do.call(rbind.data.frame, big_df_neg)
big_df_neg$direction <- 'neg'

#bind positive and negative mouse dataframes
#select top significant hits to plot, can play around with # of pos or neg you want
mouse_pos_top <- top_n(mouse_pos_df, -7, HyperBonfP)
mouse_neg_top <- top_n(mouse_neg_df, -0, HyperBonfP)
mouse_df <- rbind(mouse_neg_top, mouse_pos_top)
#transform p-value to log10
mouse_df$p.adj <- -log10(mouse_df$HyperBonfP)

#load ggplot
library(ggplot2)
#plot significant mouse phenotype pathways
#x axis is direction of association (pos or neg), y axis is pathway, make dot color based on log10 of p value
mouse_plot <- ggplot(data = mouse_df, aes(x = direction, y = Desc, colour = mouse_df$p.adj )) + 
  geom_point()+
  #choose classic theme
  theme_classic() + 
  #no x and y labels
  ylab("") + 
  xlab("") + 
  #scale dot color
  #less significant p-values will be red and more significant will be blue
  scale_color_gradient(low = "red", high = "blue") +
  #update plot title and legend title
  ggtitle("Mouse Phenotype")  + labs(color='-log10(P)') 
#view
mouse_plot

#redo same thing for any ontology dataset of interest
Go_mol_pos_top <- top_n(Go_mol_pos_df, -3, HyperBonfP)
Go_mol_neg_top <- top_n(Go_mol_neg_df, -0, HyperBonfP)
Go_mol_df <- rbind(Go_mol_neg_top, Go_mol_pos_top)
Go_mol_df$p.adj <- -log10(Go_mol_df$HyperBonfP)
Go_mol_plot <- ggplot(data = Go_mol_df, aes(x = direction, y = Desc, colour = Go_mol_df$p.adj )) + 
  geom_point()+
  theme_classic() + 
  ylab("") + 
  xlab("") + 
  scale_color_gradient(low = "red", high = "blue") +
  ggtitle("Go_mol")  + labs(color='-log10(P)') 
Go_mol_plot

MGI_Exp_Det_pos_top <- top_n(MGI_Exp_Det_pos_df, -5, HyperBonfP)
MGI_Exp_Det_neg_top <- top_n(MGI_Exp_Det_neg_df, -4, HyperBonfP)
MGI_Exp_Det_df <- rbind(MGI_Exp_Det_neg_top, MGI_Exp_Det_pos_top)
MGI_Exp_Det_df$p.adj <- -log10(MGI_Exp_Det_df$HyperBonfP)
MGI_Exp_Det_plot <- ggplot(data = MGI_Exp_Det_df, aes(x = direction, y = Desc, colour = MGI_Exp_Det_df$p.adj )) + 
  geom_point()+
  theme_classic() + 
  ylab("") + 
  xlab("") + 
  scale_color_gradient(low = "red", high = "blue") +
  ggtitle("Go_mol")  + labs(color='-log10(P)') 
MGI_Exp_Det_plot

InterPro_pos_top <- top_n(InterPro_pos_df, -2, HyperBonfP)
InterPro_neg_top <- top_n(InterPro_neg_df, -0, HyperBonfP)
InterPro_df <- rbind(InterPro_neg_top, InterPro_pos_top)
InterPro_df$p.adj <- -log10(InterPro_df$HyperBonfP)
InterPro_plot <- ggplot(data = InterPro_df, aes(x = direction, y = Desc, colour = InterPro_df$p.adj )) + 
  geom_point()+
  theme_classic() + 
  ylab("") + 
  xlab("") + 
  scale_color_gradient(low = "red", high = "blue") +
  ggtitle("InterPro")  + labs(color='-log10(P)') 
InterPro_plot

GO_bio_pos_top <- top_n(GO_bio_pos_df, -5, HyperBonfP)
GO_bio_neg_top <- top_n(GO_bio_neg_df, -0, HyperBonfP)
GO_bio_df <- rbind(GO_bio_neg_top, GO_bio_pos_top)
GO_bio_df$p.adj <- -log10(GO_bio_df$HyperBonfP)
GO_bio_plot <- ggplot(data = GO_bio_df, aes(x = direction, y = Desc, colour = GO_bio_df$p.adj )) + 
  geom_point()+
  theme_classic() + 
  ylab("") + 
  xlab("") + 
  scale_color_gradient(low = "red", high = "blue") +
  ggtitle("GO_bio")  + labs(color='-log10(P)') 
GO_bio_plot

Disease_pos_top <- top_n(Disease_pos_df, -5, HyperBonfP)
Disease_neg_top <- top_n(Disease_neg_df, -0, HyperBonfP)
Disease_df <- rbind(Disease_neg_top, Disease_pos_top)
Disease_df$p.adj <- -log10(Disease_df$HyperBonfP)
Disease_plot <- ggplot(data = Disease_df, aes(x = direction, y = Desc, colour = Disease_df$p.adj )) + 
 geom_point()+
 theme_classic() + 
  ylab("") + 
  xlab("") + 
  scale_color_gradient(low = "red", high = "blue") +
  ggtitle("Disease Ontology")  + labs(color='-log10(P)') 
Disease_plot

###final plotting
#decided to keep "GO Molecular Function", "InterPro", "MGI Expression", and "Mouse Phenotype" ontologies
#bind all four datasets together
whole_df <- rbind(Go_mol_df, InterPro_df, MGI_Exp_Det_df, mouse_df)
#make a new column with HyperBonferroni p values changed to log10
whole_df$p.adj <- -log10(whole_df$HyperBonfP)
#make new column with ratio of genes hit from pathway/total genes in pathway
whole_df$ratio <- whole_df$NumFgGenesHit/whole_df$TotalGenes

#save file as csv on computer
write.csv(whole_df,file="GREAT_df_BPI.csv")


#plot
#x axis is direction of association (pos or neg), y axis is pathway
#make dot color based on log10 of p value and dot size based on ratio of genes hit
whole_df_plot <- ggplot(data = whole_df, aes(x = direction, y = Desc, colour = whole_df$p.adj, size=whole_df$ratio )) + 
  geom_point()+
  #choose classic theme
  theme_classic() +
  #make plot title bold and centered
  theme(plot.title = element_text(hjust = 0.5, face="bold")) +
  #with str_wrap adjust the pathways' width of writing
  aes(stringr::str_wrap(Desc, 50), direction) + ylab(NULL) +
  #flip axes
  coord_flip()+
  #remove x and y labels
  ylab("") + 
  xlab("") + 
  #scale dot color
  #less significant p-values will be red and more significant will be blue
  scale_color_gradient(low = "red", high = "blue") +
  #update title and legend names
  ggtitle("Top Pathways")  + labs(color='-log10(P)') + labs(size='Gene Ratio') +
  #separate plot in 4 based on initial ontology datasets using facet_wrap, place in only one column
  facet_wrap(~Ontology, scales = "free", ncol = 1)
#view
whole_df_plot


####FINAL PLOT (FIGURE 4 IN PAPER)
#load cowplot to be able to merge all 4 plots in one figure
library(cowplot)

#to create this final figure you will need to have run the 4.2_EWAS_BPI_effect script and have the "DNAm_Location_plot" in your environment
#you will also need to have run the 5.2_Manhattan_plot_BPI_effect script and have the "manhplot" plot in your environment

#make the left part of Figure 4 with the manhattan plot and boxplots
#label manhattan plot "A" and boxplot "B", but them one on top of the other with cols = 1
#align both plots horizontally and vertically with "hv" flag, and based on left and right margind with "lr" flag
left <- plot_grid(manhplot, DNAm_Location_plot , labels = c('A', 'B'), cols = 1, align = "hv", axis = "lr")

#include pathway plot to figure and label it "C"
#basically merge the left part of the figure we just created with pathway in two different columns
#adjust space each take on the figure with rel_widths (make pathway bit smaller)
Figure_4 <- plot_grid(left, whole_df_plot , labels = c('',  'C'), cols = 2, rel_widths = c(2, 1.25))
#view final Figure 4 and manually export as pdf
Figure_4

