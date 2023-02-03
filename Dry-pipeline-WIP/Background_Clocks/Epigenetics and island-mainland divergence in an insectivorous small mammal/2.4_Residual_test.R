##R Script to test Epigenetic Clock residuals (leave-one-out)
#Main author: Marie-Laurence Cossette

#set working directory where all files are in 
setwd("~/Desktop/SCHOOL/THESIS/Chapter 2 - Methyl/PAPER/Gitlab uploads")

#load data
#only used basic clock (instead of sqrt and log clocks)
residual_df=read.csv('Clock/SCsub_ClockBasedOnAll_EpigeneticAge_PredictedValues.csv', as.is=T)
#shrew info from Steve Horvath
data = read.csv('SampleSheetAgeN90final.csv')

#remove fetal samples
data=subset(data, Age>0.1)
#add a location column to separate island from mainland (Island = BPI and Long Islane, Mainland = all others)
data<- mutate(data, Location_1 = ifelse(data$Location == "Bon Portage Island Nova Scotia" | data$Location == "Long Island Nova Scotia", "yes", "no"))  
#merge residuals data set and shrew info dat set
df <- merge(x = residual_df, y = data, by = "ExternalSampleID")

#T-test to see if different groupings epigenetically aged at different rates 
#sex
t.test(AgeAccelbasedOnAll.x ~ Sex, data=df)

#location
t.test(AgeAccelbasedOnAll.x ~ Location_1, data=df)

