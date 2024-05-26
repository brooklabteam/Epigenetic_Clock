library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)

rm(list=ls())

homewd = "/Users/sophiahorigan/Desktop"
setwd(paste0(homewd, "/ClockFoundationSamples"))

#first load the recaps
recap.df <- read.csv(file = "pre2020recaps.csv", header = T, stringsAsFactors = F)
head(recap.df)
#subset the recaps to those after 2017 and before 2021 (these are the WPs on hand)
recap.df = subset(recap.df, processing_date < "2021-01-01" & processing_date > "2017-01-01")
nrow(recap.df) #235 samples. 37K
#what about not including rousettus?
nrow(recap.df[recap.df$bat_species!="Rousettus madagascariensis",]) 

eid.recap = subset(recap.df, bat_species != "Rousettus madagascariensis")
#121 samples. 19K. better. Do all these Eidolon for sure!

# subtotal: 121 recap eidolon


#load the tooth or juvenile aged ones
age.df <- read.csv(file = "CF-toothjuv.csv", header = T, stringsAsFactors = F)

#add in link to file with date and remove those from pre 2017 (no WPs)
catch.dat <- read.csv(file = "catching_data_clean.csv", header = T, stringsAsFactors = F)
head(catch.dat)

unique(catch.dat$processing_date)
catch.merge <- dplyr::select(catch.dat, sampleid, processing_date)

#merge with age
age.df <- merge(age.df, catch.merge, by="sampleid", all.x = T)
head(age.df)
sort(unique(age.df$processing_date))

#you only have wing punches or DNA on hand for bats caught after 2017 and before 2021
age.df = subset(age.df, processing_date < "2021-01-01" & processing_date > "2017-01-01")

#how many total?
age.sum <- ddply(age.df, .(species), summarise, N=length(sampleid))
age.sum

#207 Eidolon, 96 Pteropus, 133 Rousettus - this is about half the estimate from the meeting
#that is 303 Eidolon + Pteropus (48K)

#15K just for the Pteropus. How many are juveniles?
pter.age = subset(age.df, species == "Pteropus_rufus")
length(pter.age$sampleid[pter.age$aging_method=="juvenile"]) #40 juvenile + 56 adult tooth
pter.age.adult = subset(pter.age, aging_method=="tooth")
pter.age.juv = subset(pter.age, aging_method=="juvenile")

#9K for the Pteropus tooth aged bats. Do ALL OF THESE!

#then, for the juveniles, subset... we've got 40... maybe do 20?

#which we do should be based on infection history.


# new subtotal:
# 121 recap eidolon
# 56 tooth-aged + 20 juvenile-aged pteropus


#are there any overlap between aged bats and recaps?
recap.id <- dplyr::select(recap.df, sampleid)
recap.id$recap = 1

age.recap.df <- merge(age.df, recap.id, by="sampleid", all.x = T)
head(age.recap.df)

#here is the high priority subset of bats that are recaps and known age!
high.priority = subset(age.recap.df, recap==1)
length(high.priority$sampleid[high.priority$species=="Eidolon_dupreanum" & high.priority$aging_method=="tooth"]) #24

sum.highp <- ddply(high.priority, .(species), summarise, N=length(sampleid))
sum.highp #28 Eidolon and 9 Rousettus - definitely do these! 
# about 5K saved between overlap of known age Eidolon and recap Eidolon
rou.age.recap <- subset(high.priority, species=="Rousettus_madagascariensis")

#new subtotal:
# 121 recap eidolon (28 of which are known age)
# 56 tooth-aged + 20 juvenile-aged pteropus
# 9 aged + recap rousettus


#last Q is how many more of the remaining 179 known age eidolon to do?

#how many are juvenile vs. adult?
length(age.df$sampleid[age.df$species == "Eidolon_dupreanum" & age.df$aging_method=="juvenile"]) #78
length(age.df$sampleid[age.df$species == "Eidolon_dupreanum" & age.df$aging_method=="tooth"]) 
#129. but 28 of these were already accounted for in the recaps
# so an additional 17K for the eidolon adults

eid.aged.adult = subset(age.df, species=="Eidolon_dupreanum" & aging_method =="tooth")
eid.aged.juv = subset(age.df, species=="Eidolon_dupreanum" & aging_method =="juvenile")

#new subtotal:
# 121 recap eidolon (28 of which are known age)
# 56 tooth-aged + 20 juvenile-aged pteropus
# 9 aged/recap rousettus (are these juveniles?)
# 105 tooth-aged eidolon (not including the 24 already included in the recaps)
#### total = 291 = 46K


#then, the question is which of the eidolon juveniles to send?

#might make sense to look for even representation by sex. 
#also look at season (a distribution of ages across early life)
#can also look at infection data (below)


#next, look at the infection history for the aged bats
head(age.df)

#and link the luminex data to the aged subset
#note that this only has seropositives for Rousettus and Eidolon
lum.dat <- read.csv(file = "lumdat_clust80.csv", header = T, stringsAsFactors = F)
head(lum.dat)

#first, drop MFI column and make into a wide dataset
lum.dat <- dplyr::select(lum.dat, -(MFI))
head(lum.dat)
lum.wide <- dcast(data=melt(lum.dat), formula = sampleid~antigen)
head(lum.wide) 
#lum.wide$species <- NA
#lum.wide$species[is.na(lum.wide$BDBV) & is.na(lum.wide$RAVV)] <- "Eidolon_dupreanum"
#lum.wide$species[is.na(lum.wide$species)] <- "Rousettus_madagascariensis"

head(lum.wide) 

#then, join to your aged dataset above
age.df <- merge(age.df, lum.wide, by="sampleid", all.x = T)
head(age.df)

#and look at the juvenile eidolon you are trying to decide on
age.sub.eid = subset(age.df, species=="Eidolon_dupreanum" & age.df$aging_method=="juvenile") #78
length(age.sub.eid$sampleid[age.sub.eid$extracted=="Y" & age.sub.eid$AngV==1]) #25 AngV +
length(age.sub.eid$sampleid[age.sub.eid$extracted=="Y" & age.sub.eid$AngV==0]) #46 AngV -

length(age.sub.eid$sampleid[age.sub.eid$extracted=="Y" & age.sub.eid$GhV==1]) #34 GhV +
length(age.sub.eid$sampleid[age.sub.eid$extracted=="Y" & age.sub.eid$GhV==0]) #37 AngV -

#maybe take 20 of these juveniles in addition? brings us to about an even 50K


#FINAL COUNTS
#121 recap Eidolon (SET) - recap.df
eid.recap$category = 'recap'
#56 tooth-aged pteropus (SET) 
pter.age.adult$category = 'toothaged'
#129 (101) tooth-aged eidolon (including the 28 already included in the recaps) (SET) - age.df
eid.aged.adult$category = "toothaged" #contains duplicates when merged with recaps

#merge eidolon recap and tooth aged to remove duplicates

#ALL ADULTS  = 282
merge1 <- merge(eid.recap,rou.age.recap, by=c("sampleid", "processing_date"), all = TRUE)
merge2 <- merge(merge1, pter.age.adult, by=c("sampleid", "processing_date"), all = TRUE)
CF.final.adult <- merge(merge2, eid.aged.adult, by=c("sampleid","processing_date"), all = TRUE) #this merge removes the duplicates
#get just sampleid row and category
CF.final.adult.id <- CF.final.adult[c("sampleid","processing_date")]
#and merge with catching data to all metadata
CF.final.adult.metadata <- merge(CF.final.adult.id, catch.dat, by=c("sampleid", "processing_date"), all.x=TRUE)
#move to excel for fine tuning and removing duplicates
write.csv(CF.final.adult.metadata, file=paste0(homewd,"/ClockFoundationSamples/CFadults.csv"), row.names = T)

#JUVENILES TO CHOOSE FROM (BY HAND) = 118
#40 juvenile-aged pteropus (CHOOSE 20) - age.df
pter.age.juv$category = "juvenile"
#78 juvenile eidolon (CHOOSE 21) - age.sub.eid
eid.aged.juv$category = "juvenile"
#merge
CF.final.juv <- merge(pter.age.juv, eid.aged.juv, by=c("sampleid", "processing_date"), all=TRUE)
#get just sample id and category
CF.final.juv.id <- CF.final.juv[c("sampleid", "processing_date")]
#merge with catching data to get all metadata
CF.final.juv.metadata <- merge(CF.final.juv.id, catch.dat, by=c("sampleid", "processing_date"), all.x=TRUE)
#export to excel
write.csv(CF.final.juv.metadata, file=paste0(homewd, "/ClockFoundationSamples/CFjuv.csv"), row.names=T)
#FOR JUVENILES
#distribution of ages, sex
#gwen's list of virus positive bats (filo for example)
#locations: 5 from ankarana, 15 from angavokely (eidolon), cross reference with adult locations (reflect pattern in adults) (121 and 105)
#if all else equal, take extracted

## WHAT WE KNOW FROM CLOCK FOUNDATION, WHAT OTHER QUESTIONS DO WE HAVE? EMAIL BOBBY

#Next steps:
#1. create final list of samples (TODAY)
#2. locate and compile TO BE EXTRACTED (FRIDAY/MONDAY) (based on google stocking sheet)
#3. locate and compile EXTRACTED (FRIDAY/MONDAY) (based on )
#4. Extract non-extracted (put into final plates) (TUESDAY/WEDNESDAY)
#5. Move extracted into plates (TUESDAY/WEDESDAY)
#6. Seal plates and ship! (THURSDAY/FRIDAY)


