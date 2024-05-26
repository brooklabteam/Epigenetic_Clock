library(dplyr)
library(ggplot2)

# add ID 
mammalianarray_default_batpredictions$BrookLabID <- sub1_IDlist$sampleid
colnames(mammalianarray_default_batpredictions)[29] <- "sampleid"

# add age 
tmp <- inner_join(mammalianarray_default_batpredictions, aged.since.2018, by = "sampleid")

# regress
# bat skin clocks - all aged bats
par(mfrow = c(3,1))
ggplot(data = tmp, aes(x=age, y = UniversalClock2Pan)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) + ylim(0,17) + xlim(0,17) 

ggplot(data = tmp, aes(x=age, y = DNAmAgeBat.Eidolon.helvum_Skin)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) + ylim(0,13) + xlim(0,13) +
  annotate("text",x=10,y=10,label=(paste0("slope==",coef(lm(tmp$DNAmAgeBat.Eidolon.helvum_Skin~tmp$age))[2])),parse=TRUE)


fit <- lm(DNAmAgeBat.BatSkin~age, tmp)

ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") + xlab("Estimated Age - Dentition") + ylab("Estimated Age - Epigenetic Clock") +
  ylim(0,20) + xlim(0,20) +
  labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 5),
                     "Intercept =",signif(fit$coef[[1]],5 ),
                     " Slope =",signif(fit$coef[[2]], 5),
                     " P =",signif(summary(fit)$coef[2,4], 5)))


# pteropus only
ggplot(subset(tmp, species %in% "Pteropus_rufus"), aes(x=age, y = DNAmAgeBat.BatPteropusSkin)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) + ylim(0,17) + xlim(0,17)  

ggplot(subset(tmp, species %in% "Pteropus_rufus"), aes(x=age, y = DNAmAgeBat.BatSkin)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) + ylim(0,20) + xlim(0,20) 

# eidolon only
ggplot(subset(tmp, species %in% "Eidolon_dupreanum"), aes(x=age, y = DNAmAgeBat.Eidolon.helvum_Skin)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) + ylim(0,13) + xlim(0,13) 

ggplot(subset(tmp, species %in% "Eidolon_dupreanum"), aes(x=age, y = DNAmAgeBat.BatSkin)) +
  geom_point() +
  geom_smooth(method = "lm", se = F) + ylim(0,20) + xlim(0,20) 




# add ID 
mammalianarray_default_predictions$BrookLabID <- sub1_IDlist$sampleid
colnames(mammalianarray_default_predictions)[112] <- "sampleid"

# add age 
tmp <- inner_join(mammalianarray_default_predictions, aged.since.2018, by = "sampleid")

# universal clocks
ggplot(data = tmp, aes(x=age, y = DNAmAgePanFinal)) +
  geom_point() +
  geom_smooth(method = "lm")


