##################################################################################
# TRAINING CLOCK
# This is code to train an epigenetic clock
# it is currently built to use Myotis samples (see buildDataframes_WIP)
# email shorigan@uchicago.edu, tlaverty@nmsu.edu, or cbrook@uchicago.edu with questions


# Note from Griffin et al
# The dataframe (df) should be organized with rows = samples, columns = CpGs,... 
# ...and a few metadata columns that need to be removed before clock training.

# A rough outline of the function is:
# (1) Randomly select samples from discrete age groups for the model training set (≈80%).
# (2) Filter training set to get testing set matrix.
# (3) Apply elastic net regression to train age prediction model using training data
# (4) Predict ages with the model in the testing set
# (5) Gather stats, write clock coefficeints, plot results.

##################################################################################
# PART 1: SET UP ENVIRONMENT
# Clear work environment
rm(list=ls())

# Set working directory - add yours here
#homewd= "/Users/theresalaverty/Documents/R/R_repositories/Epigenetic_Clock"
homewd= "/Users/shorigan/Documents/GitHub/Epigenetic_Clock"
# Folder where the clock site csvs are stored
setwd(paste0(homewd,"/Dry-pipeline-WIP/Dataframes/"))

# load packages
#install.packages("glmnet")
#install.packages("ggpubr")
library(glmnet)
library(readr)
library(ggpubr)
library(dplyr)

# read in files - one for each type of site (for now)
# CPG clock sites
cpg_clock_sites <- read.csv(file= paste0(homewd,"/Dry-pipeline-WIP/Dataframes/cpg_clock.csv"), header = T)
# CHG clock sites
chg_clock_sites <- read.csv(file= paste0(homewd,"/Dry-pipeline-WIP/Dataframes/chg_clock.csv"), header = T)
# CHH clock sites
chh_clock_sites <- read.csv(file= paste0(homewd,"/Dry-pipeline-WIP/Dataframes/chh_clock.csv"), header = T)

# set clock parameters
# random seed
seed = 777 # does this need to be generated new each time, probably yes
# alpha for elastic net regression (0-1) # Wilkinson used 0.5, Griffin used 0.05 or 0.1 # can also try to fit this value
a = 0.5

df <- cpg_clock_sites
# NEED TO CHANGE
#just for now to get the clock working-- removing columns (sites) with 7 or more NAs
df <- df[,colSums(is.na(df)) < nrow(df)-2]

#trainTest_clock <- function(df,seed,a) {
  
  print(seed)
  set.seed(seed)
  print(a)

  # create subsets of the data to be used for differential clock building # actually, I think it's just to ensure random sampling for clock building
  # optional! # make one for CPG, CHH, CGH in future
  tissue_names <- paste(c(df$tissue)) #we don't have tissues so this is irrelevant?
  age_methods <- paste(c(df$age_method))

  ##DISCRETIZE AGE GROUPS
  # see code below, not using for now
  
  tr_fr <- 0.8 # the training samples randomly selected from each age_group
  ts_fr <- (1-tr_fr) #testing fraction
  
  training_samples <- df %>% filter(id=="none") # creates empty dataframe of correct size
  
  ##SUBSET TRAINING SET
  for (j in  unique(df$tissue)) { #don't need this loop
    #selecting samples from initial age group to start building training matrix
    filtered <- df %>%  # filtered is age group
      filter(known_age==unique(df$known_age)[1] & tissue==j) #filter from age group 1, tissue type 1
    tr_samps <- filtered %>% #tr_samps is training 80% of age group
      sample_n(as.integer((tr_fr*length(known_age)))) #sample training proportion from group (based on tr_fr)
    
    #now selecting samples from each age group
    for (i in unique(df$known_age)[2:length(unique(df$known_age))]){ # now subset the rest of the age groups one by one
      f <- df %>%  # f is age group
        filter(known_age==i & tissue==j)
      s <- f %>% # s is 80% of f
        sample_n(as.integer((tr_fr*length(known_age))+0.5)) #and randomly sample the rest of the training data
      tr_samps <- bind_rows(tr_samps,s) #bind s (training for looped age group) to tr_samps (training for age group 1)
    }
    training_samples <- bind_rows(training_samples,tr_samps) # bind everything to make final df of training samples
  }
  
  #FORMAT DATA SETS FOR CLOCK ANALYSIS
  ### ISSUE WHEN SPLITTING TRAIN AND TEST-- the split leaves train data with some cpgs that have all NA's
  # need to either filter secondarily (then go back and re-sample test.. seems like a hassle)
  # or remove all sites where only one sample hit the probe and hope the split/test works
  # this will be as an issue with larger sample sizes as well. Creating a more stringent cutoff is necessary I think.?
  
  #training data

  y.train <- training_samples %>% select(known_age)
  y.train <- as.matrix(y.train)
  
  x.train <- training_samples %>% select(-id,-known_age,-tissue,-age_method) #drop sample, Age, age group, Tissue #which leaves CPGS?
  x.train <- as.data.frame(x.train)
  x.train <- makeX(x.train, na.impute=TRUE) # this generates an average for the NA sites. maybe better to do elsewhere
  
  #testing data
  df.ts <- df %>% filter(id %in% training_samples$id == F) # need to make sure samples are individually id'd
  y.test <- df.ts %>% select(known_age)
  y.test <- as.matrix(y.test) 
  
  x.test <- df.ts %>% select(-id,-known_age,-tissue,-age_method)
  x.test <- as.data.frame(x.test)
  #x.test <- model.matrix( ~ ., x.test)
  x.test <- makeX(x.test, na.impute = TRUE)  # not sure exactly what sparse does
  

  
  ## MODELING
  a <- 0.5 # alpha parameter for elastic net regression
  set.seed(42) # do we need this part twice?
  
  
  #fit a generalized linear model via penalized maximum likelihood
  #using cross validation (cv)
  cv.glmnet.fit <- cv.glmnet(x.train,y.train, #observations, response
                             type.measure = "mae", #loss to use for cross-validation #mean absolute error- measures the deviation from the fitted mean to the response
                             alpha=a, 
                             family="gaussian") #is our family gaussian
  #predicts fitted values, logits, coefficients from fitted glmnet object
  #why is he doing this here if he does it below
  cv.glmnet.predicted <- predict(cv.glmnet.fit,
                                 s = cv.glmnet.fit$lambda.1se, #values of penalty parameter lambda (default is use entire sequence used to create the model, what he's doing here)
                                 newx = x.test ) #matrix of new values for x at which predictions are to be made
  
  ##IS THIS JUST HOW I USE AN ALREADY FIT CLOCK???
  ##just needs the output of the fit (coefficients?) and lambda-- is that somewhere?
  
  
  
  #predicting for train and test using the model
  predicted.train <- predict(cv.glmnet.fit,
                             s = cv.glmnet.fit$lambda.1se, 
                             newx = x.train)
  predicted.test <- predict(cv.glmnet.fit, 
                            s = cv.glmnet.fit$lambda.1se,  
                            newx = x.test)
  
  
  # compile: model prediction of age, actual age, metadata (sample, tissue)
  train <- tibble(prediction = predicted.train,
                  AGE=training_samples$known_age, 
                  sample=training_samples$id, 
                  set="train",
                  Tissue=training_samples$tissue) 
  test <- tibble(prediction = predicted.test, #predicted ages
                 AGE= y.test, 
                 sample=df.ts$id,
                 set="test",
                 Tissue=df.ts$tissue)
  # combine trained and tested predictions
  all_prds <- as.data.frame(bind_rows(train,test) %>% mutate(set=factor(set,levels = c("train","test"))))
  #export data
  write.csv(all_prds, file= paste0(homewd,'/Dry-Pipeline-WIP/clock_optimization/predictions_',
                                   tissue_names,
                                   '_seed=',seed,'_alpha=',a,
                                   '_train=',100*tr_fr,"_test=",
                                   100*ts_fr,".tsv",sep=""))
  
  #calculating stats from prediction
  p_mse <- mean((y.test - cv.glmnet.predicted)^2)
  p_mae <- mean(abs(y.test - cv.glmnet.predicted))
  p_med <- median(abs(y.test - cv.glmnet.predicted))
  c <- cor(test$prediction,test$AGE,method = "pearson")[,1]
  rh <- cor(test$prediction,test$AGE,method ="spearman")[,1]
  
  #gathering the stats into a tibble
  predictions <- tibble(seed=NA, alpha= NA, mse = NA,mae = NA, cor=NA, rho=NA, mdae=NA)
  p <- tibble(seed=seed, alpha=a, mse = p_mse, mae = p_mae, cor=c, rho=rh, mdae=p_med)
  predictions <- bind_rows(predictions,p) %>% filter(is.na(alpha) == F) 
  tmp_coeffs <- coef(cv.glmnet.fit, s="lambda.1se")
  clock1 <- data.frame(name=tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1],
                        coefficient = tmp_coeffs@x)
  
  write.csv(clock1, file= paste0(homewd,'/Dry-Pipeline-WIP/clock_optimization/',tissue_names,
                                 '_seed=',seed,'_alpha=',a,
                                 '_train=',100*tr_fr,"_test=",
                                 100*ts_fr,".tsv",sep=""))
  
  
  g <- all_prds %>%
    ggplot() +
    geom_smooth(aes(AGE,prediction,color=Tissue), alpha=0.5,method="lm",se=F) +
    geom_point(aes(AGE,prediction,color=Tissue), alpha=0.5) +
    theme_bw() +
    facet_wrap(.~set) +
    geom_abline(slope=1,linetype="dashed") +
    stat_cor(aes(AGE,prediction),vjust=0.1) +
    stat_cor(aes(AGE,prediction,color=Tissue)) +
    ggtitle(paste("MdAE=", p_med, " N_CpGs=",length(clock1$coefficient))) +
    xlim(0,8) +
    ylim(0,8)
  print(g)
  ggsave(plot=g,
         height=5,
         width=10,
         filename=paste(wd,'clock_optimization/',tissue_names,'_seed=',seed,'_alpha=',a,
                        '_train=',100*tr_fr,"_test=",
                        100*ts_fr,".pdf",sep=""), useDingbats=FALSE)
  
  return(predictions)
  print(predictions)
#}
  
  
  
  ## DISCRETIZE AGE GROUPS
  #filtering for the tissues of interest and assigning discrete age-groups.
  #we need to just assign discrete age-groups
  x <- df %>%
    filter(Tissue %in% tissues) %>% #shouldn't need to do this either, just do it by age year
    mutate(age_group = case_when(Age_months<2.5 ~ "3", #we do have indls less than one year
                                 Age_months>=3 & Age_months<6 ~"3-6",
                                 Age_months>=6 & Age_months<10 ~"6-10",
                                 Age_months>=10 & Age_months<14 ~"10-14",
                                 Age_months>=14 & Age_months<18 ~"14-18",
                                 Age_months>=18 & Age_months<24 ~"18-24",
                                 Age_months>=24 & Age_months<28 ~"24-28",
                                 Age_months>=28 & Age_months<32 ~"28-32",
                                 Age_months>=32 & Age_months<34 ~"32-34",
                                 Age_months>=34~">34")) 
