
# The dataframe (df) should be organized with rows = samples, columns = CpGs,... 
# ...and a few metadata columns that need to be removed before clock training.

# A rough outline of the function is:
# (1) Randomly select samples from discrete age groups for the model training set (≈80%).
# (2) Filter training set to get testing set matrix.
# (3) Apply elastic net regression to train age prediction model using training data
# (4) Predict ages with the model in the testing set
# (5) Gather stats, write clock coefficeints, plot results.

trainTest_clock <- function(df,seed,alpha,tissues) {
  
  print(seed) #random number 
  set.seed(seed)
  print(alpha) #elastic net regression parameter (between lasso and ridge) #Wilkinson used 0.5, Griffin used 0.05 or 0.1
  #you can refine and train alpha -- see notes in paper
  
  tissue_names <- paste(c(tissues), collapse='_' ) #we don't have tissues so this is irrelevant?
  #they built clocks based on different tissues, so our 'tissue' of interest is just dentition
  
  ##DISCRETIZE AGE GROUPS
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
  
  tr_fr <- 0.8 # the training samples randomly selected from each age_group
  ts_fr <- (1-tr_fr) #testing fraction
  
  training_samples <- x %>% filter(sample=="none") #keeps rows with sample == none #bc those are un-dated? #or already subset?
  
  ##SUBSET TRAINING SET
  for (j in  unique(x$Tissue)) { #don't need this loop
    #selecting samples from initial age group to start building training matrix
    filtered <- x %>%  # filtered is age group
      filter(age_group==unique(x$age_group)[1] & Tissue==j) #filter from age group 1, tissue type 1
    tr_samps <- filtered %>% #tr_samps is training 80% of age group
      sample_n(as.integer((tr_fr*length(age_group)))) #sample training proportion from group (based on tr_fr)
    
    #now selecting samples from each age group
    for (i in unique(x$age_group)[2:length(unique(x$age_group))]){ # now subset the rest of the age groups one by one
      f <- x %>%  # f is age group
        filter(age_group==i & Tissue==j)
      s <- f %>% # s is 80% of f
        sample_n(as.integer((tr_fr*length(age_group))+0.5)) #and randomly sample the rest of the training data
      tr_samps <- bind_rows(tr_samps,s) #bind s (training for looped age group) to tr_samps (training for age group 1)
    }
    training_samples <- bind_rows(training_samples,tr_samps) # bind everything to make final df of training samples
  }
  
  #FORMAT DATA SETS FOR CLOCK ANALYSIS
  #training data
  #y is predictor, which age Age in months
  y.train <- training_samples %>% select(Age_months) #takes only the age 
  y.train <- as.matrix(y.train) 
  #x is independent, which are methylation levels at CpG sites
  x.train <- training_samples %>% select(-sample,-Age_months,-age_group,-Tissue) #drop sample, Age, age group, Tissue #which leaves CPGS?
  x.train <- model.matrix(  ~ ., x.train) #makes a 'model matrix', which is a matrix of values of explanatory variables #basically all the explanatory variables and observations
  
  #testing data
  x.ts <- x %>% filter(sample %in% training_samples$sample == F) #keep the data that wasn't used in the training set #still don't understand if he put 'none' in there or what
  y.test <- x.ts %>% select(Age_months)
  y.test <- as.matrix(y.test) 
  x.test <- x.ts %>% select(-sample,-Age_months,-age_group,-Tissue) 
  x.test <- model.matrix( ~ ., x.test)
  
  ## MODELING
  a <- alpha # alpha parameter for elastic net regression
  
  set.seed(42)
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
                  AGE=training_samples$Age_months, #what were the actual ages
                  sample=training_samples$sample, #ID? or CpG site?
                  set="train",
                  Tissue=training_samples$Tissue) #don't need
  test <- tibble(prediction = predicted.test, #predicted ages
                 AGE= y.test, #actual ages
                 sample=x.ts$sample,
                 set="test",
                 Tissue=x.ts$Tissue)
  # combine trained and tested predictions
  all_prds <- bind_rows(train,test) %>% mutate(set=factor(set,levels = c("train","test")))
  #export data
  write_tsv(all_prds, path=paste(wd,'clock_optimization/predictions_',
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
  clock_1 <- data.frame(name=tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1],
                        coefficient = tmp_coeffs@x)
  
  write_tsv(clock_1, path=paste(wd,'clock_optimization/',tissue_names,
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
    ggtitle(paste("MdAE=", p_med, " N_CpGs=",length(clock_1$coefficient))) +
    xlim(-10,40) +
    ylim(-10,40)
  print(g)
  ggsave(plot=g,
         height=5,
         width=10,
         filename=paste(wd,'clock_optimization/',tissue_names,'_seed=',seed,'_alpha=',a,
                        '_train=',100*tr_fr,"_test=",
                        100*ts_fr,".pdf",sep=""), useDingbats=FALSE)
  
  return(predictions)
  print(predictions)
}
