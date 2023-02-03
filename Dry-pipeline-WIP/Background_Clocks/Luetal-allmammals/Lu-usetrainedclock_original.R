#
#Rcode
# Lu et al
# all mammal clock
# SH: I am going to use this much easier to read script to figure out how to use 
# epigeneitc clocks that are already made, so we can apply to bat clocks
# uses glmnet - elastic net regression

## STEP 1: READ IN DATA
# read in data
info=readRDS('mydata.Rds')#This dataset include species characters and CpGs ##but what specifically?
#read in trained clocks
glmnet.csv=c('clock1.csv',
             'clock2.csv',
             'clock3.csv')#S3.1 to S3.3 tables in SupplmentaryData for Lu paper
#create lists for each clock
beta.name=c('beta_clock1','beta_clock2','beta_clock3')
y.name=c('Y.pred1','Y.pred2','Y.pred3')
age.name=c('DNAmAgeClock1','DNAmAgeClock2','DNAmAgeClock3') #predicted age 

## STEP 2: CONVERT AGE DATA TO CLOCK FORMATS
# CLOCK 1: basic clock
# directly regresses log-transformed chronological age on DNA methylation levels

# no transformations necessary ##..log?

# CLOCK 2: universal relative age clock
# defines individual age relative to the maximum lifespan of its species
# generating relative age estiamtes between 0 and 1
F2_antitrans_clock2<-function(y,y.maxAge,y.gestation,const=1){
  x0=const*exp(-exp(-1*y))
  x1=x0*(y.maxAge+y.gestation)
  x=x1-y.gestation
  x
}

# CLOCK 3: universal log-linear transformed age clock
# substitutes maximum lifespan for age at sexual maturity and gestation time
F1_logli <- function(age1, m1, m2 = m1, c1=1){
  ifelse(age1 >= m1, (age1-m1)/m2 , c1*log((age1-m1)/m2/c1 +1) )
}
#RelativeAdultAge
F2_revtrsf_clock3 <- function(y.pred, m1, m2 = m1, c1=1){
  ifelse(y.pred<0, (exp(y.pred/c1)-1)*m2*c1 + m1, y.pred*m2+m1 )
}
# The `loglifn` function shows how to calculate m1 for the transformation
# It is the `a_Logli` in the function
F3_loglifn = function(dat1,b1=1,max_tage = 4,c1=5, c2 = 0.38, c0=0){
  n=nrow(dat1)
  
  age1 = (dat1$maxAgeCaesar+dat1$GestationTimeInYears)/(dat1$averagedMaturity.yrs+dat1$GestationTimeInYears)
  
  a1 = age1/(1+max_tage)
  dat1$a1_Logli = a1 #x/m1 in manuscript
  
  a2 = (dat1$GestationTimeInYears + c0)/(dat1$averagedMaturity.yrs)
  dat1$a_Logli = a_Logli = c1*a2^c2
  #m=5*(G/ASM)^0.38 from regression analysis/formula(7)
  
  x = dat1$Age + dat1$GestationTimeInYears
  t2 = dat1$averagedMaturity.yrs*b1 + dat1$GestationTimeInYears
  x2 = x/t2 #### log(x/t2)
  y = F1_logli(x2, a_Logli, a_Logli)
  
  dat1$LogliAge <- y
  return(dat1)
}
# renaming? 
info$Intercept=1
MYMAX=1.3
info$HighmaxAgeCaesar=MYMAX*info$maxAgeCaesar
info$HighmaxAgeCaesar[info$SpeciesLatinName=='Homosapiens']=info$maxAgeCaesar[info$SpeciesLatinName=='Homo sapiens']
info$HighmaxAgeCaesar[info$SpeciesLatinName=='Musmusculus']=info$maxAgeCaesar[info$SpeciesLatinName=='Mus musculus']

## STEP 3: USE CLOCK TO PREDICT AGE
#predict RelativeAge#
for(k in 1:3){ # all 3 clocks
  glmnet=read.csv(glmnet.csv[k])
  glmnet$beta=glmnet[,beta.name[k]]
  glmnet$var[1]=ifelse(glmnet$var[1]=="(Intercept)",'Intercept',glmnet$var[1])
  temp=as.matrix(subset(info,select=as.character(glmnet$var)))
  info[,y.name[k]]=as.numeric(as.matrix(subset(info,select=as.character(glmnet$var)))%*%glmnet$beta)
}
#(1) Clock 1
info[,age.name[1]]=exp(info[,y.name[k]])-2
#(2) Clock 2
info$DNAmRelativeAge=exp(-exp(-1*info[,y.name[2]]))
info[,age.name[2]]=F2_antitrans_clock2(info[,y.name[2]],info$HighmaxAgeCaesar,info$GestationTimeInYears,const=1)
#(3) Clock 3
info=F3_loglifn(info)#to compute m estimate for tuning point in the log-linear transformation
info$m1=info$a_Logli
info$DNAmRelativeAdultAge=F2_revtrsf_clock3(info[,y.name[3]], info$m1)
info[,age.name[3]]<-
  info$DNAmRelativeAdultAge *(info$averagedMaturity.yrs + info$GestationTimeInYears) -
  info$GestationTimeInYears
#
#final output
#
output=subset(info,select=c('SampleID','Age','DNAmRelativeAge','DNAmRelativeAdultAge',age.name))