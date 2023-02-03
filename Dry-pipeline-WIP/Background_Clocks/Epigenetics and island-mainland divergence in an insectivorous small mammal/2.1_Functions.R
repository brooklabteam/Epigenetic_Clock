##R Script for helper functions for building Epigenetic Clocks
##Author: Joseph A. Zoller (from Steve Horvath's lab)

############################################
##### Defining simple helper functions #####
############################################
matches.any_in_gene_list <- Vectorize(function(x, genes.vec, split=";") {
  return(sum(match(genes.vec,strsplit(x,split=split)[[1]], nomatch=0)) > 0)
}, vectorize.args="x", USE.NAMES=F)
starts.any_in_gene_list <- Vectorize(function(x, genes.vec, split=";") {
  return(sum(sapply(sapply(genes.vec,grepstart,strsplit(x,split=split)[[1]]),length)) > 0)
}, vectorize.args="x", USE.NAMES=F)
grepstart <- function(prefix, x) {return(which(startsWith(x, prefix)))}
contains.any_in_gene_list <- Vectorize(function(x, genes.vec, split=";") {
  return(sum(sapply(sapply(genes.vec,grep,strsplit(x,split=split)[[1]]),length)) > 0)
}, vectorize.args="x", USE.NAMES=F)
contains.cg <- Vectorize(function(x) {return(length(grep("cg",x)) > 0)}, USE.NAMES=F)
get.folder_full_names <- Vectorize(function(x, folders.vec) {
  return(grep(x,folders.vec,value=T))
}, vectorize.args="x", USE.NAMES=F) #used in Primate scripts
panel.cor <- function(x, y, cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits=3)
  txt <- paste0("r = ", r)
  text(0.5, 0.5, txt, cex=cex.cor)
}
#' @export
rbind_merge <- function(x, y) {
  x.diff <- setdiff(colnames(x), colnames(y))
  y.diff <- setdiff(colnames(y), colnames(x))
  x[, as.character(y.diff)] <- NA
  y[, as.character(x.diff)] <- NA
  return(rbind(x, y))
}
#' @export
rbind_complete <- function(x, y) {
  xy.intersect <- intersect(colnames(x), colnames(y))
  x <- x[, as.character(xy.intersect)]
  y <- y[, as.character(xy.intersect)]
  return(rbind(x, y))
}
#' @export
get.square_limits <- function(vec_x, vec_y) {
  if (length(which(!is.na(vec_x)))==0) {return(NULL)}
  if (length(which(!is.na(vec_y)))==0) {return(NULL)}
  l_lim = min(1.03*min(vec_x,na.rm=T)-0.03*max(vec_x,na.rm=T),
              1.03*min(vec_y,na.rm=T)-0.03*max(vec_y,na.rm=T))
  u_lim = max(1.03*max(vec_x,na.rm=T)-0.03*min(vec_x,na.rm=T),
              1.03*max(vec_y,na.rm=T)-0.03*min(vec_y,na.rm=T))
  lim = c(l_lim, u_lim)
  return(lim)
}
#' @export
get.limits <- function(vec) {
  if (length(which(!is.na(vec)))==0) {return(NULL)}
  l_lim = 1.03*min(vec,na.rm=T)-0.03*max(vec,na.rm=T)
  u_lim = 1.03*max(vec,na.rm=T)-0.03*min(vec,na.rm=T)
  lim = c(l_lim, u_lim)
  return(lim)
}
# col2hsv <- Vectorize(function(col) {
#   return(rgb2hsv(r=col2rgb(col)[1,1], g=col2rgb(col)[2,1], b=col2rgb(col)[3,1]))
# })
#' @export
getSpeciesTissue <- Vectorize(function(sln, tissue, ignore=c("Homo sapiens")) {
  if (!sln %in% ignore) {return(paste(sln, tissue, sep="."))}
  else {return(paste(sln))}
})


###########################################
##### Defining Transformations of Age #####
###########################################
#' Identity Transformation
#'
#' @param x A numeric, representing age
#' @return The untransformed value
#' @examples
#' fun_identity(1)
#' @export
fun_identity <- function(x, ...) x
#' @export
fun_log_trans <- function(x, ...) log(x+1)
#' @export
fun_log_inv <- function(y, ...) exp(y)-1
#' @export
fun_sqrt_trans <- function(x, ...) sqrt(x+1)
#' @export
fun_sqrt_inv <- function(y, ...) y^2-1
### Applies the RelativeAge transformation to the input vector x
#' @export
fun_rel1_trans <- Vectorize(function(x, gestation, max) {
  if (is.na(x) | is.na(max)) {return(NA)}
  y <- x/max
  return(y)
})
### Applies the inverse RelativeAge transformation to the input vector y
#' @export
fun_rel1_inv <- Vectorize(function(y, gestation, max) {
  if (is.na(y) | is.na(max)) {return(NA)}
  x <- max*y
  return(x)
})
### Applies the LLin transformation to the input vector x
#' @export
fun_llin_trans <- Vectorize(function(x, maturity, max) {
  if (is.na(x) | is.na(maturity) | is.na(max)) {return(NA)}
  k <- fun_llin_trans_param(maturity, max)
  if (is.na(k)) {return(NA)}
  y <- 0
  if (x < maturity) {y = (maturity+k)*log((x+k)/(maturity+k))}
  else {y = x-maturity}
  return(y)
})
### Applies the inverse LLin transformation to the input vector y
#' @export
fun_llin_inv <- Vectorize(function(y, maturity, max) {
  if (is.na(y) | is.na(maturity) | is.na(max)) {return(NA)}
  k <- fun_llin_trans_param(maturity, max)
  if (is.na(k)) {return(NA)}
  x <- 0
  if (y < 0) {x = (maturity+k)*exp(y/(maturity+k))-k}
  else {x = y+maturity}
  return(x)
})
### Calculates the parameter k using the values of maturity and max
#' @export
fun_llin_trans_param <- Vectorize(function(maturity, max) {
  if (is.na(maturity) | is.na(max)) {return(NA)}
  N = max/(maturity)
  # max > 3*maturity+2 for all but 2 species in data (Maui's dolphin, Star-nosed mole)
  if (is.na(N)) {print("max or maturity is NA")}
  if (max <= 3*maturity+2) {
    print("k does not exist; max/maturity condition is not met")
    return(NA)
  } else {
    k_lower = 1 + maturity/exp(N-1)
    k_upper = 1 + 2*(maturity+1)^2/(maturity*(N-3)-2)
    k_interval = c(k_lower, k_upper)
    k = uniroot(function(k) {(k+maturity)*log((k+maturity)/(k-1))-(max-maturity)/2},k_interval)$root
    return(k)
  }
})
### Applies the LLin2 transformation to the input vector x
#' @export
fun_llin2_trans <- Vectorize(function(x, maturity, ...) {
  if (is.na(x) | is.na(maturity)) {return(NA)}
  k <- 1.5
  y <- 0
  if (x < maturity) {y = (maturity+k)*log((x+k)/(maturity+k))}
  else {y = x-maturity}
  return(y)
})
### Applies the inverse LLin2 transformation to the input vector y
#' @export
fun_llin2_inv <- Vectorize(function(y, maturity, ...) {
  if (is.na(y) | is.na(maturity)) {return(NA)}
  k <- 1.5
  x <- 0
  if (y < 0) {x = (maturity+k)*exp(y/(maturity+k))-k}
  else {x = y+maturity}
  return(x)
})
### Applies the LLin3 transformation to the input vector x
#' @export
fun_llin3_trans <- Vectorize(function(x, maturity, ...) {
  if (is.na(x) | is.na(maturity)) {return(NA)}
  k <- 1.5
  y <- 0
  if (x < maturity) {y = log((x+k)/(maturity+k))}
  else {y = (x-maturity)/(maturity+k)}
  return(y)
})
### Applies the inverse LLin3 transformation to the input vector y
#' @export
fun_llin3_inv <- Vectorize(function(y, maturity, ...) {
  if (is.na(y) | is.na(maturity)) {return(NA)}
  k <- 1.5
  x <- 0
  if (y < 0) {x = (maturity+k)*exp(y)-k}
  else {x = (maturity+k)*y+maturity}
  return(x)
})
### Applies the LLin4 transformation to the input vector x
#' @export
fun_llin4_trans <- Vectorize(function(x, maturity, gestation) {
  if (is.na(x) | is.na(maturity) | is.na(gestation)) {return(NA)}
  k <- 2*gestation
  y <- 0
  if (x <= -k) {return(NA)}
  if (x < maturity) {y = log((x+k)/(maturity+k))}
  else {y = (x-maturity)/(maturity+k)}
  return(y)
})
### Applies the inverse LLin4 transformation to the input vector y
#' @export
fun_llin4_inv <- Vectorize(function(y, maturity, gestation) {
  if (is.na(y) | is.na(maturity) | is.na(gestation)) {return(NA)}
  k <- 2*gestation
  x <- 0
  if (y < 0) {x = (maturity+k)*exp(y)-k}
  else {x = (maturity+k)*y+maturity}
  return(x)
})
### Applies the LLin5 transformation to the input vector x
#' @export
fun_llin5_trans <- Vectorize(function(x, maturity, ...) {
  if (is.na(x) | is.na(maturity)) {return(NA)}
  k <- 0.5*maturity
  y <- 0
  if (x <= -k) {return(NA)}
  if (x < maturity) {y = log((x+k)/(maturity+k))}
  else {y = (x-maturity)/(maturity+k)}
  return(y)
})
### Applies the inverse LLin5 transformation to the input vector y
#' @export
fun_llin5_inv <- Vectorize(function(y, maturity, ...) {
  if (is.na(y) | is.na(maturity)) {return(NA)}
  k <- 0.5*maturity
  x <- 0
  if (y < 0) {x = (maturity+k)*exp(y)-k}
  else {x = (maturity+k)*y+maturity}
  return(x)
})

#####################################
##### Defining helper functions #####
#####################################
### Loads an RData file
### Returns: A list
#' @export
loadRData <- function(file.rdata) {
  load(file.rdata)
  get(ls()[ls() != "file.rdata"])
}

### Transposes data frame where colname and first column contain labels
### Returns: a data frame with original colnames as the first column
###          and all original rows as remaining columns
#' @export
transposeDat <- function(dat_tp, VAR.label) {
  #VAR.label=String, name of the new first column containing original colnames
  dat_tp <- as.data.frame(dat_tp)
  dat <- dat_tp %>% dplyr::select(-1) %>% t() %>% as.data.frame()
  colnames(dat)=dat_tp[,1]
  dat <- cbind(V1=rownames(dat), dat)
  rownames(dat)=NULL
  dat$V1 <- as.character(dat$V1)
  colnames(dat) <- c(VAR.label, colnames(dat)[-1])
  return(dat)
}

### Loads dataset from Rdata file, then runs transposeDat
#' @export
transposeDat_rdata <- function(dat_tp.rdata, VAR.label) {
  dat_tp <- loadRData(dat_tp.rdata)
  #dat_tp <- as.data.frame(dat_tp)
  dat <- transposeDat(dat_tp, VAR.label)
  return(dat)
}

### Aligns Info and Covariate Datasets according to an ID variable
### Covariate Dataset must have the ID variable and all covariates as columns
### Returns: a 2-element list of the new info and covariate datasets, where the new
###          covariate dataset is a matrix that no longer contains the ID variable
#' @export
alignInfoDat <- function(info,dat,ALIGNVAR.info,ALIGNVAR.dat) {
  #ALIGNVAR.info=String, name of variable in info by which to align
  #ALIGNVAR.dat=String, name of variable in dat by which to align
  if (!(ALIGNVAR.info %in% colnames(info))) {stop(paste("ALIGNVAR.info=",ALIGNVAR.info," not found",sep=""))}
  if (!(ALIGNVAR.dat %in% colnames(dat))) {stop(paste("ALIGNVAR.dat=",ALIGNVAR.dat," not found",sep=""))}
  info_alignvar=info[,ALIGNVAR.info]
  dat_alignvar=dat[,ALIGNVAR.dat]
  matchid=match(info_alignvar,dat_alignvar)
  info.new=info[which(!is.na(matchid)),]
  matchid.new=matchid[which(!is.na(matchid))]
  dat.new=dat[matchid.new,]
  if (dim(info.new)[1] != dim(dat.new)[1]) {stop ("check")}
  ys=info.new
  xs=as.matrix(dat.new[,which(colnames(dat.new) != ALIGNVAR.dat)])
  xs=xs[,which(contains.cg(colnames(xs)))]
  if (is.vector(xs)) {xs=t(as.matrix(xs))}
  yxs=list(ys,xs)
  return(yxs)
}

### Partitions the dataset (randomly) into a training and a testing set
### Returns: a 4-four element list containing xs.train,ys.train,xs.test,ys.test
#' @export
splitDataTrainTest <- function(xs,ys,prop.test=0.3) {
  permu=sample(dim(xs)[1])
  xs=xs[permu,]
  ys=ys[permu,]

  N=dim(xs)[1]
  Nt=N-floor(N*prop.test)
  N1=Nt+1
  xs.train=xs[1:Nt,]
  ys.train=ys[1:Nt,]
  xs.test=xs[N1:N,]
  ys.test=ys[N1:N,]

  return(list(xs.train,ys.train,xs.test,ys.test))
}

### Computes rank p-values for each CpG based on single CpG-Age correlations, stratified
###          by a given variable (SPECVAR)
### Finds the most significant CpG sites (based on rank p-values) for both positive and negative correlation
### Returns: a character vector, cpg_names_rankp_sig, of length 2*halfsize,
###          containing the names of the significant CpG sites
#' @export
selectCPGRankPval <- function(info,dat,OUTVAR="Age",SPECVAR="SpeciesLatinName",halfsize=2000,fun_trans=fun_identity,fun_VAR1=NULL,fun_VAR2=NULL) {
  #OUTVAR=String, name of dependent variable in info (default is "Age")
  #SPECVAR=String, name of stratifying variable (FACTOR) in info (default is "SpeciesLatinName")
  #halfsize=Integer, number of CpGs with most significant positive/negative correlation with Age to select (default is 1000)
  rows.keep <- which(!is.na(fun_trans(info[,OUTVAR],info[,fun_VAR1],info[,fun_VAR2])))
  if (length(rows.keep) < nrow(info)) {
    print(paste(nrow(info)-length(rows.keep),"rows were ignored"))
    info <- info[rows.keep,]
    info[,SPECVAR] <- factor(info[,SPECVAR])
  }
  info.species=info[,SPECVAR]
  ### Generating correlation table
  cpg_names <- colnames(dat)[which(contains.cg(colnames(dat)))]
  spec_names <- levels(info.species)
  corr_cpgxspec <- data.frame(matrix(NA,nrow=length(cpg_names),ncol=length(spec_names)))
  rownames(corr_cpgxspec) <- cpg_names
  colnames(corr_cpgxspec) <- spec_names
  for (spec.num in 1:length(spec_names)) {
    indcs <- which(as.numeric(info.species) == spec.num)
    if (length(indcs) == 0) {corr_cpgxspec[,spec.num] <- NULL}
    else {
      info_subspec <- info[indcs,]
      yxs_subspec.list <- alignInfoDat(info_subspec,dat,"Basename","Basename")
      ys_subspec <- yxs_subspec.list[[1]]
      xs_subspec <- yxs_subspec.list[[2]]
      ys_subspec.outcome=ys_subspec[,OUTVAR]
      ys_subspec.funvar1=ys_subspec[,fun_VAR1]
      ys_subspec.funvar2=ys_subspec[,fun_VAR2]
      corr_cpgxspec[,spec.num] <- cor(fun_trans(ys_subspec.outcome,ys_subspec.funvar1,ys_subspec.funvar2), xs_subspec)[1,]
    }
  }
  ### Generating rank p-values and selecting CpGs
  rankp_obj <- rankPvalue(corr_cpgxspec)
  cpg_names_rankp_hi <- rownames(rankp_obj[order(rankp_obj$pValueHighRank)[1:halfsize],])
  cpg_names_rankp_lo <- rownames(rankp_obj[order(rankp_obj$pValueLowRank)[1:halfsize],])
  cpg_names_rankp_sig <- sort(c(cpg_names_rankp_hi,cpg_names_rankp_lo))
  return(cpg_names_rankp_sig)
}

### Computes S sets of rank p-values for each CpG based on single CpG-Age correlations, stratified
###          by a given variable (SPECVAR), where in each set, one level of SPECVAR is removed
### Finds the most significant CpG sites (based on rank p-values) for both positive and negative correlation
###          within each of the S sets
### Returns: a list of character vectors (of length 2*halfsize each), cpg_names_rankp_sig.list,
###          of length S=levels(SPECVAR), containing the names of the significant
###          CpG sites for each leave-one-SPECVAR-out group
#' @export
selectCPGRankPvalLOSO <- function(info,dat,OUTVAR="Age",SPECVAR="SpeciesLatinName",halfsize=2000,fun_trans=fun_identity,fun_VAR1=NULL,fun_VAR2=NULL) {
  #OUTVAR=String, name of dependent variable in info (default is "Age")
  #SPECVAR=String, name of stratifying variable (FACTOR) in info (default is "SpeciesLatinName")
  #halfsize=Integer, number of CpGs with most significant positive/negative correlation with Age to select (default is 1000)
  rows.keep <- which(!is.na(fun_trans(info[,OUTVAR],info[,fun_VAR1],info[,fun_VAR2])))
  if (length(rows.keep) < nrow(info)) {
    print(paste(nrow(info)-length(rows.keep),"rows were ignored"))
    info <- info[rows.keep,]
    info[,SPECVAR] <- factor(info[,SPECVAR])
  }
  info.species=info[,SPECVAR]
  ### Generating correlation table
  ##xs=dat[,which(contains.cg(colnames(dat)))]
  cpg_names <- colnames(dat)[which(contains.cg(colnames(dat)))]
  spec_names <- levels(info.species)
  corr_cpgxspec <- data.frame(matrix(NA,nrow=length(cpg_names),ncol=length(spec_names)))
  rownames(corr_cpgxspec) <- cpg_names
  colnames(corr_cpgxspec) <- spec_names
  for (spec.num in 1:length(spec_names)) {
    indcs <- which(as.numeric(info.species) == spec.num)
    if (length(indcs) == 0) {corr_cpgxspec[,spec.num] <- NULL}
    else {
      info_subspec <- info[indcs,]
      yxs_subspec.list <- alignInfoDat(info_subspec,dat,"Basename","Basename")
      ys_subspec <- yxs_subspec.list[[1]]
      xs_subspec <- yxs_subspec.list[[2]]
      ys_subspec.outcome=ys_subspec[,OUTVAR]
      ys_subspec.funvar1=ys_subspec[,fun_VAR1]
      ys_subspec.funvar2=ys_subspec[,fun_VAR2]
      corr_cpgxspec[,spec.num] <- cor(fun_trans(ys_subspec.outcome,ys_subspec.funvar1,ys_subspec.funvar2), xs_subspec)[1,]
    }
  }
  ### Generating rank p-values and selecting CpGs
  cpg_names_rankp_sig.list <- vector("list",length(levels(info.species)))
  for (spec.num in 1:length(levels(info.species))) {
    rankp_obj.num <- rankPvalue(corr_cpgxspec[,-spec.num])
    cpg_names_rankp_hi.num <- rownames(rankp_obj.num[order(rankp_obj.num$pValueHighRank)[1:halfsize],])
    cpg_names_rankp_lo.num <- rownames(rankp_obj.num[order(rankp_obj.num$pValueLowRank)[1:halfsize],])
    cpg_names_rankp_sig.num <- sort(c(cpg_names_rankp_hi.num,cpg_names_rankp_lo.num))
    cpg_names_rankp_sig.list[[spec.num]] <- cpg_names_rankp_sig.num
  }
  return(cpg_names_rankp_sig.list)
}

### Generates elastic net by training on a separate dataset
### Saves the coefficients, the info with predicted values, and performance plots
### Returns: a data.frame, ys.output, containing ys.train+ys.test with 4 more columns, being
###          "train" (logical), PREDVAR, RESVAR, and RESinTestVAR
#' @export
saveNetTrained <- function(xs.train,ys.train,xs.test,ys.test,OUTVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinTestVAR,ALPHA=0.5,NFOLD=10,fun_trans=fun_identity,fun_inv=fun_identity,fun_VAR1=NULL,fun_VAR2=NULL,loglambda.seq=NULL,parallel=FALSE) {
  #OUTVAR=String, name of dependent variable in ys.train and ys.test
  #out.csv <- table of coefficients
  #output.csv <- ys.output
  #out.png <- image of model diagnostics and results plot
  #out.png.title=String, title used in model diagnostics and results plot
  #PREDVAR=String, desired name of the predicted value for the dependent variable
  #RESVAR=String, desired name of the predicted value's residual on the outcome
  #RESinTestVAR=String, desired name of the predicted value's residual on the outcome, within the testing set only
  #ALPHA=1 for lasso, ALPHA=0 for ridge
  #fun_trans=Function, transformation applied to OUTVAR before fitting the model
  #fun_inv=Function, topological inverse of fun_trans
  #fun_VAR1=String, name of variable to be used as the 1st parameter of fun_inv (default is NULL)
  #fun_VAR2=String, name of variable to be used as the 2nd parameter of fun_inv (default is NULL)
  rows.keep_train <- which(!is.na(fun_trans(ys.train[,OUTVAR],ys.train[,fun_VAR1],ys.train[,fun_VAR2])))
  if (length(rows.keep_train) < nrow(ys.train)) {
    print(paste(nrow(ys.train)-length(rows.keep_train),"rows were ignored"))
  }
  ys.train <- ys.train[rows.keep_train,]
  xs.train <- xs.train[rows.keep_train,]
  if (nrow(ys.test)==0) {
    ys.test <- NULL
    xs.test <- NULL
  }
  ys.train_outcome=ys.train[,OUTVAR]
  ys.test_outcome=ys.test[,OUTVAR]
  ys.train_funvar1=ys.train[,fun_VAR1] # 0 columns iff VAR=NULL
  ys.train_funvar2=ys.train[,fun_VAR2]
  ys.test_funvar1=ys.test[,fun_VAR1]
  ys.test_funvar2=ys.test[,fun_VAR2]
  if (is.null(fun_VAR1)) {
    ys.train_funvar1=NULL
    ys.test_funvar1=NULL
  }
  if (is.null(fun_VAR2)) {
    ys.train_funvar2=NULL
    ys.test_funvar2=NULL
  }
  if (parallel) {
    library(doMC)
    registerDoMC()
  }

  if (is.null(loglambda.seq)) {lambda=NULL} else {lambda=exp(loglambda.seq)}
  glmnet.Training.CV = cv.glmnet(xs.train,fun_trans(ys.train_outcome,ys.train_funvar1,ys.train_funvar2),
                                 nfolds=NFOLD,alpha=ALPHA,family="gaussian",type.measure="mse",lambda=lambda,parallel=parallel)
  lambda.glmnet.Training = glmnet.Training.CV$lambda.1se
  glmnet.Training = glmnet(xs.train,fun_trans(ys.train_outcome,ys.train_funvar1,ys.train_funvar2),
                           alpha=ALPHA,family="gaussian",lambda=lambda,nlambda=100)
  # check glmnet.Training with print.glmnet
  # a three-column matrix:
  # (1)Df (the number of nonzero coefficients when using pure lasso)
  # (2)%dev (the percent deviance explained, relative to the null deviance)
  #    %dev = pseudoR2 / pseudoR2_max ~= pseudoR2
  # (3)Lambda (the corresponding value of lambda, the sparsity penalty)
  # print(glmnet.Training)
  print("Debug checkpoint 1")

  if (!is.null(out.csv)) {
    glmnet.final=data.frame(as.matrix(coef(glmnet.Training,s=lambda.glmnet.Training)))
    names(glmnet.final)='beta'
    glmnet.final$var=rownames(glmnet.final)
    glmnet.final=subset(glmnet.final,select=c(var,beta))
    glmnet.final1=subset(glmnet.final,abs(beta)>0)
    write.table(glmnet.final1,out.csv,sep=',',row.names=F,quote=F)
  }
  print("Debug checkpoint 2")

  train=data.frame(ys.train,train=1)
  ys.train_prediction=fun_inv(as.numeric(predict(glmnet.Training,xs.train,type='response',s=lambda.glmnet.Training)),
                              ys.train_funvar1,ys.train_funvar2)
  train[,PREDVAR]=ys.train_prediction
  attr(train[,PREDVAR],'dimnames')<-NULL
  train[,RESinTestVAR]=NA
  attr(train[,RESinTestVAR],'dimnames')<-NULL
  if (!is.null(ys.test)) {
    test=data.frame(ys.test,train=0)
    ys.test_prediction=fun_inv(as.numeric(predict(glmnet.Training,xs.test,type='response',s=lambda.glmnet.Training)),
                               ys.test_funvar1,ys.test_funvar2)
    test[,PREDVAR]=ys.test_prediction
    attr(test[,PREDVAR],'dimnames')<-NULL
    test[,RESinTestVAR]=NA
    if (nrow(ys.test) > 1 && length(which(!is.na(ys.test_outcome) & !is.na(ys.test_prediction))) > 2) {
      test[which(!is.na(ys.test_outcome) & !is.na(ys.test_prediction)),RESinTestVAR]=residuals(lm(ys.test_prediction~ys.test_outcome))
    }
    attr(test[,RESinTestVAR],'dimnames')<-NULL
    ys.output=rbind_merge(train,test)
  }
  if (is.null(ys.test)) {
    ys.output <- train
  }
  ys.output[,RESVAR]=NA
  if (!is.null(PREDVAR)) {
    ys.output[which(!is.na(ys.output[,OUTVAR]) & !is.na(ys.output[,PREDVAR])),RESVAR]=residuals(lm(ys.output[,PREDVAR]~ys.output[,OUTVAR]))
  }
  attr(ys.output[,RESVAR],'dimnames')<-NULL
  if (!is.null(output.csv)) {
    write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
  }
  print("Debug checkpoint 3")

  if (is.null(out.png)) {
    return(ys.output)
  }
  png(out.png,width=9,height=10,units='in',res=300)
  par(mfrow=c(2,2))
  par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,0))
  plot(glmnet.Training.CV,
       main=paste0('Optimal Lambda=',round(lambda.glmnet.Training,digits=3)),
       cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
  #title('A',adj=0,font=2,cex.main=1.4)
  plot(glmnet.Training,xvar='lambda',label=T,
       main=paste0(nrow(glmnet.final1),' CpGs with optimal lambda'),
       cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
  #title('B',adj=0,font=2,cex.main=1.4)
  PANEL_str <- "Training data"
  N <- length(which(!is.na(ys.train_outcome) & !is.na(ys.train_prediction)))
  if (var(ys.train_outcome,na.rm=T) > 0) {
    COR <- cor(ys.train_prediction,ys.train_outcome,
               use='pairwise.complete.obs')
    COR_str <- paste0("cor=",signif(COR,2))
  } else {
    COR_str <- ''
  }
  plot(y=ys.train_prediction,x=ys.train_outcome,
       main=paste0(PANEL_str,' (N=',N,')'),
       ylab=PREDVAR,xlab=OUTVAR,
       cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
  if (var(ys.train_outcome, na.rm=T) > 0 && !is.na(var(ys.train_prediction, na.rm=T))) {
    abline(lm(ys.train_prediction~ys.train_outcome))
  }
  abline(0,1,lty="dashed")
  title(COR_str,outer=F,line=0.4,cex.main=1.5)
  #title('C',adj=0,font=2,cex.main=1.4)
  PANEL_str <- "Testing data"
  if (!is.null(ys.test)) {
    N <- length(which(!is.na(ys.test_outcome) & !is.na(ys.test_prediction)))
  } else {
    N <- 0
  }
  if (N > 2) {
    if (var(ys.test_outcome,na.rm=T) > 0) {
      COR <- cor(ys.test_prediction,ys.test_outcome,
                 use='pairwise.complete.obs')
      COR_str <- paste0("cor=",signif(COR,2))
    } else {
      COR_str <- ''
    }
    plot(y=ys.test_prediction,x=ys.test_outcome,
         main=paste0(PANEL_str,' (N=',N,')'),
         ylab=PREDVAR,xlab=OUTVAR,
         cex.main=1.5,cex.lab=1.5,cex.axis=1.5)
    if (var(ys.test_outcome, na.rm=T) > 0 && !is.na(var(ys.test_prediction, na.rm=T))) {
      abline(lm(ys.test_prediction~ys.test_outcome))
    }
    abline(0,1,lty="dashed")
    title(COR_str,outer=F,line=0.4,cex.main=1.5)
    #title('D',adj=0,font=2,cex.main=1.4)
  }
  title(paste0(out.png.title,'\n',ncol(xs.train),' CpG sites'),outer=T,line=-1,cex.main=1.5)
  dev.off()
  return(ys.output)
}

### Generates a set of elastic nets by training on a separate dataset, stratified by a variable
### Predicts outcome on a new data-info pair, using the set of elastic nets
### Saves the matrix of coefficients, the info with predicted values and the correlation panel plot
### Returns: a data.frame, ys.output, containing ys.train+ys.test with 4 more columns, being
###          "train" (logical), PREDVAR, RESVAR, and RESinTestVAR
#' @export
saveNetTrainedStrat <- function(xs.train,ys.train,xs.test,ys.test,OUTVAR,STRATVAR,out.csv,output.csv,out.png,out.png.title,PREDVAR,RESVAR,RESinTestVAR,ALPHA=0.5,NFOLD=10,fun_trans=fun_identity,fun_inv=fun_identity,fun_VAR1=NULL,fun_VAR2=NULL,mfrow=c(3,3),width=13,height=14,COLVAR=NULL,show.original=T,loglambda.seq=NULL) {
  #OUTVAR=String, name of dependent variable in ys.train and ys.test
  #STRATVAR=String, name of variable (FACTOR) in ys.train and ys.test by which to stratify the elastic nets
  #         NOTE: CANNOT be NULL
  #out.csv <- table of coefficients
  #output.csv <- ys.output
  #out.png <- image of model diagnostics and results plot
  #out.png.title=String, title used in model diagnostics and results plot
  #PREDVAR=String, desired name of the predicted value for the dependent variable
  #RESVAR=String, desired name of the predicted value's residual on the outcome
  #RESinTestVAR=String, desired name of the predicted value's residual on the outcome, within the testing set only
  #ALPHA=1 for lasso, ALPHA=0 for ridge
  #fun_trans=Function, transformation applied to OUTVAR before fitting the model
  #fun_inv=Function, topological inverse of fun_trans
  #fun_VAR1=String, name of variable to be used as the 1st parameter of fun_inv (default is NULL)
  #fun_VAR2=String, name of variable to be used as the 2nd parameter of fun_inv (default is NULL)
  #mfrow=2-element Vector, desired rxc dimensions of the panel plot
  #COLVAR=String, name of variable (FACTOR) in ys.output by which to color points
  #       NOTE: CAN be NULL (default is NULL, uses STRATVAR)
  #show.original=Logical, whether to show the original full-data plot in the first panel (default is TRUE)
  if (is.null(STRATVAR)) {
    stop('STRATVAR Cannot be NULL')
  }
  if (is.null(COLVAR)) {
    COLVAR = STRATVAR
  }

  ys.output <- data.frame()
  glmnet.final <- data.frame(var=colnames(xs.train))
  out.strat.csv='temp_saveNetTrainedStrat.csv' # temporary file where coefficients are saved, per strata
  for (strat in levels(ys.train[,STRATVAR])) {
    print(paste("Debug checkpoint: stratum",strat))
    xs.train.strat <- xs.train[which(ys.train[,STRATVAR] == strat),]
    ys.train.strat <- ys.train[which(ys.train[,STRATVAR] == strat),] %>%
      dplyr::mutate(train = 1)
    xs.test.strat <- xs.test[which(ys.test[,STRATVAR] == strat),]
    ys.test.strat <- ys.test[which(ys.test[,STRATVAR] == strat),] %>%
      dplyr::mutate(train = 0)
    if (length(which(ys.train[,STRATVAR] == strat)) == 1) {
      stop(paste("Training set is too small in",strat))
    }
    if (length(which(ys.test[,STRATVAR] == strat)) == 1) {
      xs.test.strat <- t(as.matrix(xs.test[which(ys.test[,STRATVAR] == strat),]))
    }
    ys.output.strat <- saveNetTrained(xs.train.strat,ys.train.strat,xs.test.strat,ys.test.strat,
                                      OUTVAR,out.csv=out.strat.csv,output.csv=NULL,out.png=NULL,out.png.title=NULL,PREDVAR=PREDVAR,RESVAR=RESVAR,RESinTestVAR=RESinTestVAR,
                                      ALPHA=ALPHA,NFOLD=NFOLD,fun_trans=fun_trans,fun_inv=fun_inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,loglambda.seq=loglambda.seq)
    ys.output <- rbind(ys.output, ys.output.strat)
    glmnet.strat <- read.csv(out.strat.csv)
    file.remove(out.strat.csv)
    colnames(glmnet.strat) <- c("var",strat)
    glmnet.final <- base::merge(glmnet.final, glmnet.strat, "var", all=T, sort=F)
  }
  glmnet.final$num_shared <- rowSums(sapply(as.data.frame(!is.na(glmnet.final)), as.numeric))-1
  glmnet.final <- dplyr::filter(glmnet.final, num_shared >= 1)
  ys.output[,STRATVAR] <- factor(ys.output[,STRATVAR])

  if (!is.null(out.csv)) {
    write.table(glmnet.final,out.csv,sep=',',row.names=F,quote=F)
  }
  if (!is.null(output.csv)) {
    write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
  }
  if (!is.null(out.png)) {
    TITLE_str=paste0(out.png.title,'\n')
    saveValidationPanelPlot(ys.output,OUTVAR,PREDVAR,PANELVAR=STRATVAR,out.png,TITLE_str,mfrow=mfrow,width=width,height=height,COLVAR=COLVAR,show.original=show.original)
  }

  return(ys.output)
}

### Predicts outcome on a new data-info pair, using a given table of coefficients
### Saves the info with predicted values and the correlation plot
### Returns: a data.frame, ys.output, containing ys.new with 2 more columns, being
###          PREDVAR and RESVAR
#' @export
saveNetNewTest <- function(in.valbeta,xs.new,ys.new,OUTVAR,output.csv,out.png,out.png.title,PREDVAR,RESVAR,fun_inv=fun_identity,fun_VAR1=NULL,fun_VAR2=NULL,COLVAR=NULL,show.legend=T,oma.right=9) {
  #in.valbeta <- table of coefficients to be used (column1=name, column2=value)
  #OUTVAR=String, name of dependent variable in ys.new
  #output.csv <- ys.output
  #out.png <- image of validation plot
  #out.png.title=String, title used in validation plot
  #PREDVAR=String, desired name of the predicted value for the dependent variable
  #RESVAR=String, desired name of the predicted value's residual on the outcome
  #fun_inv=Function, topological inverse of fun_trans
  #fun_VAR1=String, name of variable to be used as the 1st parameter of fun_inv (default is NULL)
  #fun_VAR2=String, name of variable to be used as the 2nd parameter of fun_inv (default is NULL)
  #COLVAR=String, name of variable (FACTOR) by which to color the points (default is NULL)
  if (nrow(ys.test)==0) {
    ys.test <- NULL
    xs.test <- NULL
  }
  ys.new_outcome=ys.new[,OUTVAR]
  ys.new_funvar1=ys.new[,fun_VAR1] # 0 columns iff VAR=NULL
  ys.new_funvar2=ys.new[,fun_VAR2]
  if (is.null(fun_VAR1)) {ys.funvar1=NULL}
  if (is.null(fun_VAR2)) {ys.funvar2=NULL}

  ones.column=data.frame(V1=rep(1,nrow(xs.new)))
  colnames(ones.column) <- as.character(in.valbeta[1,1])

  xs.new <- t(na.omit(t(xs.new))) #only keep probes with values in all samples
  in.valbeta_trimmed <- in.valbeta[sort(c(1,which(as.character(in.valbeta[,1]) %in% colnames(xs.new)))),]
  num_beta.not_in_New <- nrow(in.valbeta) - nrow(in.valbeta_trimmed)
  xs.new_trimmed <- xs.new[,which(colnames(xs.new) %in% as.character(in.valbeta_trimmed[,1]))]
  xs.new_trimmed = cbind(ones.column, xs.new_trimmed)
  xs.new_trimmed <- dplyr::select(xs.new_trimmed, as.character(in.valbeta_trimmed[,1])) #sort xs.new to match in.valbeta
  ys.new_prediction.raw = as.vector(as.matrix(xs.new_trimmed) %*% in.valbeta_trimmed[,2])
  if (length(ys.new_prediction.raw) > 0) {
    ys.new_prediction = fun_inv(ys.new_prediction.raw,ys.new_funvar1,ys.new_funvar2)
  } else {
    ys.new_prediction = ys.new_prediction.raw
  }
  ys.output=ys.new
  ys.output[,PREDVAR]=ys.new_prediction
  attr(ys.output[,PREDVAR],'dimnames')<-NULL
  if (!is.null(RESVAR)) {
    ys.output[,RESVAR]=NA
    ys.output[which(!is.na(ys.new_outcome) & !is.na(ys.new_prediction)),RESVAR]=residuals(lm(ys.new_prediction~ys.new_outcome))
    attr(ys.output[,RESVAR],'dimnames')<-NULL
  }
  if (!is.null(output.csv)) {
    write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
  }
  if (!is.null(out.png)) {
    TITLE_str <- paste0(out.png.title,"\n",ncol(xs.new_trimmed)," CpG sites used, ",
                        num_beta.not_in_New," unused")
    if (show.legend) {
      width=9
      height=8
    } else {
      width=7
      height=8
    }
    saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str,width=width,height=height,show.legend=show.legend,oma.right=oma.right)
  } else {
    print(paste0(ncol(xs.new_trimmed)," CpG sites used, ",
                 num_beta.not_in_New," unused"))
  }

  return(ys.output)
}

### Estimates outcome on a new data-info pair, using LOO glmnet models
### Saves the info with predicted values and the correlation plot
### Returns: a data.frame, ys.output, containing ys with 3 more columns, being
###          PREDVAR, RESVAR, and "cpg_count_used"
#' @export
saveLOOEstimation <- function(xs,ys,OUTVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA=0.5,NFOLD=10,fun_trans=fun_identity,fun_inv=fun_identity,fun_VAR1=NULL,fun_VAR2=NULL,COLVAR=NULL,show.legend=T,oma.right=9,NUMVAR=NULL,xs.add_train=NULL,ys.add_train=NULL,loglambda.seq=NULL,parallel=FALSE) {
  #OUTVAR=String, name of dependent variable in ys
  #out.rdata <- list containing table of coefficients for each LOO model
  #output.csv <- ys.output
  #out.png <- image of validation plot
  #out.png.title=String, title used in validation plot
  #PREDVAR=String, desired name of the predicted value for the dependent variable
  #RESVAR=String, desired name of the predicted value's residual on the outcome
  #ALPHA=1 for lasso, ALPHA=0 for ridge
  #fun_trans=Function, transformation applied to OUTVAR before fitting the model
  #fun_inv=Function, topological inverse of fun_trans
  #fun_VAR1=String, name of variable to be used as the 1st parameter of fun_trans and fun_inv (default is NULL)
  #fun_VAR2=String, name of variable to be used as the 2nd parameter of fun_trans and fun_inv (default is NULL)
  #COLVAR=String, name of variable (FACTOR) by which to color the points (default is NULL)
  #NUMVAR=String, name of variable (FACTOR) by which to number the points (default is NULL)
  #xs.add_train/ys.add_train=Additional data to be included in the training set at every iteration
  rows.keep <- which(!is.na(fun_trans(ys[,OUTVAR],ys[,fun_VAR1],ys[,fun_VAR2])))
  if (length(rows.keep) < nrow(ys)) {
    print(paste(nrow(ys)-length(rows.keep),"rows were ignored"))
  }
  ys <- ys[rows.keep,]
  xs <- xs[rows.keep,]
  ys.outcome=ys[,OUTVAR]
  ys.funvar1=ys[,fun_VAR1] # 0 columns iff VAR=NULL
  ys.funvar2=ys[,fun_VAR2]
  if (is.null(fun_VAR1)) {ys.funvar1=NULL}
  if (is.null(fun_VAR2)) {ys.funvar2=NULL}
  if (!is.null(xs.add_train) & !is.null(ys.add_train)) {
    rows.keep <- which(!is.na(fun_trans(ys.add_train[,OUTVAR],ys.add_train[,fun_VAR1],ys.add_train[,fun_VAR2])))
    ys.add_train <- ys.add_train[rows.keep,]
    xs.add_train <- xs.add_train[rows.keep,]
    ys.add_train_outcome=ys.add_train[,OUTVAR]
    ys.add_train_funvar1=ys.add_train[,fun_VAR1] # 0 columns iff VAR=NULL
    ys.add_train_funvar2=ys.add_train[,fun_VAR2]
    if (is.null(fun_VAR1)) {ys.add_train_funvar1=NULL}
    if (is.null(fun_VAR2)) {ys.add_train_funvar2=NULL}
  }

  if (parallel) {
    library(doMC)
    registerDoMC()
  }

  if (is.null(loglambda.seq)) {lambda=NULL} else {lambda=exp(loglambda.seq)}
  glmnet.Training.CV = cv.glmnet(xs,fun_trans(ys.outcome,ys.funvar1,ys.funvar2),
                                 nfolds=NFOLD,alpha=ALPHA,family="gaussian",type.measure="mse",lambda=lambda,parallel=parallel)
  lambda.glmnet.Training = glmnet.Training.CV$lambda.1se
  print("Debug checkpoint 1")
  print(paste("log(Lambda) =",log(lambda.glmnet.Training)))

  ### Building LOO Clocks for each Sample
  ys.prediction <- 0*ys.outcome
  glmnet.final.list <- vector("list",length(ys.prediction))
  cpg.count <- 0*ys.outcome
  for (i in 1:length(ys.prediction)) {
    xs.trn <- xs[-i,]
    ys.trn_outcome <- ys.outcome[-i]
    ys.trn_funvar1 <- ys.funvar1[-i]
    ys.trn_funvar2 <- ys.funvar2[-i]
    if (!is.null(xs.add_train) & !is.null(ys.add_train)) {
      xs.trn <- rbind(xs.trn, xs.add_train)
      ys.trn_outcome <- c(ys.trn_outcome, ys.add_train_outcome)
      ys.trn_funvar1 <- c(ys.trn_funvar1, ys.add_train_funvar1)
      ys.trn_funvar2 <- c(ys.trn_funvar2, ys.add_train_funvar2)
    }
    glmnet.Training <- glmnet(xs.trn,fun_trans(ys.trn_outcome,ys.trn_funvar1,ys.trn_funvar2),
                              alpha=ALPHA,family="gaussian",lambda=lambda.glmnet.Training)
    ys.prediction[i] <- fun_inv(as.numeric(predict(glmnet.Training,t(as.matrix(xs[i,])),type="response")),
                                ys.funvar1[i],ys.funvar2[i])
    glmnet.final=data.frame(as.matrix(coef(glmnet.Training,s=lambda.glmnet.Training)))
    names(glmnet.final)='beta'
    glmnet.final$var=rownames(glmnet.final)
    glmnet.final=subset(glmnet.final,select=c(var,beta))
    glmnet.final.list[[i]] <- subset(glmnet.final,abs(beta)>0)
    cpg.count[i] <- nrow(subset(glmnet.final,abs(beta)>0))
  }
  save(glmnet.final.list, file=out.rdata)
  print("Debug checkpoint 2")
  ys.output <- ys
  ys.output[,PREDVAR]=ys.prediction
  attr(ys.output[,PREDVAR],'dimnames')<-NULL
  ys.output[,RESVAR]=residuals(lm(ys.prediction~ys.outcome))
  attr(ys.output[,RESVAR],'dimnames')<-NULL
  ys.output[,"cpg_count_used"]=cpg.count
  print("Debug checkpoint 3")
  if (!is.null(output.csv)) {
    write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
  }
  if (!is.null(out.png)) {
    TITLE_str <- paste0(out.png.title,'\n')
    if (show.legend) {
      width=9
      height=8
    } else {
      width=7
      height=8
    }
    saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str,width=width,height=height,show.legend=show.legend,oma.right=oma.right,NUMVAR=NUMVAR)
  }

  return(ys.output)
}

### Estimates outcome on a new data-info pair, using LOSO glmnet models
### Saves the info with predicted values and the correlation plot
### Returns: a data.frame, ys.output, containing ys with 4 more columns, being
###          PREDVAR, RESVAR, "cpg_count_used", and "log_lambda_hat"
#' @export
saveLOSOEstimation <- function(xs,ys,OUTVAR,SPECVAR,out.rdata,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA=0.5,NFOLD=10,fun_trans=fun_identity,fun_inv=fun_identity,fun_VAR1=NULL,fun_VAR2=NULL,COLVAR=NULL,show.legend=T,oma.right=9,NUMVAR=NULL,cpg_names.list=NULL,xs.add_train=NULL,ys.add_train=NULL,loglambda.seq=NULL,parallel=FALSE) {
  #OUTVAR=String, name of dependent variable in ys
  #SPECVAR=String, name of grouping variable (FACTOR) in ys
  #out.rdata <- list containing table of coefficients for each LOSO model
  #output.csv <- ys.output
  #out.png <- image of validation plot
  #out.png.title=String, title used in validation plot
  #PREDVAR=String, desired name of the predicted value for the dependent variable
  #RESVAR=String, desired name of the predicted value's residual on the outcome
  #ALPHA=1 for lasso, ALPHA=0 for ridge
  #fun_trans=Function, transformation applied to OUTVAR before fitting the model
  #fun_inv=Function, topological inverse of fun_trans
  #fun_VAR1=String, name of variable to be used as the 1st parameter of fun_trans and fun_inv (default is NULL)
  #fun_VAR2=String, name of variable to be used as the 2nd parameter of fun_trans and fun_inv (default is NULL)
  #COLVAR=String, name of variable (FACTOR) by which to color the points (default is NULL)
  #NUMVAR=String, name of variable (FACTOR) by which to number the points (default is NULL)
  #cpg_names.list=List, string vectors containing CpG names for subsetting at each iteration
  #               NOTE: MUST HAVE length(cpg_names.list)==levels(ys[,SPECVAR])
  #xs.add_train/ys.add_train=Additional data to be included in the training set at every iteration
  rows.keep <- which(!is.na(fun_trans(ys[,OUTVAR],ys[,fun_VAR1],ys[,fun_VAR2])))
  if (length(rows.keep) < nrow(ys)) {
    print(paste(nrow(ys)-length(rows.keep),"rows were ignored"))
  }
  ys <- ys[rows.keep,]
  xs <- xs[rows.keep,]
  ys.outcome=ys[,OUTVAR]
  ys.funvar1=ys[,fun_VAR1] # 0 columns iff VAR=NULL
  ys.funvar2=ys[,fun_VAR2]
  if (is.null(fun_VAR1)) {ys.funvar1=NULL}
  if (is.null(fun_VAR2)) {ys.funvar2=NULL}
  ys.species=ys[,SPECVAR]
  if (!is.null(xs.add_train) & !is.null(ys.add_train)) {
    rows.keep <- which(!is.na(fun_trans(ys.add_train[,OUTVAR],ys.add_train[,fun_VAR1],ys.add_train[,fun_VAR2])))
    ys.add_train <- ys.add_train[rows.keep,]
    xs.add_train <- xs.add_train[rows.keep,]
    ys.add_train_outcome=ys.add_train[,OUTVAR]
    ys.add_train_funvar1=ys.add_train[,fun_VAR1] # 0 columns iff VAR=NULL
    ys.add_train_funvar2=ys.add_train[,fun_VAR2]
    if (is.null(fun_VAR1)) {ys.add_train_funvar1=NULL}
    if (is.null(fun_VAR2)) {ys.add_train_funvar2=NULL}
  }
  print("Debug checkpoint 1")

  if (parallel) {
    library(doMC)
    registerDoMC()
  }

  ### Building LOSO Clocks for each level of SPECVAR
  ys.prediction <- 0*ys.outcome
  glmnet.final.list <- vector("list",length(ys.prediction))
  cpg.count <- 0*ys.outcome
  log_lambda.hat <- 0*ys.outcome
  if (is.null(loglambda.seq)) {lambda=NULL} else {lambda=exp(loglambda.seq)}
  for (spec.num in 1:length(levels(ys.species))) {
    print(paste0(round((spec.num-1)/length(levels(ys.species))*100),"%"))
    indcs <- which(as.numeric(ys.species) == spec.num)
    if (is.null(cpg_names.list)) {cpg_names.num=colnames(xs)} else {cpg_names.num=cpg_names.list[[spec.num]]}
    if (length(indcs) > 0) {
      xs.trn <- xs[-indcs,]
      ys.trn_outcome <- ys.outcome[-indcs]
      ys.trn_funvar1 <- ys.funvar1[-indcs]
      ys.trn_funvar2 <- ys.funvar2[-indcs]
      if (!is.null(xs.add_train) & !is.null(ys.add_train)) {
        xs.trn <- rbind(xs.trn, xs.add_train)
        ys.trn_outcome <- c(ys.trn_outcome, ys.add_train_outcome)
        ys.trn_funvar1 <- c(ys.trn_funvar1, ys.add_train_funvar1)
        ys.trn_funvar2 <- c(ys.trn_funvar2, ys.add_train_funvar2)
      }
      glmnet.Training.CV = cv.glmnet(xs.trn[,cpg_names.num],fun_trans(ys.trn_outcome,ys.trn_funvar1,ys.trn_funvar2),
                                     nfolds=NFOLD,alpha=ALPHA,family="gaussian",type.measure="mse",lambda=lambda,parallel=parallel)
      #print("Debug checkpoint 1.1")
      lambda.glmnet.Training = glmnet.Training.CV$lambda.1se
      #print(log(lambda.glmnet.Training))
      glmnet.Training <- glmnet(xs.trn[,cpg_names.num],fun_trans(ys.trn_outcome,ys.trn_funvar1,ys.trn_funvar2),
                                alpha=ALPHA,family="gaussian",lambda=lambda.glmnet.Training)
      #print("Debug checkpoint 1.2")
      xs.test <- xs[indcs,cpg_names.num]
      if (length(indcs)==1) {
        xs.test <- t(as.matrix(xs[indcs,cpg_names.num]))
      }
      ys.prediction[indcs] <- fun_inv(as.numeric(predict(glmnet.Training,xs.test,type="response")),
                                      ys.funvar1[indcs],ys.funvar2[indcs])
      glmnet.final=data.frame(as.matrix(coef(glmnet.Training,s=lambda.glmnet.Training)))
      names(glmnet.final)='beta'
      glmnet.final$var=rownames(glmnet.final)
      glmnet.final=subset(glmnet.final,select=c(var,beta))
      for (i in indcs) {
        glmnet.final.list[[i]] <- subset(glmnet.final,abs(beta)>0)
      }
      cpg.count[indcs] <- nrow(subset(glmnet.final,abs(beta)>0))
      log_lambda.hat[indcs] <- log(lambda.glmnet.Training)
    }
  }
  print("100%")
  save(glmnet.final.list, file=out.rdata)
  print("Debug checkpoint 2")
  ys.output <- ys
  ys.output[,PREDVAR]=ys.prediction
  attr(ys.output[,PREDVAR],'dimnames')<-NULL
  ys.output[,RESVAR]=residuals(lm(ys.prediction~ys.outcome))
  attr(ys.output[,RESVAR],'dimnames')<-NULL
  ys.output[,"cpg_count_used"]=cpg.count
  ys.output[,"log_lambda_hat"]=log_lambda.hat
  print("Debug checkpoint 3")
  if (!is.null(output.csv)) {
    write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
  }
  if (!is.null(out.png)) {
    TITLE_str <- paste0(out.png.title,'\n')
    if (show.legend) {
      width=9
      height=8
    } else {
      width=7
      height=8
    }
    saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str,width=width,height=height,show.legend=show.legend,oma.right=oma.right,NUMVAR=NUMVAR)
  }

  return(ys.output)
}

### Estimates outcome on a new data-info pair, using a set of LOO glmnet models, stratified by a variable
### Saves the info with predicted values and the correlation plot
### Returns: a data.frame, ys.output, containing ys with 3 more columns, being
###          PREDVAR, RESVAR, and "cpg_count_used"
#' @export
saveLOOEstimationStrat <- function(xs,ys,OUTVAR,STRATVAR,out.rdata.PREFIX,output.csv,out.png,out.png.title,PREDVAR,RESVAR,ALPHA=0.5,NFOLD=10,fun_trans=fun_identity,fun_inv=fun_identity,fun_VAR1=NULL,fun_VAR2=NULL,COLVAR=NULL,show.legend=T,oma.right=9,NUMVAR=NULL,xs.add_train=NULL,ys.add_train=NULL,parallel=FALSE) {
  #OUTVAR=String, name of dependent variable in ys
  #STRATVAR=String, name of variable (FACTOR) in ys by which to stratify the elastic nets
  #         NOTE: CANNOT be NULL
  #out.rdata.PREFIX <- set of lists containing table of coefficients for each LOO model
  #output.csv <- ys.output
  #out.png <- image of validation plot
  #out.png.title=String, title used in validation plot
  #PREDVAR=String, desired name of the predicted value for the dependent variable
  #RESVAR=String, desired name of the predicted value's residual on the outcome
  #ALPHA=1 for lasso, ALPHA=0 for ridge
  #fun_trans=Function, transformation applied to OUTVAR before fitting the model
  #fun_inv=Function, topological inverse of fun_trans
  #fun_VAR1=String, name of variable to be used as the 1st parameter of fun_trans and fun_inv (default is NULL)
  #fun_VAR2=String, name of variable to be used as the 2nd parameter of fun_trans and fun_inv (default is NULL)
  #COLVAR=String, name of variable (FACTOR) by which to color the points (default is NULL)
  #NUMVAR=String, name of variable (FACTOR) by which to number the points (default is NULL)
  #xs.add_train/ys.add_train=Additional data to be included in the training set at every iteration
  if (is.null(STRATVAR)) {
    stop('STRATVAR Cannot be NULL')
  }

  rows.keep <- which(!is.na(fun_trans(ys[,OUTVAR],ys[,fun_VAR1],ys[,fun_VAR2])))
  if (length(rows.keep) < nrow(ys)) {
    print(paste(nrow(ys)-length(rows.keep),"rows were ignored"))
  }
  ys <- ys[rows.keep,]
  xs <- xs[rows.keep,]
  ys.outcome=ys[,OUTVAR]
  ys.funvar1=ys[,fun_VAR1] # 0 columns iff VAR=NULL
  ys.funvar2=ys[,fun_VAR2]
  if (is.null(fun_VAR1)) {ys.funvar1=NULL}
  if (is.null(fun_VAR2)) {ys.funvar2=NULL}
  ys.strata=ys[,STRATVAR]
  if (!is.null(xs.add_train) & !is.null(ys.add_train)) {
    rows.keep <- which(!is.na(fun_trans(ys.add_train[,OUTVAR],ys.add_train[,fun_VAR1],ys.add_train[,fun_VAR2])))
    ys.add_train <- ys.add_train[rows.keep,]
    xs.add_train <- xs.add_train[rows.keep,]
    ys.add_train_outcome=ys.add_train[,OUTVAR]
    ys.add_train_funvar1=ys.add_train[,fun_VAR1] # 0 columns iff VAR=NULL
    ys.add_train_funvar2=ys.add_train[,fun_VAR2]
    if (is.null(fun_VAR1)) {ys.add_train_funvar1=NULL}
    if (is.null(fun_VAR2)) {ys.add_train_funvar2=NULL}
  }
  print("Debug checkpoint 1")

  if (parallel) {
    library(doMC)
    registerDoMC()
  }

  ### Building set of LOO Clocks for each level of STRATVAR
  ys.output <- data.frame()
  for (strat.num in 1:length(levels(ys.strata))) {
    print(paste0(round((strat.num-1)/length(levels(ys.strata))*100),"%"))
    indcs <- which(as.numeric(ys.strata) == strat.num)
    strat.name <- levels(ys.strata)[strat.num]
    stratName <- to_upper_camel_case(strat.name)
    # TODO: Refine If Statement #
    if (length(unique(ys.outcome[indcs])) <= 2) {
      print(paste(stratName,"could not be fit: Outcome values are equal"))
    } else {
      out.rdata=paste0(out.rdata.PREFIX,stratName,'.RData')
      ys.output.num <- saveLOOEstimation(xs[indcs,],ys[indcs,],OUTVAR,out.rdata,output.csv=NULL,out.png=NULL,out.png.title=NULL,PREDVAR,RESVAR,ALPHA,NFOLD,fun_trans=fun_trans,fun_inv=fun_inv,fun_VAR1=fun_VAR1,fun_VAR2=fun_VAR2,parallel=parallel)
      ys.output <- rbind(ys.output,ys.output.num)
    }
  }
  print("100%")
  ys.output[,STRATVAR] <- factor(ys.output[,STRATVAR]) # Correction, in case strata were skipped
  if (!is.null(output.csv)) {
    write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
  }
  if (!is.null(out.png)) {
    TITLE_str <- paste0(out.png.title,'\n')
    if (show.legend) {
      width=9
      height=8
    } else {
      width=7
      height=8
    }
    saveValidationPlot(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str,width=width,height=height,show.legend=show.legend,oma.right=oma.right,NUMVAR=NUMVAR)
  }

  return(ys.output)
}

### Computes the correlation of each CpG site with OUTVAR
### Saves a table with the correlations, z scores, p-values, q-values, and more
### Returns: a data.frame, ewas.table
#' @export
saveEWAS <- function(xs,ys,OUTVAR,ewas.csv,fun_trans=fun_identity,fun_VAR1=NULL,fun_VAR2=NULL) {
  #OUTVAR=String, name of dependent variable in ys
  #ewas.csv <- ewas.table
  #fun_trans=Function, transformation applied to OUTVAR before computing the correlations
  #fun_VAR1=String, name of variable to be used as the 1st parameter of fun_trans (default is NULL)
  #fun_VAR2=String, name of variable to be used as the 2nd parameter of fun_trans (default is NULL)
  rows.keep <- which(!is.na(fun_trans(ys[,OUTVAR],ys[,fun_VAR1],ys[,fun_VAR2])))
  if (length(rows.keep) < nrow(ys)) {
    print(paste(nrow(ys)-length(rows.keep),"rows were ignored"))
  }
  ys <- ys[rows.keep,]
  xs <- xs[rows.keep,]
  ys.outcome=ys[,OUTVAR]
  ys.funvar1=ys[,fun_VAR1] # 0 columns iff VAR=NULL
  ys.funvar2=ys[,fun_VAR2]
  if (is.null(fun_VAR1)) {ys.funvar1=NULL}
  if (is.null(fun_VAR2)) {ys.funvar2=NULL}

  ewas.table <- standardScreeningNumericTrait(xs,fun_trans(ys.outcome,ys.funvar1,ys.funvar2),
                                              areaUnderROC=F)
  if (!is.null(ewas.csv)) {
    write.table(ewas.table,ewas.csv,sep=',',row.names=F,quote=F)
  }

  return(ewas.table)
}

### Computes the correlation of each CpG site with OUTVAR, stratified by STRATVAR
### Saves a set of tables, each with the correlations, z scores, p-values, q-values, and more
### Saves an aggregate table with strata z scores, and meta-analysis z scores and p-values
### Saves a pairwise plot of the strata z scores with pairwise correlations
### Returns: a data.frame, ewas.meta.table
#' @export
saveEWASStrat <- function(xs,ys,OUTVAR,STRATVAR,ewas.csv.PREFIX,ewas.png,fun_trans=fun_identity,fun_VAR1=NULL,fun_VAR2=NULL,ewas.png.size=9) {
  #OUTVAR=String, name of dependent variable in ys
  #STRATVAR=String, name of variable (FACTOR) in ys by which to stratify the analysis
  #         NOTE: CANNOT be NULL
  #ewas.csv.PREFIX <- ewas.table.list AND ewas.meta.table
  #ewas.png <- pair(ewas.meta.table) with correlations
  #fun_trans=Function, transformation applied to OUTVAR before computing the correlations
  #fun_VAR1=String, name of variable to be used as the 1st parameter of fun_trans (default is NULL)
  #fun_VAR2=String, name of variable to be used as the 2nd parameter of fun_trans (default is NULL)
  #ewas.png.size=Numeric, width and height (in inches, at 300psi) of the image ewas.png (default is 10)
  rows.keep <- which(!is.na(fun_trans(ys[,OUTVAR],ys[,fun_VAR1],ys[,fun_VAR2])))
  if (length(rows.keep) < nrow(ys)) {
    print(paste(nrow(ys)-length(rows.keep),"rows were ignored"))
  }
  ys <- ys[rows.keep,]
  xs <- xs[rows.keep,]
  ys.outcome=ys[,OUTVAR]
  ys.funvar1=ys[,fun_VAR1] # 0 columns iff VAR=NULL
  ys.funvar2=ys[,fun_VAR2]
  if (is.null(fun_VAR1)) {ys.funvar1=NULL}
  if (is.null(fun_VAR2)) {ys.funvar2=NULL}
  ys.strata=ys[,STRATVAR]

  ewas.table.list <- vector("list",length(levels(ys.strata)))
  ewas.Z.list <- vector("list",length(levels(ys.strata)))
  strat.names <- vector("character",0)
  for (strat.num in 1:length(levels(ys.strata))) {
    indcs <- which(as.numeric(ys.strata) == strat.num)
    if (!is.na(var(ys.outcome[indcs], na.rm=T)) && var(ys.outcome[indcs], na.rm=T) > 0) {
      ewas.table.list[[strat.num]] <- standardScreeningNumericTrait(xs[indcs,],fun_trans(ys.outcome[indcs],ys.funvar1[indcs],ys.funvar2[indcs]),
                                                                    areaUnderROC=F)
      if (!is.null(ewas.csv.PREFIX)) {
        ewas.csv <- paste0(ewas.csv.PREFIX,'_',levels(ys.strata)[strat.num],'.csv')
        write.table(ewas.table.list[[strat.num]],ewas.csv,sep=',',row.names=F,quote=F)
      }
      ewas.Z.list[[strat.num]] <- standardScreeningNumericTrait(xs[indcs,],fun_trans(ys.outcome[indcs],ys.funvar1[indcs],ys.funvar2[indcs]),
                                                                areaUnderROC=F)$Z
      strat.names <- c(strat.names, paste0("Z.",levels(ys.strata)[strat.num]))
    } else {
      ewas.table.list[[strat.num]] <- NULL
      ewas.Z.list[[strat.num]] <- NULL
    }
  }
  ewas.meta.table <- as.data.frame(do.call(cbind,ewas.Z.list))
  colnames(ewas.meta.table) <- strat.names
  ewas.meta.table$metaZ.Age <- rowSums(ewas.meta.table)/sqrt(ncol(ewas.meta.table))
  ewas.meta.table$metaP.Age <- 2*(1-pnorm(abs(ewas.meta.table$metaZ.Age)))
  ewas.meta.table <- cbind(ID=colnames(xs), ewas.meta.table)
  if (!is.null(ewas.csv.PREFIX)) {
    ewas.meta.csv <- paste0(ewas.csv.PREFIX,'.csv')
    write.table(ewas.meta.table,ewas.meta.csv,sep=',',row.names=F,quote=F)
  }

  if (!is.null(ewas.png) && ncol(ewas.meta.table) > 4) {
    png(ewas.png,width=ewas.png.size,height=ewas.png.size,units='in',res=300)
    pairs(ewas.meta.table[,2:(ncol(ewas.meta.table)-2)],
          lower.panel=panel.cor,upper.panel=points,pch=16,
          cex.cor=1.8,cex.labels=1.4)
    dev.off()
  }

  return(ewas.meta.table)
}

### Maps each CpG in CPGVAR to its nearest gene in Human and Sample Species
### Augments df with the gene annotation, and Saves the output
### Returns: a data.frame, df.geneann, containing df with 6 more columns, being
###          "ensemblHuman", "ensemblSampSpecies", "geneSymbol",
###          "chrSampSpecies", "startSampSpecies", and "endSampSpecies"
#' @export
saveGeneAnnotation <- function(df,CPGVAR,sln,geneann.csv) {
  #CPGVAR=String, name of variable in df that contains CpG ID
  #SLN=String, name of variable (FACTOR) in df that contains species latin name
  #            OR species latin name (length 1 character vector)
  #geneann.csv <- df.geneann
  df.cpg=df[,CPGVAR]
  if (length(sln)==1) {df.sln=factor(rep(sln,nrow(df)))} else {df.sln=df[,SLN]}

  df.geneann <- df
  ### TODO: Make an if statement to check when df.cpg[i] does not begin with "cg"

  if (!is.null(geneann.csv)) {
    write.table(df.geneann,geneann.csv,sep=',',row.names=F,quote=F)
  }

  return(df.geneann)
}

# saveFunctionalAnnotation <- function(df,CPGVAR,funann.csv) {
#   ### TODO
#   return(NULL)
# }

### Generates random forest by training on a separate dataset
### Saves the info with predicted values and performance plot
### Returns: a data.frame, ys.output, containing ys.train+ys.test with 2 more columns, being
###          "train" (logical) and PREDVAR
#' @export
saveRFTrained <- function(xs.train,ys.train,xs.test,ys.test,OUTVAR,output.csv,out.png,out.png.title,PREDVAR) {
  #OUTVAR=String, name of dependent variable in ys.train and ys.test
  #output.csv <- ys.output
  #out.png <- image of diagnostic plot
  #out.png.title=String, title used in diagnostic plot
  #PREDVAR=String, desired name of the predicted value for the dependent variable
  ys.train_outcome=ys.train[,OUTVAR]
  ys.test_outcome=ys.test[,OUTVAR]

  RF.Training = randomForest(x=xs.train,y=ys.train_outcome,xtest=xs.test)
  print("Debug checkpoint 1")

  train=data.frame(ys.train,train=1)
  ys.train_prediction=ys.train_outcome
  train[,PREDVAR]=ys.train_prediction
  attr(train[,PREDVAR],'dimnames')<-NULL
  if (!is.null(ys.test)) {
    test=data.frame(ys.test,train=0)
    ys.test_prediction=RF.Training$test$predicted
    test[,PREDVAR]=ys.test_prediction
    attr(test[,PREDVAR],'dimnames')<-NULL
    ys.output=rbind_merge(train,test)
  }
  if (is.null(ys.test)) {
    ys.output <- train
  }
  if (!is.null(output.csv)) {
    write.table(ys.output,output.csv,sep=',',row.names=F,quote=F)
  }
  print("Debug checkpoint 2")

  png(out.png,width=7,height=7,units='in',res=300)
  plot(RF.Training,cex.axis=1.5,cex.lab=1.5)
  dev.off()
  return(ys.output)
}

### Generates Cox PH regression model on a given set of covariates
### Saves the GEE summary table of coefficients and survival curve plot
### Returns nothing
#' @export
saveSurvAnalysis <- function(ys.output,ENDVAR,STATUSVAR,COVAR,coxphtable.tsv,coxphtable.png,coxph.png,coxphcombined.png,coxph.png.title,STRATVAR=NULL,CLUSTVAR=NULL,coxph.png.xlab=NULL) {
  #ys.output <- df used in model, containing ENDVAR,STATUSVAR,COVAR,STRATVAR,CLUSTVAR
  #ENDVAR=String, name of endpoint variable in model formula
  #STATUSVAR=String, name of endpoint status variable (LOGICAL) in model formula
  #COVAR=String, name of variable(s) in model formula that are to be used as covariates
  #coxphtable.tsv <- GEE summary table
  #coxphtable.png <- image of GEE summary table
  #coxph.png <- image of survival curve
  #coxphcombined.png <- image of survival curve and GEE summary table
  #coxph.png.title=String, title used in all plots
  #STRATVAR=String, name of variable(s) (FACTOR) in model formula by which to stratify the analysis (default is NULL)
  #CLUSTVAR=String, name of variable(s) (FACTOR) in model formula by which to cluster the data (default is NULL)
  if (is.null(COVAR)) {
    COVAR="1"
  }
  fm.str <- paste0("Surv(",ENDVAR,",",STATUSVAR,")~", paste(COVAR,collapse="+"))
  if (!is.null(STRATVAR)) {
    fm.str <- paste0(fm.str,"+strata(",paste(STRATVAR,collapse=","),")")
  }
  if (!is.null(CLUSTVAR)) {
    fm.str <- paste0(fm.str,"+cluster(",paste(CLUSTVAR,collapse=","),")")
  }

  coxph.object <- coxph(as.formula(fm.str), data=ys.output, na.action=na.exclude)
  #coxph.png.subtitle = paste0("Concordance = ",round(coxph.object$concordance["concordance"],3)," (",round(coxph.object$concordance["std"],3),")")
  coxph.table <- tabcoxph(coxph.object, factor.compression=2, decimals=3)
  if (!is.null(coxphtable.tsv)) {
    write.table(coxph.table,coxphtable.tsv,sep='\t',row.names=F,quote=F)
  }

  if (!is.null(coxph.png.xlab)) {xlab=coxph.png.xlab} else {xlab=ENDVAR}
  p <- autoplot(survfit(coxph.object), xlab=xlab, ylab="Survival", main=coxph.png.title) +
    theme_bw(base_size=13) + theme(plot.title=element_text(size=10, hjust=0.5))
  tg <- tableGrob(coxph.table, rows=NULL,
                  theme=ttheme_default(base_size=13, padding=unit(c(0.25,0.15),"in")))
  tg.title <- textGrob(coxph.png.title, gp=gpar(fontsize=7))
  tg2 <- gtable::gtable_add_rows(tg, heights=grobHeight(tg.title)+unit(0.15,"in"), pos=0)
  tg2 <- gtable::gtable_add_grob(tg2, tg.title, 1, 1, 1, ncol(tg))

  ## PLOTS grDevices::dev.size() is the default size when length(grDevices::dev.list()) != 0
  ## PLOTS c(7,7) is the default size when length(grDevices::dev.list()) == 0
  ggsave(coxph.png, p, width=7, height=7)
  ggsave(coxphtable.png, grid.draw(tg2),
         width=sum(convertX(tg2$widths,"in",valueOnly=T))+0.15,
         height=sum(convertX(tg2$heights,"in",valueOnly=T))+0.15)
  dev.off() # Use dev.off() after using convertX()
  ggsave(coxphcombined.png, grid.arrange(p,tg,heights=unit.c(unit(1,"null"),sum(tg$heights,unit(0.15,"in")))),
         width=7, height=7)
}

### Plots outcome versus prediction, colored by a factor variable
### Saves the plot as an image
### Returns nothing
#' @export
saveValidationPlot <- function(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str,width=7,height=8,show.legend=T,oma.right=9,NUMVAR=NULL,ys.colors=NULL,ys.numbers=NULL) {
  #ys.output <- df containing OUTVAR,PREDVAR,COLVAR
  #OUTVAR=String, name of dependent variable in ys.output
  #PREDVAR=String, name of the predicted value for the dependent variable in ys.output
  #COLVAR=String, name of variable (FACTOR) in ys.output by which to color the points
  #       NOTE: CAN be NULL
  #out.png <- image of validation plot
  #TITLE_str=String, title used in validation plot
  #NUMVAR=String, name of variable (FACTOR) in ys.output by which to number the points
  #       NOTE: CAN be NULL
  #ys.colors=String vector, colors by which to color the points, according to levels(COLVAR)
  #          NOTE: CAN be NULL (default is NULL, uses Dark2() or rainbow())
  #ys.numbers=String vector, numbers by which to number the points, according to levels(NUMVAR)
  #           NOTE: CAN be NULL (default is NULL, uses as.numeric(factor))
  if (is.null(out.png)) {
    return()
  }

  ys.outcome <- ys.output[,OUTVAR]
  ys.prediction <- ys.output[,PREDVAR]
  if (!is.null(COLVAR)) {
    ys.colfactor <- ys.output[,COLVAR]
    # ys.colors contains the palette of distinct colors
    # ys.colors_vec contains the color assignment for every row in the data set
    if (is.null(ys.colors)) {
      if (length(levels(ys.colfactor)) <= 2) {
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
      } else if (length(levels(ys.colfactor)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor)))
      }
    }
  } else {
    ys.colfactor <- rep(1, nrow(ys.output))
    ys.colors <- 'black'
  }
  if (!is.null(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.null(ys.numbers)) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]

  png(out.png,width=width,height=height,units='in',res=300)
  par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
  lim <- get.square_limits(ys.prediction, ys.outcome)
  l_lim=lim[1]
  u_lim=lim[2]
  MAE <- median(abs(ys.prediction-ys.outcome), na.rm=T)
  MAE_str <- paste0("MAE=",signif(MAE,3))
  N <- length(which(!is.na(ys.outcome) & !is.na(ys.prediction)))
  if (!is.na(var(ys.outcome, na.rm=T)) && var(ys.outcome, na.rm=T) > 0) {
    COR <- cor(ys.prediction,ys.outcome,
               use='pairwise.complete.obs')
    COR_str <- paste0("cor=",signif(COR,2))
    MAE_str <- paste0(MAE_str,', ',COR_str)
  }
  plot(y=ys.prediction,x=ys.outcome,
       type=type,main=paste0('(N=',N,')'),
       ylab=PREDVAR,xlab=OUTVAR,pch=ys.numbers_vec,
       cex.main=1.5,cex.lab=1.5,cex.axis=1.5,
       col=ys.colors_vec,xlim=lim,ylim=lim)
  if (type=='n') {
    text(y=ys.prediction,x=ys.outcome,cex=0.8,
         labels=ys.numbers_vec,col=ys.colors_vec)
  }
  if (!is.na(var(ys.outcome, na.rm=T)) && var(ys.outcome, na.rm=T) > 0 && !is.na(var(ys.prediction, na.rm=T))) {
    abline(lm(ys.prediction~ys.outcome))
  }
  abline(0,1,lty="dashed")
  title(MAE_str,outer=F,line=0.4,cex.main=1.5)
  if (show.legend && !is.null(COLVAR)) {
    legend(y=u_lim,x=u_lim+0.05*(u_lim-l_lim),
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.1,pt.cex=2,xpd=NA)
  }
  title(TITLE_str,outer=T,line=0,cex.main=1.0)
  dev.off()
}

### Plots outcome versus prediction, colored by a factor variable
###       Also: makes x-axis wider, and limits legend to top 9 factor levels
### Saves the plot as an image
### Returns nothing
#' @export
saveValidationPlotManyLevels <- function(ys.output,OUTVAR,PREDVAR,COLVAR,out.png,TITLE_str,width=7,height=8,show.legend=F,xlim=NULL,NUMVAR=NULL,ys.colors=NULL,ys.numbers=NULL) {
  #ys.output <- df containing OUTVAR,PREDVAR,COLVAR
  #OUTVAR=String, name of dependent variable in ys.output
  #PREDVAR=String, name of the predicted value for the dependent variable in ys.output
  #COLVAR=String, name of variable (FACTOR) in ys.output by which to color the points
  #       NOTE: CANNOT be NULL
  #out.png <- image of validation plot
  #TITLE_str=String, title used in validation plot
  #xlim=2-element Vector, desired limits of x-axis (default is NULL)
  #NUMVAR=String, name of variable (FACTOR) in ys.output by which to number the points
  #       NOTE: CAN be NULL
  #ys.colors=String vector, colors by which to color the points, according to levels(COLVAR)
  #          NOTE: CAN be NULL (default is NULL, uses Dark2() or rainbow())
  #ys.numbers=String vector, numbers by which to number the points, according to levels(NUMVAR)
  #           NOTE: CAN be NULL (default is NULL, uses as.numeric(factor))
  if (is.null(out.png)) {
    return()
  }

  ys.outcome <- ys.output[,OUTVAR]
  ys.prediction <- ys.output[,PREDVAR]
  ys.colfactor <- ys.output[,COLVAR]
  # ys.colors contains the palette of distinct colors
  # ys.colors_vec contains the color assignment for every row in the data set
  if (is.null(ys.colors)) {
    if (length(levels(ys.colfactor)) <= 2) {
      ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
    } else if (length(levels(ys.colfactor)) <= 8) {
      ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
    } else {
      ys.colors <- rainbow(length(levels(ys.colfactor)))
    }
  }
  if (!is.null(NUMVAR)) {
    type='n'
    ys.numfactor <- ys.output[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (is.null(ys.numbers)) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor <- rep(1, nrow(ys.output))
    ys.numbers <- 16
  }
  ys.colors_vec <- ys.colors[ys.colfactor]
  ys.numbers_vec <- ys.numbers[ys.numfactor]
  idx.colfactor_levels.sorted <- match(names(sort(table(ys.colfactor),decreasing=T)), levels(ys.colfactor))

  png(out.png,width=width,height=height,units='in',res=300)
  par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2.5,0)) #oma_top extended 0.5 because MAE_str is 2 lines
  MAE <- median(abs(ys.prediction-ys.outcome), na.rm=T)
  MAE_str <- paste0("MAE=",signif(MAE,3))
  N <- length(which(!is.na(ys.outcome) & !is.na(ys.prediction)))
  if (!is.na(var(ys.outcome, na.rm=T)) && var(ys.outcome, na.rm=T) > 0) {
    COR <- cor(ys.prediction,ys.outcome,
               use='pairwise.complete.obs')
    COR_str <- paste0("cor=",signif(COR,2))
    cor.within_table <- as.data.frame(ys.output %>% transmute(val=ys.outcome,pred=ys.prediction,group=ys.colfactor) %>%
                                        group_by(group) %>% transmute(cor.within=cor(val,pred), grp.size=length(val)))
    COR.within_vec <- unique(dplyr::filter(cor.within_table, grp.size >= 10))$cor.within
    MEDICOR <- median(COR.within_vec, na.rm=T)
    MEDICOR_str <- paste0("Median within-group cor=",signif(MEDICOR,3))
    MEANCOR <- mean(COR.within_vec, na.rm=T)
    MEANCOR_str <- paste0("Mean within-group cor=",signif(MEANCOR,3))
    MAE_str <- paste0(MAE_str,', ',COR_str,', ',MEDICOR_str,'\n',MEANCOR_str)
  }
  plot(y=ys.prediction,x=ys.outcome,
       type=type,main=NULL,
       ylab=PREDVAR,xlab=OUTVAR,pch=ys.numbers_vec,
       cex.main=1.5,cex.lab=1.5,cex.axis=1.5,
       col=ys.colors_vec,xlim=xlim)
  title(paste0('(N=',N,')'),line=3.5,cex.main=1.5) #title 'line' changed from default (centered) because MAE_str is 2 lines
  if (type=='n') {
    text(y=ys.prediction,x=ys.outcome,cex=0.8,
         labels=ys.numbers_vec,col=ys.colors_vec)
  }
  if (!is.na(var(ys.outcome, na.rm=T)) && var(ys.outcome, na.rm=T) > 0 && !is.na(var(ys.prediction, na.rm=T))) {
    abline(lm(ys.prediction~ys.outcome))
  }
  abline(0,1,lty="dashed")
  title(MAE_str,outer=F,line=0.4,cex.main=1.5)
  if (show.legend) {
    legend("bottomright",
           legend=levels(ys.colfactor)[idx.colfactor_levels.sorted[1:9]],
           col=ys.colors[idx.colfactor_levels.sorted[1:9]],
           pch=16,cex=1.1,pt.cex=2,xpd=NA)
  }
  title(TITLE_str,outer=T,line=0.5,cex.main=1.0) #TITLE_str moved up 0.5 because MAE_str is 2 lines
  dev.off()
}

### Plots original outcome versus prediction, colored by a factor variable
### Plots individual sub-groups as well, creating a panel plot
### Saves the plot as an image
### Returns nothing
#' @export
saveValidationPanelPlot <- function(ys.output,OUTVAR,PREDVAR,PANELVAR,out.png,TITLE_str,mfrow,width=13,height=14,COLVAR=NULL,show.original=T,y.axis.labs=PREDVAR,x.axis.labs=OUTVAR,panel.labs=LETTERS,show.legend=F,oma.right=0,show.panel_main=T,panel.mains=NULL,equal.axes=T,ys.colors=NULL) {
  #ys.output <- df containing OUTVAR,PREDVAR,PANELVAR
  #OUTVAR=String, name of dependent variable in ys.output
  #PREDVAR=String, name of the predicted value for the dependent variable in ys.output
  #PANELVAR=String, name of variable (FACTOR) in ys.output by which to partition into panels
  #         NOTE: CANNOT be NULL
  #out.png <- image of validation panel plot
  #TITLE_str=String, title used in validation panel plot
  #mfrow=2-element Vector, desired rxc dimensions of the panel plot
  #COLVAR=String, name of variable (FACTOR) in ys.output by which to color points
  #       NOTE: CAN be NULL (default is NULL, uses PANELVAR)
  #show.original=Logical, whether to show the original full-data plot in the first panel (default is TRUE)
  #y.axis.labs=String Vector, labels to put on y-axis of each panel (default is PREDVAR)
  #x.axis.labs=String Vector, labels to put on x-axis of each panel (default is OUTVAR)
  #panel.labs=String Vector, labels to put in upper-left of each panel (default is LETTERS)
  #show.panel_main=Logical, whether to print PANELVAR at the top of each panel (default is TRUE)
  #panel.mains=String Vector, titles to put on each panel (default is levels(PANELVAR))
  #            NOTE: Must be length 1 OR length(levels(PANELVAR))
  #equal.axes=Logical, whether to force all panels to have the same x- and y-axis limits (default is TRUE)
  #ys.colors=String vector, colors by which to color the points, according to levels(COLVAR)
  #          NOTE: CAN be NULL (default is NULL, uses Dark2() or rainbow())
  if (is.null(out.png)) {
    return()
  }
  if (is.null(COLVAR)) {
    COLVAR = PANELVAR
  }

  ys.outcome <- ys.output[,OUTVAR]
  ys.prediction <- ys.output[,PREDVAR]
  ## ASSERT(!is.null(PANELVAR)) ##
  if (is.null(PANELVAR)) {stop("PANELVAR cannot be NULL")}
  ys.panelfactor <- ys.output[,PANELVAR]
  if (mfrow[1]*mfrow[2] < length(levels(ys.panelfactor))) {
    stop("Error in mfrow: All levels of PANELVAR cannot fit within a single panel")
  }
  if (length(y.axis.labs)==1) {
    y.axis.labs <- rep(y.axis.labs, length(levels(ys.panelfactor))+1)
  }
  if (length(x.axis.labs)==1) {
    x.axis.labs <- rep(x.axis.labs, length(levels(ys.panelfactor))+1)
  }
  if (is.null(panel.mains)) {
    panel.mains = levels(ys.panelfactor)
  }
  if (length(panel.mains) != 1 & length(panel.mains) != length(levels(ys.panelfactor))) {
    print("Error in panel.mains: Must be length 1 or match levels of PANELVAR")
    return()
  }
  ys.colfactor <- ys.output[,COLVAR]
  # ys.colors contains the palette of distinct colors
  # ys.colors_vec contains the color assignment for every row in the data set
  if (is.null(ys.colors)) {
    if (length(levels(ys.colfactor)) <= 2) {
      ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor))]
    } else if (length(levels(ys.colfactor)) <= 8) {
      ys.colors <- brewer.pal(length(levels(ys.colfactor)), "Dark2")
    } else {
      ys.colors <- rainbow(length(levels(ys.colfactor)))
    }
  }
  ys.colors_vec <- ys.colors[ys.colfactor]

  png(out.png,width=width,height=height,units='in',res=300)
  par(mfrow=mfrow)
  par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,oma.right))
  lim <- get.square_limits(ys.prediction, ys.outcome)
  l_lim=lim[1]
  u_lim=lim[2]
  if (show.original) {
    MAE <- median(abs(ys.prediction-ys.outcome), na.rm=T)
    MAE_str <- paste0("MAE=",signif(MAE,3))
    PANEL_str <- 'All'
    N <- length(which(!is.na(ys.outcome) & !is.na(ys.prediction)))
    ylab=y.axis.labs[1]
    xlab=x.axis.labs[1]
    plab=panel.labs[1]
    if (!is.na(var(ys.outcome, na.rm=T)) && var(ys.outcome, na.rm=T) > 0) {
      COR <- cor(ys.prediction,ys.outcome,
                 use='pairwise.complete.obs')
      COR_str <- paste0("cor=",signif(COR,2))
      MAE_str <- paste0(MAE_str,', ',COR_str)
    }
    plot(y=ys.prediction,x=ys.outcome,
         main=paste0(PANEL_str,' (N=',N,')'),
         ylab=ylab,xlab=xlab,pch=16,
         cex.main=1.5,cex.lab=1.5,cex.axis=1.5,
         col=ys.colors_vec,xlim=lim,ylim=lim)
    if (!is.na(var(ys.outcome, na.rm=T)) && var(ys.outcome, na.rm=T) > 0 && !is.na(var(ys.prediction, na.rm=T))) {
      abline(lm(ys.prediction~ys.outcome))
    }
    abline(0,1,lty="dashed")
    title(MAE_str,outer=F,line=0.4,cex.main=1.5)
    mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  }
  for (i in 1:length(levels(ys.panelfactor))) {
    rows_i <- which(as.numeric(ys.panelfactor)==i)
    if (!equal.axes) {
      lim <- get.square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
      l_lim=lim[1]
      u_lim=lim[2]
      if (is.null(l_lim)) {
        lim <- get.square_limits(ys.prediction, ys.outcome)
        l_lim=lim[1]
        u_lim=lim[2]
      }
    }
    MAE <- median(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
    MAE_str <- paste0("MAE=",signif(MAE,3))
    if (show.panel_main) {PANEL_str <- panel.mains[i]}
    else if (!show.panel_main) {PANEL_str <- ""}
    N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
    if (show.original) {
      ylab=y.axis.labs[i+1]
      xlab=x.axis.labs[i+1]
      plab=panel.labs[i+1]
    } else {
      ylab=y.axis.labs[i]
      xlab=x.axis.labs[i]
      plab=panel.labs[i]
    }
    if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0) {
      COR <- cor(ys.prediction[rows_i],ys.outcome[rows_i],
                 use='pairwise.complete.obs')
      COR_str <- paste0("cor=",signif(COR,2))
      MAE_str <- paste0(MAE_str,', ',COR_str)
    }
    plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
         main=paste0(PANEL_str,' (N=',N,')'),
         ylab=ylab,xlab=xlab,pch=16,
         cex.main=1.5,cex.lab=1.5,cex.axis=1.5,
         col=ys.colors_vec[rows_i],xlim=lim,ylim=lim)
    if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0 && !is.na(var(ys.prediction[rows_i], na.rm=T))) {
      abline(lm(ys.prediction[rows_i]~ys.outcome[rows_i]))
    }
    abline(0,1,lty="dashed")
    title(MAE_str,outer=F,line=0.4,cex.main=1.5)
    mtext(plab,at=l_lim,adj=1,font=2,cex=1.4)
  }
  if (show.legend) {
    #creates legend to the right of the final panel
    #aligned to the bottom
    y_leg_row_height=1.7*(u_lim - l_lim)*mfrow[1]/height*0.166
    y_leg_height=y_leg_row_height*length(levels(ys.colfactor))+y_leg_row_height
    legend(y=l_lim+1.6*y_leg_height,x=u_lim+0.05*(u_lim-l_lim),
           legend=levels(ys.colfactor),
           col=ys.colors,
           pch=16,cex=1.6,pt.cex=3,xpd=NA)
  }
  title(TITLE_str,outer=T,line=-1,cex.main=1.5)
  dev.off()
}

### Plots outcome versus prediction for each species group,
###       colored by a factor variable, creating a series of panel plots
### Saves the plots as a set of images
### Returns nothing
#' @export
saveValidationPanelPlotSeriesForUniversal <- function(ys.output_original,OUTVAR,PREDVAR,COLVAR=NULL,out.png.PREFIX,TITLE_str.PREFIX,mfrow,NUMVAR=NULL,PANELVAR="SpeciesLatinName",sort.by_size=T,ys.colors=NULL,ys.numbers=NULL) {
  #ys.output_original <- original df containing OUTVAR,PREDVAR,COLVAR
  #OUTVAR=String, name of dependent variable in ys.output
  #PREDVAR=String, name of the predicted value for the dependent variable in ys.output
  #COLVAR=String, name of variable (FACTOR) in ys.output by which to color the points
  #       NOTE: CAN be NULL (default is NULL)
  #out.png.PREFIX <- PREFIX for image of validation panel plots
  #TITLE_str.PREFIX=String, PREFIX for title used in validation panel plots
  #mfrow=2-element Vector, desired rxc dimensions of the panel plots
  #NUMVAR=String, name of variable (FACTOR) in ys.output by which to number the points
  #       NOTE: CAN be NULL
  #PANELVAR=String, name of variable (FACTOR) in ys.output by which to partition into panels
  #         NOTE: CANNOT be NULL (default is "SpeciesLatinName")
  #sort.by_size=Logical, whether to sort panels by number of samples in each panel (default is TRUE)
  #ys.colors=String vector, colors by which to color the points, according to levels(COLVAR)
  #          NOTE: CAN be NULL (default is NULL, uses Dark2() or rainbow())
  #ys.numbers=String vector, numbers by which to number the points, according to levels(NUMVAR)
  #           NOTE: CAN be NULL (default is NULL, uses as.numeric(factor))
  if (is.null(out.png.PREFIX)) {
    return()
  }
  panel.size = mfrow[1]*mfrow[2]
  if (!is.null(COLVAR)) {
    ys.colfactor_original <- ys.output_original[,COLVAR]
    # ys.colors contains the palette of distinct colors
    if (is.null(ys.colors)) {
      if (length(levels(ys.colfactor_original)) <= 2) {
        ys.colors <- brewer.pal(3, "Dark2")[1:length(levels(ys.colfactor_original))]
      } else if (length(levels(ys.colfactor_original)) <= 8) {
        ys.colors <- brewer.pal(length(levels(ys.colfactor_original)), "Dark2")
      } else {
        ys.colors <- rainbow(length(levels(ys.colfactor_original)))
      }
    }
  } else {
    ys.colfactor_original <- rep(1, nrow(ys.output_original))
    ys.colors <- 'black'
  }
  if (!is.null(NUMVAR)) {
    type='n'
    ys.numfactor_original <- ys.output_original[,NUMVAR]
    # ys.numbers contains the palette of distinct numbers
    if (is.null(ys.numbers)) {
      ys.numbers <- as.character(1:length(levels(ys.numfactor)))
    }
  } else {
    type='p'
    ys.numfactor_original <- rep(1, nrow(ys.output_original))
    ys.numbers <- 16
  }
  if (sort.by_size) {
    ys.output_sorted <- ys.output_original %>%
      dplyr::group_by_at(PANELVAR) %>% dplyr::mutate(panel.sample_size = n()) %>% dplyr::ungroup() %>%
      dplyr::arrange(desc(panel.sample_size)) %>% as.data.frame() %>% dplyr::select(-panel.sample_size)
    ys.output_sorted[,PANELVAR] <- factor(ys.output_sorted[,PANELVAR],
                                          levels=unique(ys.output_sorted[,PANELVAR]))
    ys.output_sorted[,NUMVAR] <- factor(ys.output_sorted[,NUMVAR],
                                        levels=unique(ys.output_sorted[,NUMVAR]))
  } else {
    ys.output_sorted <- ys.output_original
  }
  for (NUM in 1:ceiling(length(levels(ys.output_sorted[,PANELVAR]))/panel.size)) {
    idx.vec = panel.size*(NUM-1)+1:panel.size
    ys.output <- ys.output_sorted[which(as.numeric(ys.output_sorted[,PANELVAR]) %in% idx.vec),]
    ys.output[,PANELVAR] <- as.character(ys.output[,PANELVAR])
    ys.output[,PANELVAR] <- factor(ys.output[,PANELVAR],levels=unique(ys.output[,PANELVAR]))
    out.png=paste0(out.png.PREFIX,'_part',NUM,'.png')
    TITLE_str=paste0(TITLE_str.PREFIX,'\n','(part ',NUM,')')

    #ys.output <- subset of ys.output_sorted used for a single panel (changes as loop iterates)
    ys.outcome <- ys.output[,OUTVAR]
    ys.prediction <- ys.output[,PREDVAR]
    ys.panelfactor <- ys.output[,PANELVAR]
    # ys.colors_vec contains the color assignment for every row in the data set
    if (!is.null(COLVAR)) {
      ys.colfactor <- ys.output[,COLVAR]
    } else {
      ys.colfactor <- rep(1, nrow(ys.output))
    }
    # ys.numbers_vec contains the number assignment for every row in the data set
    if (!is.null(NUMVAR)) {
      ys.numfactor <- ys.output[,NUMVAR]
    } else {
      ys.numfactor <- rep(1, nrow(ys.output))
    }
    ys.colors_vec <- ys.colors[ys.colfactor]
    ys.numbers_vec <- ys.numbers[ys.numfactor]

    png(out.png,width=13,height=14,units='in',res=300)
    par(mfrow=mfrow)
    par(mar=c(5,5,5,2)+0.1, oma=c(1,0,2,0))
    for (i in 1:length(levels(ys.panelfactor))) {
      rows_i <- which(as.numeric(ys.panelfactor)==i)
      lim <- get.square_limits(ys.prediction[rows_i], ys.outcome[rows_i])
      l_lim=lim[1]
      u_lim=lim[2]
      MAE <- median(abs(ys.prediction[rows_i]-ys.outcome[rows_i]), na.rm=T)
      MAE_str <- paste0("MAE=",signif(MAE,3))
      PANEL_str <- levels(ys.panelfactor)[i]
      N <- length(which(!is.na(ys.outcome[rows_i]) & !is.na(ys.prediction[rows_i])))
      if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0) {
        COR <- cor(ys.prediction[rows_i],ys.outcome[rows_i])
        COR_str <- paste0("cor=",signif(COR,2))
        MAE_str <- paste0(MAE_str,', ',COR_str)
      }
      plot(y=ys.prediction[rows_i],x=ys.outcome[rows_i],
           type=type,main=paste0(PANEL_str,' (N=',N,')'),
           ylab=PREDVAR,xlab=OUTVAR,pch=ys.numbers_vec[rows_i],
           cex.main=1.5,cex.lab=1.5,cex.axis=1.5,
           col=ys.colors_vec[rows_i],xlim=lim,ylim=lim)
      if (type=='n') {
        text(y=ys.prediction[rows_i],x=ys.outcome[rows_i],cex=0.8,
             labels=ys.numbers_vec[rows_i],col=ys.colors_vec[rows_i])
      }
      if (!is.na(var(ys.outcome[rows_i], na.rm=T)) && var(ys.outcome[rows_i], na.rm=T) > 0 && !is.na(var(ys.prediction[rows_i], na.rm=T))) {
        abline(lm(ys.prediction[rows_i]~ys.outcome[rows_i]))
      }
      abline(0,1,lty="dashed")
      title(MAE_str,outer=F,line=0.4,cex.main=1.5)
      mtext(LETTERS[i],at=l_lim,adj=1,font=2,cex=1.4)
    }
    title(TITLE_str,outer=T,line=-1,cex.main=1.5)
    dev.off()
  }
}

# ### Generates elastic nets by training on a separate dataset
# ### Trains on a single specified subset of covariates
# ### Returns: a numeric equal to the correlation on the testing dataset
# generateCorrSingleSubsetTrained <- function(xs.train,ys.train,subset_idx.vec,xs.test,ys.test,OUTVAR,ALPHA=0.5,NFOLD=10,fun_trans=fun_identity,fun_inv=fun_identity,fun_VAR1=NULL,fun_VAR2=NULL,loglambda.seq=NULL) {
#   #subset_idx.vec=Integer Vector, column indices indicating the random subset
#   #OUTVAR=String, name of dependent variable in ys.train and ys.test
#   #ALPHA=1 for lasso, ALPHA=0 for ridge
#   #fun_trans=Function, transformation applied to OUTVAR before fitting the model
#   #fun_inv=Function, topological inverse of fun_trans
#   #fun_VAR1=String, name of variable to be used as the 1st parameter of fun_trans and fun_inv (default is NULL)
#   #fun_VAR2=String, name of variable to be used as the 2nd parameter of fun_trans and fun_inv (default is NULL)
#   xs.train_sub = xs.train[,subset_idx.vec]
#   xs.test_sub = xs.test[,subset_idx.vec]
#   ys.train_outcome=ys.train[,OUTVAR]
#   ys.test_outcome=ys.test[,OUTVAR]
#   ys.train_funvar1=ys.train[,fun_VAR1] # 0 columns iff VAR=NULL
#   ys.train_funvar2=ys.train[,fun_VAR2]
#   ys.test_funvar1=ys.test[,fun_VAR1]
#   ys.test_funvar2=ys.test[,fun_VAR2]
#   if (is.null(loglambda.seq)) {lambda=NULL} else {lambda=exp(loglambda.seq)}
#   glmnet.Training.CV = cv.glmnet(xs.train_sub,fun_trans(ys.train_outcome,ys.train_funvar1,ys.train_funvar2),
#                                  nfolds=NFOLD,alpha=ALPHA,family="gaussian",type.measure="mse",lambda=lambda)
#   lambda.glmnet.Training = glmnet.Training.CV$lambda.1se
#   glmnet.Training = glmnet(xs.train_sub,fun_trans(ys.train_outcome,ys.train_funvar1,ys.train_funvar2),
#                            alpha=ALPHA,family="gaussian",lambda=lambda,nlambda=100)
#   Y.pred_sub.test_outcome = fun_inv(as.numeric(predict(glmnet.Training,xs.test_sub,type="response",s=lambda.glmnet.Training)),
#                                     ys.test_funvar1,ys.test_funvar2)
#   return(cor(Y.pred_sub.test_outcome, ys.test_outcome))
# }
