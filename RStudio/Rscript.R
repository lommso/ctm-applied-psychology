#--- STARTUP -------------------------------------------------------------------------------------------------------------------------------------------------------------------------#

# PREREQUISITES
setwd("C:/Users/schneluc/OneDrive - adidas/5. Master Program/5. Masterarbeit/implementation/RStudio")
library(ctsemOMX)
library(OpenMx)
require(lavaan)
require(tidyverse)
library(hash)
library(plyr)

# LOAD DATA
data1 = read.csv(file = '../../data/samples/data1.csv') # All anchors (incl. singles & refreshment) and Tpoints with a value for sat6 and per1i6
data2 = read.csv(file = '../../data/samples/data2.csv') # All subjects (incl. partners) in a relationship who have at least 1 non-NA variable of interest
data3 = read.csv(file = '../../data/samples/data3.csv') # Only anchors in a relationship who have at least 1 non-NA variable of interest
data4 = read.csv(file = '../../data/samples/data4.csv') # Same like data1, just without the refreshment sample
data5 = read.csv(file = '../../data/samples/data5.csv') # All subjects (incl. partners) in a relationship who have at least 1 non-NA variable of interest, 
data6 = read.csv(file = '../../data/samples/data6.csv') # All subjects (incl. partners) in a relationship who have at least 1 non-NA variable of interest (revised)
data7 = read.csv(file = '../../data/samples/data7.csv') # All subjects (incl. partners) in a relationship who have at least 1 non-NA variable of interest (revised)




#--- USEFUL FUNCTIONS ---------------------------------------------------------------------------------------------------------------------------------------------------------------#

intervalise <- function(df_wide, manifestNames, startoffset=-1, individualRelativeTime=FALSE) {
  df_intervalised = ctIntervalise(datawide=df_wide, 
                          Tpoints=ncol(df_wide)/(length(manifestNames)+1), 
                          n.manifest=length(manifestNames), 
                          manifestNames=manifestNames, 
                          individualRelativeTime=individualRelativeTime, 
                          startoffset=startoffset)
  return(as.data.frame(df_intervalised))
}

# FUNCTION TO PREPRocess THE DATASET
createWideDataSet <- function(df_long, manifestNames,  timeVar="wave", Tpoints=NULL, minTpoints=2, exclWaves=c(), keepOnlyWaves=NULL, exclNANs=FALSE) {
  
  # Exclude certain waves from the dataset
  exclWaves = if(is.null(keepOnlyWaves)) exclWaves else setdiff(1:11, keepOnlyWaves)
  for(wave in exclWaves) {
    df_long = df_long[df_long[,'wave']!=wave,]
  }
  
  # Keep only subjects with data for at least minTpoints time points
  frequencies = data.frame(table(df_long['id']))
  df_long = merge(x=df_long, y=frequencies, by.x="id", by.y='Var1', all.x = TRUE)
  df_long = df_long[df_long[,'Freq']>=minTpoints,]
  
  # Transform to wide dataset
  df_wide = ctLongToWide(datalong=df_long, 
                         id='id', 
                         time=timeVar, 
                         manifestNames=manifestNames)
  
  # Intervalise for continuous time models
  if(!is.null(Tpoints)) {
    df_wide = intervalise(df_wide, manifestNames)
  }
  
  # Remove NANs
  if(exclNANs) {
    df_wide = df_wide[complete.cases(df_wide), ]
  }
  
  cat('Number of subjects:           ', nrow(df_wide), '\n')
  cat('Average number of time points:', nrow(df_wide), '\n')
  cat('Variance of time intervals:   ', nrow(df_wide), '\n')
  cat('Percentage of NANs:           ', nrow(df_wide), '\n')
  
  return(as.data.frame(df_wide))
}

# FUNCTION TO DEFINE, TRAIN AND ANALYZE A CONTINUOUS TIME MODEL
buildModel <- function(data, manifestNames, Tpoints, timeVar=NULL, TIpredNames=c(), TDpredNames=c(), latentNames=manifestNames, modelType='omx', nSamples=NULL, minTpoints=2, exclWaves=c(), exclNANs=FALSE) {
  
  # Retrieve parameters
  n.manifest = length(manifestNames)
  n.latent = length(latentNames)
  n.TIpred=length(TIpredNames)
  n.TDpred=length(TDpredNames)
  
  # Construct wide data set
  if(!is.null(timeVar)) {
    data = createWideDataSet(data, timeVar, manifestNames, Tpoints, TIpredNames, TDpredNames, minTpoints, exclWaves, exclNANs)
    if(timeVar=='age') {Tpoints=Tpoints+1}
  }
  if(!is.null(nSamples)) {
    data=data[1:nSamples,]
  }
  
  # Define model
  model = ctModel(type=modelType,
                  Tpoints=Tpoints, 
                  n.manifest=n.manifest, 
                  n.TIpred=n.TIpred, 
                  n.TDpred=n.TDpred, 
                  manifestNames=manifestNames, 
                  TIpredNames=TIpredNames, 
                  TDpredNames=TDpredNames, 
                  MANIFESTMEANS = matrix(0, nrow=2, ncol=1),
                  CINT = matrix(paste0("cint_",manifestNames),ncol=1),
                  TRAITVAR = "auto",
                  LAMBDA=diag(length(manifestNames)))
  
  # Fit the model
  start_time = Sys.time()
  fit = ctFit(dat=data, ctmodelobj=model, dataform="wide", retryattempts=10)
  end_time = Sys.time()
  
  return(fit)
}





#--- SECTION 1: TRY SOME SIMPLE CONTINUOUS TIME MODELS  -----------------------------------------------------------------------------------------------------------------------------#

fit1 = buildModel(data1, manifestNames=c("sat6", "per1i6"), Tpoints=8, timeVar='age', minTpoints=1)
fit2 = buildModel(data1, manifestNames=c("sat6", "per1i6"), Tpoints=8, timeVar='age', nSamples=1000)
fit3 = buildModel(data2, manifestNames=c("sat6", "per1i6"), Tpoints=11, timeVar='age', nSamples=1000, minTpoints=1)
fit4 = buildModel(data2, manifestNames=c("sat6", "per1i6"), Tpoints=11, timeVar='age', nSamples=1000)
fit5 = buildModel(data3, manifestNames=c("sat6", "per1i6"), Tpoints=11, timeVar='age', nSamples=1000)
fit6 = buildModel(data1[data1[,'relevant']==1,], manifestNames=c("sat6", "per1i6"), Tpoints=8, timeVar='age')
fit7 = buildModel(data1[data1[,'relevant1']>0,], manifestNames=c("sat6", "per1i6"), Tpoints=8, timeVar='age', nSamples=1000)
fit8 = buildModel(data1[data1[,'relevant2']>0,], manifestNames=c("sat6", "per1i6"), Tpoints=8, timeVar='age', nSamples=1000)
fit9 = buildModel(data1[data1[,'relevant3']>0,], manifestNames=c("sat6", "per1i6"), Tpoints=8, timeVar='age', nSamples=1000)
fit10 = buildModel(data2, manifestNames=c("atts", "sat6"), Tpoints=11, timeVar='age', nSamples=3000)
fit11 = buildModel(data2, manifestNames=c("sat6", "atts"), Tpoints=11, timeVar='age', nSamples=3000)
fit12 = buildModel(data2, manifestNames=c("per1i6", "atto"), Tpoints=11, timeVar='age', nSamples=3000)
fit13 = buildModel(data2, manifestNames=c("per1i6", "atts"), Tpoints=11, timeVar='age', nSamples=3000)
fit14 = buildModel(data2, manifestNames=c("sat6", "patto"), Tpoints=11, timeVar='age', nSamples=3000)
fit15 = buildModel(data2, manifestNames=c("sat6", "patts"), Tpoints=11, timeVar='age', nSamples=3000)
fit16 = buildModel(data2, manifestNames=c("per1i6", "patto"), Tpoints=11, timeVar='age', nSamples=3000)
fit17 = buildModel(data2, manifestNames=c("per1i6", "patts"), Tpoints=11, timeVar='age', nSamples=3000)
fit18 = buildModel(data2, manifestNames=c("atto", "patto"), Tpoints=11, timeVar='age', nSamples=3000)
fit19 = buildModel(data2, manifestNames=c("atto", "patts"), Tpoints=11, timeVar='age', nSamples=3000)
fit20 = buildModel(data2, manifestNames=c("atts", "patto"), Tpoints=11, timeVar='age', nSamples=3000)
fit21 = buildModel(data2, manifestNames=c("atts", "patts"), Tpoints=11, timeVar='age', nSamples=3000)
fit22 = buildModel(data2, manifestNames=c("sat6", "atts", "atto", "per1i6"), Tpoints=11, timeVar='age', nSamples=2000)
fit23 = buildModel(data4, manifestNames=c("sat6", "atts"), Tpoints=8, timeVar='age', nSamples=1000)
fit24 = buildModel(data4, manifestNames=c("atts", "per1i6"), Tpoints=8, timeVar='age', nSamples=1000)
fit25 = buildModel(data4, manifestNames=c("sat6", "per1i6"), Tpoints=8, timeVar='age', nSamples=1000)
fit26 = buildModel(data4, manifestNames=c("atts", "sat6"), Tpoints=8, timeVar='age', nSamples=1000)
fit27 = buildModel(data4, manifestNames=c("sat6", "per1i6"), Tpoints=8, timeVar='age', nSamples=2000)





#--- SECTION 2: TRY SOME VSIMPLE DISCRETE TIME MODELS --------------------------------------------------------------------------------------------------------------------------------#

data = createWideDataSet(data2, exclWaves=c(2,3,6), timeVar="wave", manifestNames=c("sat6", "per1i6"))
mod = readChar("discrete_models/mod28.txt", file.info("discrete_models/mod28.txt")$size) %>% str_replace_all("\r", "")
fit28 <- lavaan(mod, data[complete.cases(data),], missing='ML', int.ov.free=F, int.lv.free=F, auto.fix.first=F, auto.fix.single=F, auto.cov.lv.x=F, auto.cov.y=F, auto.var=F)

data = createWideDataSet(data2, exclWaves=c(2,3,6), timeVar="wave", manifestNames=c("sat6", "per1i6"))
mod = readChar("discrete_models/mod29.txt", file.info("discrete_models/mod29.txt")$size) %>% str_replace_all("\r", "")
fit29 <- lavaan(mod, data[complete.cases(data[, c('sat6_T0', 'sat6_T1', 'sat6_T2', 'per1i6_T0', 'per1i6_T1', 'per1i6_T2')]),], missing='ML', int.ov.free=F, int.lv.free=F, auto.fix.first=F, auto.fix.single=F, auto.cov.lv.x=F, auto.cov.y=F, auto.var=F)

data = createWideDataSet(data2, exclWaves=c(2,3,6), timeVar="wave", manifestNames=c("sat6", "per1i6"))
mod = readChar("discrete_models/mod29.txt", file.info("discrete_models/mod29.txt")$size) %>% str_replace_all("\r", "")
fit30 = lavaan(mod, data[complete.cases(data),], missing='ML', int.ov.free=F, int.lv.free=F, auto.fix.first=F, auto.fix.single=F, auto.cov.lv.x=F, auto.cov.y=F, auto.var=F)

data = createWideDataSet(data2, exclWaves=c(2,3,6), timeVar="wave", manifestNames=c("sat6", "per1i6"))
mod = readChar("discrete_models/mod31.txt", file.info("discrete_models/mod31.txt")$size) %>% str_replace_all("\r", "")
fit31 <- lavaan(mod, data[complete.cases(data[, c('sat6_T0', 'sat6_T1', 'sat6_T2', 'sat6_T3', 'sat6_T4', 'per1i6_T0', 'per1i6_T1', 'per1i6_T2', 'per1i6_T3', 'per1i6_T4')]),], missing='ML', int.ov.free=F, int.lv.free=F, auto.fix.first=F, auto.fix.single=F, auto.cov.lv.x=F, auto.cov.y=F, auto.var=F)





#--- SECTION X: BAYESIAN MODEL ------------------------------------------------------------------------------------------------------------------------------------------------------#

# BAYESIAN MODEL

manifestNames=c("sat6", "atto")
data_long = raw_data[1:100, c("id", "age", "sat6", "atto")]
data_long[, c("sat6", "atto")] = scale(data_long[, c("sat6", "atto")])
colnames(data_long) = c("id", "time", "sat6", "atto")
n.manifest = length(c("sat6", "atto"))
model = ctModel(type='stanct', 
                n.latent=n.manifest, 
                n.manifest=n.manifest,
                manifestNames=c("sat6", "atto"), 
                latentNames=c("sat6", "atto"), 
                LAMBDA=diag(n.manifest))
fit = ctStanFit(datalong=data_long, ctstanmodel=model, optimize=FALSE, nopriors=FALSE, chains = 2, iter=300, plot=TRUE)
summary(fit)





#--- SECTION 3: COMPARE LINEAR REGRESSION, DISCRETE TIME, AND CONTINUOUS TIME MODELS ------------------------------------------------------------------------------------------------#

# Function to calculate the RMSE
calcRMSE <- function(estimates, data) {
  
  # Retrieve Variable Names
  var1 = colnames(estimates)[1]
  var2 = colnames(estimates)[2]
  
  rmse <- c()
  
  for(row in 1:nrow(estimates)) {
    
    # Make predictions
    est1 = estimates[row,][var1]
    est2 = estimates[row,][var2]
    data1 = data[,c(paste(var1,'_T0',sep=''))]
    data2 = data[,c(paste(var2,'_T0',sep=''))]
    int = estimates[row,]['intcp']
    y_pred = est1*data1 + est2*data2 + int
    
    # Compute RMSE
    y_true = data[,c(paste(rownames(estimates)[row],'_T1',sep=''))]
    rmse = append(rmse, sqrt(mean((y_pred-y_true)^2)))
  }
  
  return(rmse)
}

# Function to transform continuous time parameters to discrete time parameters
contToDisc <- function(A_cont=NULL, b_cont=NULL, b_disc=NULL, manifestNames=NULL, estimates=NULL, interval=1) {
  
  manifestNames= if(is.null(manifestNames)) rownames(estimates) else manifestNames
  m = length(manifestNames)
  A_cont = if(is.null(A_cont)) estimates[,(1:m)] else A_cont
  b_cont = if(is.null(b_cont)) estimates[,m+1] else b_cont
  b_disc = if(ncol(estimates)>m+1) estimates[,m+2] else if(!is.null(b_disc)) b_disc else numeric(m)
  
  A_disc = expm(A_cont*interval)
  b_disc = solve(A_cont) %*% (A_disc-diag(1,m)) %*% b_cont + b_disc
  estimates = matrix(cbind(A_disc, b_disc), nrow=m, dimnames=list(manifestNames, c(manifestNames, 'intcp')))
  
  return(estimates)
}

# Data preparation
manifestNames=c("sat6", "per1i6")
data = as.data.frame(createWideDataSet(data5, Tpoints=2, exclWaves=c(1,2,3,6,7,8,9,10,11), minTpoints=2, timeVar="wave", manifestNames=manifestNames, exclNANs=TRUE))
set.seed(123)
train_ind = sample(seq_len(nrow(data)), size = floor(0.7 * nrow(data)))
train = data[train_ind, ]
test = data[-train_ind, ]
results = matrix(nrow=0, ncol=8, dimnames=list(c(), c('auto_sat6', 'cross_per1i6>sat6', 'intcp_sat6', 'cross_sat6>per1i6', 'auto_per1i6', 'intcp_per1i6', 'train_rmse', 'test_rmse')))

# Naive Model
estimates_naive = matrix(c(1,0,0,0,1,0), ncol=3, byrow=TRUE, dimnames=list(manifestNames, c(manifestNames, 'intcp')))
rmse_naive_train = calcRMSE(estimates_naive, data=train)
rmse_naive_test = calcRMSE(estimates_naive, data=test)
results = rbind(results, matrix(nrow=1, c(t(estimates_naive), mean(rmse_naive_train), mean(rmse_naive_test)), dimnames=list('Naive Model', c())))

# Linear Regression
fit_lr1 = lm(sat6_T1 ~ 1 + sat6_T0 + per1i6_T0, data=train)
fit_lr2 = lm(per1i6_T1 ~ 1 + sat6_T0 + per1i6_T0, data=train)
estimates_lr = matrix(c(fit_lr1$coefficients[c(2,3,1)],fit_lr2$coefficients[c(2,3,1)]), ncol=3, dimnames=list(manifestNames, c(manifestNames, 'intcp')), byrow=TRUE)
rmse_lr_train = c(sqrt(mean((predict(fit_lr1, train)-train$sat6_T1)^2)), sqrt(mean((predict(fit_lr2, train)-train$per1i6_T1)^2)))
rmse_lr_test = c(sqrt(mean((predict(fit_lr1, test)-test$sat6_T1)^2)), sqrt(mean((predict(fit_lr2, test)-test$per1i6_T1)^2)))
results = rbind(results, matrix(nrow=1, c(t(estimates_lr), mean(rmse_lr_train), mean(rmse_lr_test)), dimnames=list('Linear Regression', c())))

# Discrete Cross Lagged
mod = '
  # Regressions
  sat6_T1 ~ a*sat6_T0 + b*per1i6_T0
  sat6_T1 ~ 1
  per1i6_T1 ~ c*sat6_T0 + d*per1i6_T0
  per1i6_T1 ~ 1 
  
  # Correlations
  sat6_T1 ~~ per1i6_T1'
fit_dCL = sem(mod, data=train)
estimates_dCL = matrix(coef(fit_dCL)[1:6], ncol=3, byrow=TRUE, dimnames=list(manifestNames, c(manifestNames, 'intcp')))
rmse_dCL_train = calcRMSE(estimates_dCL, data=train)
rmse_dCL_test = calcRMSE(estimates_dCL, data=test)
results = rbind(results, matrix(nrow=1, c(t(estimates_dCL), mean(rmse_dCL_train), mean(rmse_dCL_test)), dimnames=list('Discrete Cross Lagged', c())))

# Continuous Cross Lagged
fit_cCL = buildModel(train, manifestNames=manifestNames, Tpoints=2, minTpoints=1)
estimates_cCL = contToDisc(A_cont=summary(fit_cCL)$DRIFT, b_cont=summary(fit_cCL)$MANIFESTMEANS, manifestNames=manifestNames)
rmse_cCL_train = calcRMSE(estimates_cCL, data=train)
rmse_cCL_test = calcRMSE(estimates_cCL, data=test)
results = rbind(results, matrix(nrow=1, c(t(estimates_cCL), mean(rmse_cCL_train), mean(rmse_cCL_test)), dimnames=list('Continuous Cross Lagged', c())))

# Compare
results[,c(1,2,4,5,3,6,7,8)]

summary(fit_cCL)


#--- SECTION 4: IMPLEMENT FUNCTION FOR TRADITIONAL VOELKLE-MODEL --------------------------------------------------------------------------------------------------------------------#

# Algebra function 1
EXP_Function <- function(number, value) {
  eval(substitute(mxAlgebra(EVEC %*% (exp(replace2 %x% EVA) - tempa)%*% solve(EVEC), 
                            name = paste("EXPd", replace1, sep="")),
                  list(replace1 = number, replace2 = value)))
}

# Algebra function 2
int_Function <- function(number2, value2) {
  eval(substitute(mxAlgebra(solve(DRIFT)%*%((EVEC %*% (exp(replace4 %x% EVA) - tempa)%*% solve(EVEC))-II)%*%t(CINT), 
                            name = paste("intd", replace3, sep="")),
                  list(replace3 = number2, replace4 = value2)))
}

# Algebra function 3
Qd_Function <- function(number3, value3) {
  eval(substitute(mxAlgebra(((solve(DRIFTHATCH)%*%((EVECH%*%(exp(replace6%x%EVAH)-tempb)%*%solve(EVECH))-(II%x%II))%*%rvectorize(Q))) ,
                            name = paste("Qd", replace5, sep="")),
                  list(replace5 = number3, replace6 = value3)))
}

# Voelkle Model from 2012
voelkleModel <- function(n.manifest, n.latent, PHI1, latentM1, manifestM, LAMBDA, THETA, DRIFT, CINT, Q, delta_t, Tpoints, data, measurements) {
  
  #######################################
  #######################################
  ### MATRIX SPECIFICATION FOR OPENMX ###
  #######################################
  #######################################
  
  lambda.differ	<- FALSE
  theta.differ	<- FALSE
  
  
  ###################
  ### 1. A matrix ###
  ###################
  
  Avalues <- matrix(0,nrow=(n.manifest*Tpoints+n.latent*Tpoints), ncol=(n.manifest*Tpoints+n.latent*Tpoints))
  
  for (i in seq(n.latent,n.latent*(Tpoints-1),n.latent)){
    Avalues[(1+i):(i+n.latent), ((1+i)-n.latent):i] <- NA
  }
  
  for (j in 1:Tpoints){
    Avalues[((1+n.latent*Tpoints-n.manifest)+j*n.manifest):(n.latent*Tpoints+j*n.manifest), (j*n.latent-n.latent+1):(j*n.latent)] <- LAMBDA
  }
  
  Afree <-  ((Avalues != 0) & (Avalues != 1) & is.na(Avalues) == F)  
  
  DRIFTfree <- ((DRIFT != 0) & (DRIFT != 1))
  
  Alabels <- matrix(,nrow=(n.manifest*Tpoints+n.latent*Tpoints), ncol=(n.manifest*Tpoints+n.latent*Tpoints))
  
  for (i in 1:(Tpoints-1)){
    temp <- matrix((paste("EXPd", i,  sep = "")), nrow=n.latent, ncol=n.latent, byrow=TRUE)	
    BETA_labels <- temp 
    for (k in 1:n.latent){
      for(j in 1:n.latent){
        BETA_labels [k,j] = paste(BETA_labels[k,j], "[",k,",",j,"]",sep="") 
      }}
    
    Alabels[(1+i*n.latent):(n.latent+i*n.latent), (i*n.latent-n.latent+1):(i*n.latent)] <- BETA_labels
    
  }
  
  templambda <- matrix((paste("lambda", 1:(n.manifest*n.latent), sep = "")), nrow=n.manifest, ncol=n.latent, byrow=FALSE)
  
  for (l in 1:Tpoints){
    if(lambda.differ) LAMBDA_labels <- paste(templambda, "_", l, sep="")
    else LAMBDA_labels <- templambda 
    Alabels[((1+n.latent*Tpoints-n.manifest)+l*n.manifest):(n.latent*Tpoints+l*n.manifest), (l*n.latent-n.latent+1):(l*n.latent)] <- LAMBDA_labels
  } 
  
  ###################
  ### 2. S matrix ###
  ###################
  
  Svalues <- matrix(0,nrow=(n.manifest*Tpoints+n.latent*Tpoints), ncol=(n.manifest*Tpoints+n.latent*Tpoints))
  
  for(i in 1:Tpoints){
    Svalues[(i*n.latent-n.latent+1):(i*n.latent), (i*n.latent-n.latent+1):(i*n.latent)] <- NA
  }
  
  Svalues[1:n.latent, 1:n.latent] <- PHI1	
  
  for(j in 1:Tpoints){
    Svalues[((n.latent*Tpoints+1)+j*n.manifest-n.manifest):((n.latent*Tpoints)+j*n.manifest), ((n.latent*Tpoints+1)+j*n.manifest-n.manifest):((n.latent*Tpoints)+j*n.manifest)] <- THETA
  }
  
  Sfree <- ((Svalues != 0) & (Svalues != 1) & is.na(Svalues) == F)
  Qfree <-  ((Q != 0) & (Q != 1))
  
  Slabels <- matrix(,nrow=(n.manifest*Tpoints+n.latent*Tpoints), ncol=(n.manifest*Tpoints+n.latent*Tpoints))
  
  for (i in 1:Tpoints){
    temp <- matrix((paste("Qd", i-1,  sep = "")), nrow=n.latent, ncol=n.latent, byrow=TRUE)
    PSI_labels <- temp 
    Qd.ind <- matrix(1:n.latent**2,n.latent,n.latent)
    for (k in 1:n.latent){
      for(j in 1:n.latent){
        PSI_labels [k,j] = paste(PSI_labels[k,j], "[",Qd.ind[k,j],",",1,"]",sep="") 
      }
    }
    
    Slabels[(i*n.latent-n.latent+1):(i*n.latent), (i*n.latent-n.latent+1):(i*n.latent)] <- PSI_labels
  }
  
  tempphi <- matrix(paste("phi",1:n.latent,matrix(1:n.latent,n.latent,n.latent,byrow=T),sep=""),n.latent,n.latent)
  
  for(i in 1:n.latent){
    for(j in 1:n.latent){
      if(Svalues[i,j]==Svalues[j,i]) tempphi[i,j]=tempphi[j,i]
    }
  }
  
  Slabels[1:n.latent, 1:n.latent] <- tempphi
  
  temptheta <- matrix(paste("theta",1:n.manifest,matrix(1:n.manifest,n.manifest,n.manifest,byrow=T),sep=""),n.manifest,n.manifest)
  
  for (j in 1:Tpoints){
    if(theta.differ) THETA_labels<- paste(temptheta, "_", j, sep="")
    else THETA_labels <- temptheta 
    Slabels[((n.latent*Tpoints+1)+j*n.manifest-n.manifest):((n.latent*Tpoints)+j*n.manifest), ((n.latent*Tpoints+1)+j*n.manifest-n.manifest):((n.latent*Tpoints)+j*n.manifest)] <- THETA_labels
  }
  
  ###################
  ### 3. F matrix ###
  ###################
  
  Fvalues 	 <- cbind(matrix(0,nrow=(n.manifest*Tpoints), ncol = n.latent*Tpoints), diag(1,(n.manifest*Tpoints)))
  Fnamesy 	 <- c(c(paste("F", 1:(n.latent*Tpoints), sep="")), c(paste("V", 1:(n.manifest*Tpoints), sep="")))
  Fnamesx		 <- c(paste("V", 1:(n.manifest*Tpoints), sep=""))
  
  ###################
  ### 4. M matrix ###
  ###################
  
  Mvalues <- matrix(,nrow=(n.manifest*Tpoints+n.latent*Tpoints), ncol=1)
  
  for (i in 1:Tpoints){
    Mvalues[(i*n.latent-n.latent+1):(i*n.latent), 1] <- NA
  }
  
  Mvalues[1:n.latent, 1] <- latentM1	# first time point
  
  for (j in 1:Tpoints){
    Mvalues[((n.latent*Tpoints+1)+j*n.manifest-n.manifest):((n.latent*Tpoints)+j*n.manifest), 1] <- manifestM
  }
  
  Mfree <- ((Mvalues != 0) & (Mvalues != 1) & is.na(Mvalues)==F)
  
  CINTfree <- ((CINT!= 0))
  
  Mlabels <- matrix(,nrow=(n.manifest*Tpoints+n.latent*Tpoints), ncol=1)
  
  for (i in 1:(Tpoints-1)){
    M_labels <- matrix((paste("intd", i, sep = "")), nrow=n.latent, ncol=1, byrow=TRUE)
    
    for(j in 1:n.latent){
      M_labels [j,1] = paste(M_labels[j,1], "[",j,",",1,"]",sep="") 
    }
    
    Mlabels[(1+i*n.latent):(n.latent+i*n.latent), 1] <- M_labels
  }
  
  Mlabels[1:n.latent, 1] <- matrix((paste("m", 1:n.latent, sep = "")), nrow=n.latent, ncol=1, byrow=TRUE)	
  
  for (j in 1:Tpoints){
    Mlabels[((1+n.latent*Tpoints-n.manifest)+j*n.manifest):(n.latent*Tpoints+j*n.manifest), 1] <- NA
  }
  
  ###################################
  ### Remaining labels for Output ###
  ###################################
  
  CINTlabels 	<- paste("cint",1:n.latent,sep="")
  
  DRIFTlabels 	<- matrix(,n.latent,n.latent)
  
  for (i in 1:n.latent){
    for (j in 1:n.latent){
      DRIFTlabels[i,j] = paste("F",i,j,sep="")
    }			
  }
  
  Qlabels 	<- matrix(,n.latent,n.latent)
  
  for (i in 1:n.latent){
    for (j in i:n.latent){
      Qlabels[j,i] = paste("q",i,j,sep="")
    }
    Qlabels[i,j] <- Qlabels[j,i]			
  }		
  
  ##################
  ##################
  ### RUN OPENMX ###
  ##################
  ##################
  
  colnames(data)=Fnamesx	# rename columns of dataset
  
    model=mxModel("EDM", type="RAM",
                  mxData(observed=data, type="raw"),
                  mxMatrix(type="Full", labels=DRIFTlabels, values=DRIFT, byrow=TRUE, free=DRIFTfree, name="DRIFT"), 
                  mxMatrix(type="Full", labels=CINTlabels, values=CINT, free=CINTfree, name="CINT"),
                  mxMatrix(type="Full", labels=Qlabels, values=Q, byrow=TRUE, free=Qfree, name="Q"),
                  mxMatrix(type="Iden", nrow=n.latent, ncol=n.latent, free=FALSE, name="II"),
                  
                  mxMatrix(
                    values=Avalues,
                    free=Afree,
                    labels=Alabels,
                    name="A"),
                  mxMatrix(
                    values=Svalues,
                    free=Sfree,
                    labels=Slabels,
                    name="S"),
                  mxMatrix(
                    values=Fvalues,
                    free=FALSE,
                    dimnames=list(Fnamesx, Fnamesy),
                    name="F"),
                  mxMatrix( 
                    free=t(Mfree),
                    values=t(Mvalues),
                    labels=t(Mlabels),
                    dimnames=list(1, Fnamesy),
                    name="M"),
                  
                  ############################################################
                  ### different intervals for all timepoints (see delta_t) ###
                  ############################################################
                  
                  ############### drift matrix (Alabels) ###############
                  
                  mxAlgebra(vec2diag(eigenval(DRIFT)), name = "EVA"),
                  mxAlgebra(eigenvec(DRIFT), name = "EVEC"),
                  mxMatrix("Full", values=(matrix(1,n.latent,n.latent)-diag(n.latent)), name = "tempa"),
                  
                  EXP_algebras <- mapply(EXP_Function, measurements, delta_t),		# see Algebra 1)
                  
                  ############### intercepts (Mlabels) ###############
                  
                  int_algebras <- mapply(int_Function, measurements, delta_t),		# see Algebra 2)
                  
                  ############### Q matrix (Slabels)  ###############
                  
                  mxAlgebra(DRIFT%x%II + II%x%DRIFT, name = "DRIFTHATCH"),
                  mxAlgebra(vec2diag(eigenval(DRIFTHATCH)), name = "EVAH"),
                  mxAlgebra(eigenvec(DRIFTHATCH), name = "EVECH"),
                  mxMatrix("Full", values=(matrix(1,n.latent**2,n.latent**2)-diag(n.latent**2)), name = "tempb"),
                  
                  Qd_algebras <- mapply(Qd_Function,measurements,delta_t),			# see Algebra 3)
                  
                  ############### matrices in the model ###############
                  
                  mxRAMObjective("A","S","F","M")
    )
    
    model2=mxModel(model,EXP_algebras,int_algebras,Qd_algebras)
  
  #fit=mxRun(model2)
  
  return(model2)
}

# Fit a voelkle model
n.manifest=2
n.latent=2
Tpoints=2
A_disc = estimates_dCL[1:2,1:2]
mod_voelkle = voelkleModel(n.manifest   = n.manifest,                                                   # total number of indicators for ALL factors. 
                           n.latent     = n.latent,                                                     # number of latent variables
                           PHI1         = cov(data[,1:2], use = "pairwise.complete.obs"),               # var/cov matrix of latent variables at first time point
                           latentM1     = matrix(c(colMeans(train[,1:2])),nrow=n.latent, ncol=1),       # means of latent variables at first time point
                           manifestM    = matrix(c(0,0),nrow=n.manifest, ncol=1, byrow=TRUE),           # intercepts of manifest variables (usually fixed to zero)
                           LAMBDA       = matrix(c(1,0,0,1),nrow=n.manifest, ncol=n.latent,byrow=TRUE), # factor loading matrix
                           THETA        = matrix(c(0,0,0,0),nrow=n.manifest, ncol=n.manifest),          # var/cov matrix of measurement error
                           DRIFT        = A_disc-diag(1,2),                             # drift matrix
                           CINT         = matrix(c(estimates_dCL[1:2,3]),ncol=n.latent, nrow=1, byrow=TRUE),      # continuous time intercepts
                           Q            = t(chol(round(A_disc,3))) %*% chol(round(A_disc,3)),           # diffusion matrix	
                           delta_t      = c(1),                                                         # vector of (possibly different) time intervals
                           data         = train[,1:4],
                           Tpoints      = 2, 
                           measurements = 1:(Tpoints-1)) 
fit_voelkle=mxRun(mod_voelkle)
estimates_voelkle = contToDisc(A_cont = matrix(summary(fit_voelkle)$parameters[5:8,'Estimate'], ncol=2), b_cont = matrix(summary(fit_voelkle)$parameters[9:10,'Estimate'], ncol=1), manifestNames=manifestNames)
rmse_voelkle_train = calcRMSE(estimates_voelkle, data=train)
rmse_voelkle_test = calcRMSE(estimates_voelkle, data=test)
results = rbind(results, matrix(nrow=1, c(t(estimates_voelkle), mean(rmse_voelkle_train), mean(rmse_voelkle_test)), dimnames=list('Traditional CT Model', c())))

# Compare
results




#--- SECTION 5: IMPLEMENT CROSS-VALIDATION ------------------------------------------------------------------------------------------------------------------------------------------#

# Data preparation
set.seed(123)
manifestNames=c("sat6", "per1i6")
data = as.data.frame(createWideDataSet(data5, Tpoints=2, exclWaves=c(1,2,3,6,7,8,9,10,11), minTpoints=2, timeVar="wave", manifestNames=manifestNames, exclNANs=TRUE))
data<-data[sample(nrow(data)),]  # Randomly shuffle the data
k=5
folds <- cut(seq(1, nrow(data)), breaks=k, labels=FALSE) # Create 5 equally size folds
results = matrix(nrow=0, ncol=9, dimnames=list(c(), c('fold', 'auto_sat6', 'cross_per1i6_sat6', 'intcp_sat6', 'cross_sat6>per1i6', 'auto_per1i6', 'intcp_per1i6', 'train_rmse', 'test_rmse')))

# Perform 10-fold CV
for(fold in 1:k){
  
  # Retrieve data for current fold
  testIndexes = which(folds==fold, arr.ind=TRUE)
  test        = data[testIndexes, ]
  train       = data[-testIndexes, ]
  
  # Naive Model
  estimates_naive = matrix(c(1,0,0,0,1,0), ncol=3, byrow=TRUE, dimnames=list(manifestNames, c(manifestNames, 'intcp')))
  rmse_naive_train = calcRMSE(estimates_naive, data=train)
  rmse_naive_test = calcRMSE(estimates_naive, data=test)
  results = rbind(results, matrix(nrow=1, c(fold, t(estimates_naive), mean(rmse_naive_train), mean(rmse_naive_test)), dimnames=list('Naive Model', c())))
  
  # Linear Regression
  fit_lr1 = lm(sat6_T1 ~ 1 + sat6_T0 + per1i6_T0, data=train)
  fit_lr2 = lm(per1i6_T1 ~ 1 + sat6_T0 + per1i6_T0, data=train)
  estimates_lr = matrix(c(fit_lr1$coefficients[c(2,3,1)],fit_lr2$coefficients[c(2,3,1)]), ncol=3, dimnames=list(manifestNames, c(manifestNames, 'intcp')), byrow=TRUE)
  rmse_lr_train = c(sqrt(mean((predict(fit_lr1, train)-train$sat6_T1)^2)), sqrt(mean((predict(fit_lr2, train)-train$per1i6_T1)^2)))
  rmse_lr_test = c(sqrt(mean((predict(fit_lr1, test)-test$sat6_T1)^2)), sqrt(mean((predict(fit_lr2, test)-test$per1i6_T1)^2)))
  results = rbind(results, matrix(nrow=1, c(fold, t(estimates_lr), mean(rmse_lr_train), mean(rmse_lr_test)), dimnames=list('Linear Regression', c())))
  
  # Discrete Cross Lagged
  mod = '
  # Regressions
  sat6_T1 ~ a*sat6_T0 + b*per1i6_T0
  sat6_T1 ~ 1
  per1i6_T1 ~ c*sat6_T0 + d*per1i6_T0
  per1i6_T1 ~ 1 
  
  # Correlations
  sat6_T1 ~~ per1i6_T1'
  fit_dCL = sem(mod, data=train)
  estimates_dCL = matrix(coef(fit_dCL)[1:6], ncol=3, byrow=TRUE, dimnames=list(manifestNames, c(manifestNames, 'intcp')))
  rmse_dCL_train = calcRMSE(estimates_dCL, data=train)
  rmse_dCL_test = calcRMSE(estimates_dCL, data=test)
  results = rbind(results, matrix(nrow=1, c(fold, t(estimates_dCL), mean(rmse_dCL_train), mean(rmse_dCL_test)), dimnames=list('Discrete Cross Lagged', c())))
  
  # Continuous Cross Lagged
  fit_cCL = buildModel(train, manifestNames=manifestNames, Tpoints=2, minTpoints=1)
  estimates_cCL = contToDisc(A_cont=summary(fit_cCL)$DRIFT, b_cont=summary(fit_cCL)$CINT, manifestNames=manifestNames)
  rmse_cCL_train = calcRMSE(estimates_cCL, data=train)
  rmse_cCL_test = calcRMSE(estimates_cCL, data=test)
  results = rbind(results, matrix(nrow=1, c(fold, t(estimates_cCL), mean(rmse_cCL_train), mean(rmse_cCL_test)), dimnames=list('Continuous Cross Lagged', c())))

  # Voelkle model
  n.manifest=2
  n.latent=2
  Tpoints=2
  A_disc = estimates_dCL[1:2,1:2]
  mod_voelkle = voelkleModel(n.manifest   = n.manifest,                                                   # total number of indicators for ALL factors. 
                             n.latent     = n.latent,                                                     # number of latent variables
                             PHI1         = cov(data[,1:2], use = "pairwise.complete.obs"),               # var/cov matrix of latent variables at first time point
                             latentM1     = matrix(c(colMeans(train[,1:2])),nrow=n.latent, ncol=1),       # means of latent variables at first time point
                             manifestM    = matrix(c(0,0),nrow=n.manifest, ncol=1, byrow=TRUE),           # intercepts of manifest variables (usually fixed to zero)
                             LAMBDA       = matrix(c(1,0,0,1),nrow=n.manifest, ncol=n.latent,byrow=TRUE), # factor loading matrix
                             THETA        = matrix(c(0,0,0,0),nrow=n.manifest, ncol=n.manifest),          # var/cov matrix of measurement error
                             DRIFT        = A_disc-diag(1,2),                             # drift matrix
                             CINT         = matrix(c(estimates_dCL[1:2,3]),ncol=n.latent, nrow=1, byrow=TRUE),      # continuous time intercepts
                             Q            = t(chol(round(A_disc,3))) %*% chol(round(A_disc,3)),           # diffusion matrix	
                             delta_t      = c(1),                                                         # vector of (possibly different) time intervals
                             data         = train[,1:4],
                             Tpoints      = 2, 
                             measurements = 1:(Tpoints-1)) 
  fit_voelkle=mxRun(mod_voelkle)
  estimates_voelkle = contToDisc(A_cont = matrix(summary(fit_voelkle)$parameters[5:8,'Estimate'], ncol=2), b_cont = matrix(summary(fit_voelkle)$parameters[9:10,'Estimate'], ncol=1), manifestNames=manifestNames)
  rmse_voelkle_train = calcRMSE(estimates_voelkle, data=train)
  rmse_voelkle_test = calcRMSE(estimates_voelkle, data=test)
  results = rbind(results, matrix(nrow=1, c(fold, t(estimates_voelkle), mean(rmse_voelkle_train), mean(rmse_voelkle_test)), dimnames=list('Traditional CT Model', c())))
}

# Compare
results_mean = aggregate(results, list(row.names(results)), mean)
results_mean[order(results_mean$test_rmse),]





#--- SECTION 6: FUNCTIONS FOR FLEXIBLE EXPERIMENTS ------------------------------------------------------------------------------------------------------------------------------------------#

            createModel = function(data, manifestNames, Tpoints, type, hyperparams) {
              
              if(type=='simple') {
                mod = ''
                for(var in manifestNames) {
                  for(T in 1:(Tpoints-1)) {
                    target = paste(var, T, sep='_T')
                    inputs = paste0(var, '_T', T-1)
                    if(sum(complete.cases(data[,c(target, inputs)])) > 0) { # Only add equations if at least 1 line with data for all involved variables
                      regr = paste0(target,' ~ ', paste0(paste0('auto_', var), ' * ', inputs))
                      intcp = paste0(target,' ~ intcp_',var,'*1')
                      mod = paste(mod, regr, intcp, sep=' \n ' )
                    }
                  }
                }
                return(mod)
              }
              
              else if(type=='lm' && Tpoints==2) {
                mod = c()
                for(var in manifestNames) {
                  mod = append(mod, as.formula(paste0(paste0(var,'_T1'), ' ~ ', paste(lapply(manifestNames, paste0, "_T0"), collapse=' + '))))
                }
                return(mod)
              }
              
              else if(type=='discrete') {
                mod = ''
                for(var in manifestNames) {
                  for(T in 1:(Tpoints-1)) {
                    target = paste(var, T, sep='_T')
                    inputs = paste0(manifestNames, '_T', T-1)
                    if(sum(complete.cases(data[,c(target, inputs)])) > 0) { # Only add equations if at least 1 line with data for all involved variables
                      regr = paste0(target,' ~ ',paste(paste0(var,'_',manifestNames, '*' ,inputs), collapse=' + '))
                      intcp = paste0(target,' ~ intcp_',var,'*1')
                      mod = paste(mod, regr, intcp, sep=' \n ' )
                    }
                  }
                }
                for(i in 1:length(manifestNames)) {
                  for(j in 1:length(manifestNames)) {
                    if(j>i) {
                      for(T in 1:(Tpoints-1)) {
                        var_i = paste0(manifestNames[i],'_T',T)
                        var_j = paste0(manifestNames[j],'_T',T)
                        if(grepl(paste(var_i,'~ '), mod, fixed=TRUE) && grepl(paste(var_j,'~ '), mod, fixed=TRUE)) { # Only add correlations for variables that have been defined as target earlier
                          mod = paste(mod, paste(paste(var_i,'~~',var_j), collapse= ' \n '), sep=' \n ')
                        }
                      }
                    }
                  }
                }
                return(mod)
              }
              
              else if(type=='voelkle') {
                # The definition of the Voelkle model requires a data object. Therefore it is part of the 'fit' step
                return(NULL)
              }
              
              else if(type=='ctsem') {
                hyperparams$cint = (if (is.null(hyperparams$cint)) matrix(paste0("cint_",manifestNames),ncol=1) else if (hyperparams$cint=='auto') 'auto' else matrix(paste0("cint_",manifestNames),ncol=1))
                hyperparams$manifestmeans = (if (is.null(hyperparams$manifestmeans)) matrix(0, nrow=length(manifestNames), ncol=1) else if (hyperparams$manifestmeans=='auto') 'auto' else matrix(0, nrow=length(manifestNames), ncol=1))
                hyperparams$traitvar = (if (is.null(hyperparams$traitvar)) 'auto' else if (hyperparams$traitvar=='null') 'null' else 'auto')
                mod = ctModel(type='omx', 
                              Tpoints=Tpoints,
                              n.manifest=length(manifestNames),
                              manifestNames=manifestNames, 
                              MANIFESTMEANS = hyperparams$manifestmeans,
                              CINT = hyperparams$cint,
                              TRAITVAR = hyperparams$traitvar,
                              LAMBDA=diag(length(manifestNames)))
                return(mod)
              }
            }
            
            fit = function(mod, data, manifestNames, type, hyperparams=list(), startVals=NULL) {
              
              m = length(manifestNames)
              
              if (type=='simple') {
                if (is.null(hyperparams$missing)) { hyperparams$missing = 'listwise' }
                fit = sem(mod, data=data, missing=hyperparams$missing)
                coefs = coef(fit)[!duplicated(names(coef(fit)))]
                estimates = matrix(cbind(diag(coefs[(1:m)*2-1]), coefs[(1:m)*2]), nrow=m, dimnames=list(manifestNames, c(manifestNames, 'intcp')))
                return(estimates)
              }
              
              if (type=='discrete') {
                if (is.null(hyperparams$missing)) { hyperparams$missing = 'listwise' }
                fit = sem(mod, data=data, missing=hyperparams$missing)
                coefs = coef(fit)[!duplicated(names(coef(fit)))]
                estimates = matrix(coefs[1:(m*(m+1))], nrow=m, byrow=TRUE, dimnames=list(manifestNames, c(manifestNames, 'intcp')))
                return(estimates)
              }
              
              else if (type=='lm') {
                estimates = matrix(nrow=0, ncol=m+1, dimnames=list(c(), c(manifestNames, 'intcp')))
                for(i in 1:m) {
                  fit = lm(mod[[i]], data=data)
                  estimates = rbind(estimates, matrix(c(fit$coefficients[c(2:(m+1),1)]), nrow=1, dimnames=list(manifestNames[i], c())))
                }
                return(estimates)
              }
              
              else if (type=='ctsem') {
                retryattempts = if (is.null(hyperparams$retryattempts)) 9 else hyperparams$retryattempts
                fit = ctFit(ctmodelobj=mod, dat=data, dataform="wide", retryattempts=retryattempts)
                res = summary(fit)
                estimates = matrix(cbind(res$DRIFT, res$CINT), nrow=m, dimnames=list(manifestNames, c(manifestNames, 'intcp')))
                return(estimates)
              }
              
              else if (type=='voelkle') {
                m = length(manifestNames)
                l = m #  number of latentvariables
                Tpoints = ceiling(ncol(data)/(length(manifestNames)+1))
                intervals = c(as.matrix(data[,(m*Tpoints+1):ncol(data)]))
                delta_t = mean(intervals[intervals>=1])
                
                # Obtain starting values for parameters from discrete model
                if(is.null(startVals)) {
                  mod_disc = createModel(data, manifestNames, Tpoints, type='discrete')
                  startVals = fit(mod_disc, data, manifestNames, type='discrete', hyperparams=list(missing='fiml'))
                }
                A_disc = startVals[,1:m]
                b_disc = startVals[,m+1]
                predictions_disc = makePredictions(startVals, data, type='discrete')
                errCov_disc = evaluate(predictions_disc, data, manifestNames, getErrCov=TRUE)
                #A_disc.eigen = eigen(A_disc)
                #A_cont = if(hyperparams$approx == 'crude') (A_disc-diag(1,l)) / delta_t else (A_disc.eigen$vectors)%*%log(diag(A_disc.eigen$values)+(matrix(1,l,l)-diag(l)))%*%solve(A_disc.eigen$vectors) / (if(hyperparams$approx == 'unitIntv') 1 else delta_t)
                A_cont = logm(A_disc)/ (if(hyperparams$approx == 'unitIntv') 1 else delta_t)
                b_cont = solve(A_disc-diag(l)) %*% A_cont %*% b_disc
                b_cont[round(b_cont,2) == 0] = 0.01
                A_hash = A_cont %x% diag(l) + diag(l) %x% A_cont
                Q_cont = matrix(solve(solve(A_hash) %*% (expm(A_hash*delta_t) - diag(l*l))) %*% c(errCov_disc), nrow=l, ncol=l)
            
                # Define the model using the obtained starting values
                mod = voelkleModel(n.manifest   = m,                                                                       # total number of indicators for ALL factors. 
                                   n.latent     = l,                                                                       # number of latent variables
                                   PHI1         = cov(data[,paste0(manifestNames, '_T0')], use = "pairwise.complete.obs"), # var/cov matrix of latent variables at first time point
                                   latentM1     = matrix(colMeans(data[,paste0(manifestNames, '_T0')], na.rm=TRUE), nrow=l, ncol=1),   # means of latent variables at first time point
                                   manifestM    = matrix(0, nrow=m, ncol=1, byrow=TRUE),                                   # intercepts of manifest variables (usually fixed to zero)
                                   LAMBDA       = matrix(diag(1, m), nrow=m, ncol=l,byrow=TRUE),                           # factor loading matrix
                                   THETA        = matrix(0,nrow=m, ncol=m),                                                # var/cov matrix of measurement error
                                   DRIFT        = round(A_cont,2),                                                         # crude approximation of drift matrix
                                   CINT         = matrix(round(b_cont,2), ncol=l, nrow=1, byrow=TRUE),                     # continuous time intercepts
                                   Q            = round(Q_cont,2),                                                         # diffusion matrix	
                                   delta_t      = colMeans(data[tail(names(data), Tpoints-1)]),                            # vector of (possibly different) time intervals
                                   data         = data[,1:(l*Tpoints)],
                                   Tpoints      = Tpoints, 
                                   measurements = 1:(Tpoints-1)) 
                
                # Fit the model 
                fit = mxRun(mod)
                summary = summary(fit)$parameters[,c('name', 'matrix', 'Estimate')]
                A_cont = matrix(summary[summary[,'matrix']=='DRIFT', 'Estimate'], ncol=m)
                b_cont = matrix(summary[summary[,'matrix']=='CINT', 'Estimate'], ncol=1)
                estimates = matrix(cbind(A_cont, b_cont), nrow=m, dimnames=list(manifestNames, c(manifestNames, 'intcp')))
                
                return(estimates)
              }
            }
            
            makePredictions = function(estimates, data, type) {
              manifestNames = rownames(estimates)
              Tpoints = ceiling(ncol(data)/(length(manifestNames)+1))
              predictions = data.frame(matrix(nrow=nrow(data), ncol=0))
              
              # Predictions independent of time interval
              if(type=='ctsem' || type=='voelkle') {
                # Derive dicrete drift matrixes for all possible time intervals
                intervals = unique(as.vector(as.matrix(data[,(ncol(data)-Tpoints+2):ncol(data)])))
                estimates_disc = matrix(ncol=length(manifestNames)+1+1, nrow=0)
                for(interval in intervals) {
                  estimate_disc = contToDisc(estimates=estimates, interval=interval)
                  estimates_disc = rbind(estimates_disc, cbind(interval=interval, estimate_disc))
                }
                # Make predictions
                for (T in 1:(Tpoints-1)) {
                  for (var in manifestNames) {
                    rel_cols = cbind(data.matrix(data[,c(paste0('dT',T),paste0(manifestNames, '_T', T-1))]), 1)
                    prediction = c()
                    for (i in 1:nrow(rel_cols)) {
                      val = rel_cols[i,-1]
                      est = estimates_disc[(estimates_disc[,'interval']==rel_cols[i,1]) & (rownames(estimates_disc)==var), -1]
                      prediction = rbind(prediction, val %*% est)
                    }
                    predictions[paste0(var,'_T',T)] = prediction
                  }
                }
              }
              
              # Predictions independent of time interval
              else {
                for (T in 1:(Tpoints-1)) {
                  for (var in manifestNames) {
                    rel_cols = cbind(data.matrix(data[,paste0(manifestNames, '_T', T-1)]), 1)
                    prediction = (rel_cols %*% estimates[var,])
                    predictions[paste0(var,'_T',T)] = prediction
                  }
                }
              }
              
              return(predictions)
            }
            
            evaluate = function(predictions, data, manifestNames, getErrCov=FALSE) {
              
              residuals_all = matrix(ncol=0, nrow=nrow(data))
              for (var in manifestNames) {
                residuals_var = c()
                Tpoints = ceiling(ncol(data)/(length(manifestNames)+1))
                for (T in 1:(Tpoints-1)) {
                  col = paste0(var,'_T',T)
                  residuals = data[col] - predictions[col]
                  residuals = residuals[complete.cases(residuals),]
                  residuals_var = c(residuals_var, residuals)
                }
                residuals_all = cbind(residuals_all, residuals_var)
              }
              
              if(getErrCov) {
                return(var(residuals_all))
              }
              else {
                return(sqrt(colMeans(residuals_all^2)))
              }
            }
            
            tryModel = function(data, manifestNames, type, hyperparams=list(), k=1, valSize=0.2, startVals=NULL) {
              
              m=length(manifestNames)
              
              # Data preparation
              if (type=='ctsem' || type=='voelkle') {
                data = intervalise(data, manifestNames)
              }
              set.seed(123)
              data = data[sample(nrow(data)),]  # Randomly shuffle the data
              set.seed(123)
              if(k>1) { folds = cut(seq(1, nrow(data)), breaks=k, labels=FALSE) } # Create equally size folds 
              
              # Build the model
              Tpoints = ceiling(ncol(data)/(length(manifestNames)+1))
              mod = createModel(data, manifestNames, Tpoints, type, hyperparams)
              
              # Perform CV for parameter estimation
              estimates_sum = matrix(0, nrow=length(manifestNames), ncol=length(manifestNames)+1)
              res = data.frame(fold=integer(), target=character(), rmse_train=double(), rmse_val=double()) 
              for(fold in 1:k){
                
                # Retrieve data for current fold (if no k specified, use valSize argument to split into train and validation set)
                valIndeces = if(k>1) which(folds==fold, arr.ind=TRUE) else 1:floor(nrow(data)*valSize)
                val         = data[valIndeces, ]
                train       = data[-valIndeces, ]
                
                # Retrieve starting values (only for Voelkle model)
                if(!is.null(startVals)) {
                  startVals_fold = as.matrix(startVals[startVals[,'fold']==fold, c(manifestNames, 'intcp')])
                  rownames(startVals_fold) = manifestNames
                }
                else { startVals_fold = NULL }
                
                # Fit the model
                estimates_fold = fit(mod, train, manifestNames, type, hyperparams, startVals_fold)
                estimates_sum = estimates_sum+estimates_fold
                if (type=='ctsem' || type=='voelkle') {
                  intervals = c(as.matrix(data[,(m*Tpoints+1):ncol(data)]))
                  delta_t = mean(intervals[intervals>=1])
                  estimates_disc = contToDisc(estimates=estimates_fold, interval=delta_t) 
                }
                else {
                  estimates_disc = estimates_fold
                }
                
                res = rbind(res, data.frame(fold=fold, target=rownames(estimates_disc), estimates_disc))
              }
              estimates = estimates_sum/k
              
              # Perform CV for model evaluation
              for(fold in 1:k){
                
                # Retrieve data for current fold (if no k specified, use valSize argument to split into train and validation set)
                valIndeces = if(k>1) which(folds==fold, arr.ind=TRUE) else 1:floor(nrow(data)*valSize)
                val         = data[valIndeces, ]
                train       = data[-valIndeces, ]
                
                # Evaluate the model
                pred_train  = makePredictions(estimates, train, type)
                res_train   = evaluate(pred_train, train, manifestNames)
                pred_val  = makePredictions(estimates, val, type)
                res_val    = evaluate(pred_val, val, manifestNames)
                res[((fold-1)*m+1):(fold*(m)), c('rmse_train','rmse_val')] = cbind(res_train, res_val)
              }
              
              return(res)
            }
            
            compareModels = function(models, data, manifestNames, keepOnlyWaves=NULL, exclNANs=FALSE, k=1, split=list(train=0.6, val=0.2, test=0.2)) {
              
              cat('Transformation into wide dataset starts...\n')
              data_wide = createWideDataSet(data, manifestNames, keepOnlyWaves=keepOnlyWaves, exclNANs=exclNANs)
              cat('Transformation into wide dataset is complete.\n')
              
              # Split into train and and test set
              set.seed(123)
              data_wide = data_wide[sample(nrow(data_wide)),]  # Randomly shuffle the data
              testIndeces = 1:floor(nrow(data_wide)*split$test)
              data_test   = data_wide[testIndeces, ]
              data_train  = data_wide[-testIndeces, ]
              
              if(ncol(data_wide) != (length(manifestNames)+1)*2) { models = models[names(models)!='lm'] }
              
              res_combined = data.frame()
              for (i in (1:length(models))) {
                modelID = sort(names(models))[i]
                modelType = if(is.null(models[[modelID]]$type)) modelID else models[[modelID]]$type
                startVals = if(modelType=='voelkle' & nrow(res_combined)>0) res_combined[res_combined[,'model']=='discrete_fiml',] else NULL
                cat('Definition & Training of', modelID, 'model begins...\n')
                res = tryModel(data_wide, manifestNames, modelType, hyperparams=models[[modelID]], k=k, valSize=split$val, startVals=startVals)
                res_combined = rbind(res_combined, cbind(model=modelID, res))
                cat('Evaluation of', modelID, 'model is comple.\n')
              }
              
              return(res_combined)
            }
            
            compareResults = function(results) {
              
              # Aggregate results from CV
              aggregated = aggregate(results[,-(1:3)], by=list(model=results$model, target=results$target), FUN=mean)
              
              # Compute variances across folds
              variances = matrix(ncol=1, nrow=0, dimnames=list(c(), c('var')))
              for (i in 1:nrow(aggregated)) {
                row = aggregated[i,]
                variances = rbind(variances, mean(diag(var(results[results$model==row$model & results$target==row$target, (4:(ncol(row)-2))]))))
              }
              aggregated = cbind(aggregated, variances)
              
              # Round values to 3 decimals
              aggregated[,-(1:2)] = round(aggregated[,-(1:2)], 3)
              
              # Output
              print(aggregated[order(aggregated$target, aggregated$rmse_val),])
            }
            
            models = list(lm                = list(type='lm'),
                          voelkle_crude     = list(type='voelkle',  approx='crude'),
                          voelkle_unitIntv  = list(type='voelkle',  approx='unitIntv'),
                          voelkle_precise   = list(type='voelkle',  approx='precise'),
                          ctsem             = list(type='ctsem',    retryattempts=30),
                          xsimple           = list(type='simple',   missing='fiml'),
                          discrete_listwise = list(type='discrete', missing='listwise'),
                          discrete_fiml     = list(type='discrete', missing='fiml'))
            
            # Train the models
            results1.1 = compareModels(models, data5, manifestNames=c('sat6', 'per1i6'), keepOnlyWaves=c(4,5), exclNANs=TRUE, k=5)
            compareResults(results1.1)
            results2.1 = compareModels(models, data5, manifestNames=c('sat6', 'atts'), keepOnlyWaves=c(1,2), exclNANs=TRUE, k=5)
            compareResults(results2.1)
            results3.1 = compareModels(models, data5, manifestNames=c('sat6', 'atto'), keepOnlyWaves=c(1,3), exclNANs=TRUE, k=5)
            compareResults(results3.1)
            results4.1 = compareModels(models, data5, manifestNames=c('atts', 'atto'), keepOnlyWaves=c(1,2), exclNANs=TRUE, k=5)
            compareResults(results4.1)
            results5.1 = compareModels(models, data5, manifestNames=c('atts', 'per1i6'), keepOnlyWaves=c(4,5), exclNANs=TRUE, k=5)
            compareResults(results5.1)
            results6.1 = compareModels(models, data5, manifestNames=c('atto', 'per1i6'), keepOnlyWaves=c(3,5), exclNANs=TRUE, k=5)
            compareResults(results6.1)
            
            # Train the models
            results1.2 = compareModels(models, data5, manifestNames=c('sat6', 'per1i6'), keepOnlyWaves=c(4,5), exclNANs=FALSE, k=5)
            compareResults(results1.2)
            results2.2 = compareModels(models, data5, manifestNames=c('sat6', 'atts'), keepOnlyWaves=c(1,2), exclNANs=FALSE, k=5)
            compareResults(results2.2)
            results3.2 = compareModels(models, data5, manifestNames=c('sat6', 'atto'), keepOnlyWaves=c(1,3), exclNANs=FALSE, k=5)
            compareResults(results3.2)
            results4.2 = compareModels(models, data5, manifestNames=c('atts', 'atto'), keepOnlyWaves=c(1,2), exclNANs=FALSE, k=5)
            compareResults(results4.2)
            results5.2 = compareModels(models, data5, manifestNames=c('atts', 'per1i6'), keepOnlyWaves=c(4,5), exclNANs=FALSE, k=5)
            compareResults(results5.2)
            results6.2 = compareModels(models, data5, manifestNames=c('atto', 'per1i6'), keepOnlyWaves=c(3,5), exclNANs=FALSE, k=5)
            compareResults(results6.2)
            
            # Train the models
            results1.3 = compareModels(models, data5, manifestNames=c('sat6', 'per1i6'), exclNANs=FALSE, k=5)
            compareResults(results1.3)
            results2.3 = compareModels(models, data5, manifestNames=c('sat6', 'atts'), exclNANs=FALSE, k=5)
            compareResults(results2.3)
            results3.3 = compareModels(models, data5, manifestNames=c('sat6', 'atto'), exclNANs=FALSE, k=5)
            compareResults(results3.3)
            results4.3 = compareModels(models, data5, manifestNames=c('atts', 'atto'), exclNANs=FALSE, k=5)
            compareResults(results4.3)
            results5.3 = compareModels(models, data5, manifestNames=c('atts', 'per1i6'), exclNANs=FALSE, k=5)
            compareResults(results5.3)
            results6.3 = compareModels(models, data5, manifestNames=c('atto', 'per1i6'), exclNANs=FALSE, k=5)
            compareResults(results6.3)


            

               

#--- SECTION 7: FUNCTIONS FOR FLEXIBLE EXPERIMENTS ------------------------------------------------------------------------------------------------------------------------------------------#

prepData <- function(df_long, manifestNames, timeVar="wave", minTpoints=2, exclWaves=NULL, keepOnlyWaves=NULL, exclNANs=FALSE, subj_types=c('anchor', 'partner'), genders=c(1,2)) {
  
  # Keep only relevant subjects
  df_long = df_long[df_long$subj_type %in% subj_types,]
  df_long = df_long[df_long$gen %in% genders,]
  
  # Exclude/Keep specifed waves from the dataset (parameters: exclWaves, keepWaves)
  exclWaves = if(is.null(keepOnlyWaves)) exclWaves else setdiff(1:11, keepOnlyWaves)
  for(exclWave in exclWaves) {
    df_long = df_long[df_long$wave != exclWave,]
  }
  
  # Remove time points with only NANs for at least one variable (anchors and partners treated separately)
  aggr = aggregate(df_long[,manifestNames], by=list(wave=df_long$wave, subj=str_sub(df_long$id,-1)), FUN=mean, na.rm=TRUE) # Mean variable values for each manifest variable per wave and subject (0=anchor; 1=partner) combination
  emptyTpoints = aggr[rowSums(is.na(aggr[,manifestNames]))>0, c('wave')] # Waves without any data for at least one of the selected subj_types
  for(emptyTpoint in emptyTpoints) {
    df_long = df_long[df_long$wave!=emptyTpoint,]
  }
  
  # Keep only subjects with data for at least minTpoints time points (parameter: minTpoints)
  frequencies = data.frame(table(df_long['id']))
  df_long = merge(x=df_long, y=frequencies, by.x="id", by.y='Var1', all.x = TRUE)
  df_long = df_long[df_long[,'Freq']>=minTpoints,]
  
  
  # Transform & Intervalise long to wide data format (Default wide data format)
  df_wide = reshape(df_long[,c('id', 'wave', manifestNames)], idvar="id", timevar="wave", sep=paste0("_T"), direction="wide") # Reshape long to wide format
  rownames(df_wide) = df_wide$id # Set ID column as index
  df_wide$id = NULL # Drop ID column
  print(colnames(df_wide))
  df_wide = df_wide[, order(as.integer(sub(".*_T", "", colnames(df_wide))))] # Reorder time points
  print(colnames(df_wide))
  colnames(df_wide) = paste0(sub("\\_.*", "", colnames(df_wide)), '_T', rep(0:ncol(df_wide), each=length(manifestNames), len=ncol(df_wide))) # Rename columns
  print(colnames(df_wide))
  delta_t = diff(sort(unique(df_long[,'wave'])))
  df_intv = cbind(df_wide, matrix(delta_t, nrow=1, dimnames=list(c(),paste0('dT',(1:length(delta_t)))))) # Add time interval columns
   
  # For the ctsem model, we need a special wide & intervalised dataset
  df_wide_ctsem = ctLongToWide(datalong=df_long, id='id', time=timeVar, manifestNames=manifestNames)
  df_intv_ctsem = intervalise(df_wide_ctsem, manifestNames, individualRelativeTime=TRUE)
  
  # Remove NANs, if specified in 'exclNANs'
  if(exclNANs) {
    df_intv       = df_intv[complete.cases(df_intv), ]
    df_intv_ctsem = df_intv_ctsem[complete.cases(df_intv_ctsem), ]
  }
  
  # Transform back?!
  #Tpoints = (ncol(df_intv_ctsem)+1)/(length(manifestNames)+1)
  #df_long = ctWideToLong(datawide=df_intv, Tpoints=Tpoints, length(manifestNames), manifestNames=manifestNames)
  #df_long = ctDeintervalise(df_long)
  
  # Print dataset characteristics
  Tpoints = (ncol(df_intv_ctsem)+1)/(length(manifestNames)+1)
  intvls = as.matrix(df_intv_ctsem[,(((length(manifestNames)*Tpoints)+1):ncol(df_intv_ctsem))])
  intvls = intvls[intvls!=0.001]
  cat('Number of subjects:           ', nrow(df_intv), '\n')
  cat('Average number of time points:', if(exclNANs) length(delta_t)+1 else round(nrow(df_long)/nrow(df_intv),2), '\n')
  cat('Mean interval:                ', round(mean(intvls),2), '\n')
  cat('Variance of time intervals:   ', round(var(intvls),2), '\n')
  cat('Percentage of NANs:           ', round(sum(is.na(df_long[,manifestNames]))/(nrow(df_long)*length(manifestNames))*100,1), '%\n')
  
  # Wrap everything inside a data package
  #data_pckge = list(wide=df_intv, long=df_long)
  #return(data_pckge)
  
  return(df_intv)
}

build = function(data, manifestNames, Tpoints, type, hyperparams) {
  
  if(type=='simple') {
    mod = ''
    for(var in manifestNames) {
      for(T in 1:(Tpoints-1)) {
        target = paste(var, T, sep='_T')
        inputs = paste0(var, '_T', T-1)
        if(sum(complete.cases(data[,c(target, inputs)])) > 0) { # Only add equations if at least 1 line with data for all involved variables
          regr = paste0(target,' ~ ', paste0(paste0('auto_', var), ' * ', inputs))
          intcp = paste0(target,' ~ intcp_',var,'*1')
          mod = paste(mod, regr, intcp, sep=' \n ' )
        }
      }
    }
    return(mod)
  }
  
  else if(type=='lm' && Tpoints==2) {
    mod = c()
    for(var in manifestNames) {
      mod = append(mod, as.formula(paste0(paste0(var,'_T1'), ' ~ ', paste(lapply(manifestNames, paste0, "_T0"), collapse=' + '))))
    }
    return(mod)
  }
  
  else if(type=='discrete') {
    mod = ''
    for(var in manifestNames) {
      for(T in 1:(Tpoints-1)) {
        target = paste(var, T, sep='_T')
        inputs = paste0(manifestNames, '_T', T-1)
        if(sum(complete.cases(data[,c(target, inputs)])) > 0) { # Only add equations if at least 1 line with data for all involved variables
          regr = paste0(target,' ~ ',paste(paste0(var,'_',manifestNames, '*' ,inputs), collapse=' + '))
          intcp = paste0(target,' ~ intcp_',var,'*1')
          mod = paste(mod, regr, intcp, sep=' \n ' )
        }
      }
    }
    for(i in 1:length(manifestNames)) {
      for(j in 1:length(manifestNames)) {
        if(j>i) {
          for(T in 1:(Tpoints-1)) {
            var_i = paste0(manifestNames[i],'_T',T)
            var_j = paste0(manifestNames[j],'_T',T)
            if(grepl(paste(var_i,'~ '), mod, fixed=TRUE) && grepl(paste(var_j,'~ '), mod, fixed=TRUE)) { # Only add correlations for variables that have been defined as target earlier
              mod = paste(mod, paste(paste(var_i,'~~',var_j), collapse= ' \n '), sep=' \n ')
            }
          }
        }
      }
    }
    return(mod)
  }
  
  else if(type=='stanct' || type=='omx') {
    hyperparams$cint =          if (is.null(hyperparams$cint))             matrix(paste0("cint_",manifestNames),ncol=1) 
                                else if (hyperparams$cint=='free')         matrix(paste0("cint_",manifestNames),ncol=1)
                                else if (hyperparams$cint=='zero')         'auto'
                                else                                       'auto'
    
    hyperparams$manifestmeans = if(is.null(hyperparams$manifestmeans))     matrix(0, nrow=length(manifestNames), ncol=1)
                                else if(hyperparams$manifestmeans=='free') 'auto'  
                                else if(hyperparams$manifestmeans=='zero') matrix(0, nrow=length(manifestNames), ncol=1) 
                                else                                       matrix(0, nrow=length(manifestNames), ncol=1)
    if(is.null(hyperparams$traitvar))  { hyperparams$traitvar=NULL }
    
    mod = ctModel(type=type, 
                  Tpoints=Tpoints,
                  n.manifest=length(manifestNames),
                  manifestNames=manifestNames, 
                  MANIFESTMEANS = hyperparams$manifestmeans,
                  MANIFESTVAR = matrix(0, nrow=length(manifestNames), ncol=length(manifestNames)),
                  CINT = hyperparams$cint,
                  TRAITVAR = hyperparams$traitvar, 
                  LAMBDA=diag(length(manifestNames)))
    
    return(mod)
  }
  
  else { # The definition of the Voelkle model requires a data object. Therefore it is part of the 'fit' step
    return(NULL)
  }
}

fit = function(mod, data, manifestNames, type, hyperparams=list(), startVals=NULL) {
  
  m = length(manifestNames)
  
  if (type=='naive') { 
    estimates = matrix(cbind(diag(m), 0), nrow=m, dimnames=list(manifestNames, c(manifestNames, 'intcp')))
    return(estimates)
  }
  
  if (type=='simple') {
    if (is.null(hyperparams$missing)) { hyperparams$missing = 'listwise' }
    fit = sem(mod, data=data, missing=hyperparams$missing)
    coefs = coef(fit)[!duplicated(names(coef(fit)))]
    estimates = matrix(cbind(diag(coefs[(1:m)*2-1]), coefs[(1:m)*2]), nrow=m, dimnames=list(manifestNames, c(manifestNames, 'intcp')))
    return(estimates)
  }
  
  if (type=='discrete') {
    if (is.null(hyperparams$missing)) { hyperparams$missing = 'listwise' }
    fit = sem(mod, data=data, missing=hyperparams$missing, fixed.x=FALSE)
    coefs = coef(fit)[!duplicated(names(coef(fit)))]
    estimates = matrix(coefs[1:(m*(m+1))], nrow=m, byrow=TRUE, dimnames=list(manifestNames, c(manifestNames, 'intcp')))
    return(estimates)
  }
  
  else if (type=='lm') {
    estimates = matrix(nrow=0, ncol=m+1, dimnames=list(c(), c(manifestNames, 'intcp')))
    for(i in 1:m) {
      fit = lm(mod[[i]], data=data)
      estimates = rbind(estimates, matrix(c(fit$coefficients[c(2:(m+1),1)]), nrow=1, dimnames=list(manifestNames[i], c())))
    }
    return(estimates)
  }
  
  else if (type=='omx') {
    retryattempts = if (is.null(hyperparams$retryattempts)) 24 else hyperparams$retryattempts
    if(is.null(hyperparams$carefulFit)) { hyperparams$carefulFit=TRUE }
    set.seed(42)
    fit = ctFit(ctmodelobj=mod, dat=data, dataform="wide", carefulFit=hyperparams$carefulFit, retryattempts=retryattempts)
    res = summary(fit)
    estimates = matrix(cbind(res$DRIFT, res$CINT, res$MANIFESTMEANS), nrow=m, dimnames=list(manifestNames, c(manifestNames, 'intcp', 'mmeans')))
    return(estimates)
  }
  
  else if (type=='stanct') {
    retryattempts = if (is.null(hyperparams$retryattempts)) 24 else hyperparams$retryattempts
    set.seed(42)
    fit = ctStanFit(datalong=data, ctstanmodel=mod, cores=6, verbose=1)
    print(ctStanContinuousPars(fit))
  }
  
  else if (type=='voelkle') {
    m = length(manifestNames)
    l = m #  number of latentvariables
    Tpoints = (ncol(data)+1)/(length(manifestNames)+1)
    intvls = as.vector(as.matrix(data[1,(m*Tpoints+1):ncol(data)]))
    mean_dt = mean(intvls)
    retryattempts = if (is.null(hyperparams$retryattempts)) 24 else hyperparams$retryattempts
    
    # Obtain starting values for parameters from discrete model
    if(is.null(startVals)) {
      mod_disc = build(data, manifestNames, Tpoints, type='discrete')
      startVals = fit(mod_disc, data, manifestNames, type='discrete', hyperparams=list(missing='fiml'))
    }
    A_disc = startVals[,1:m]
    b_disc = startVals[,m+1]
    predictions_disc = pred(startVals, data, type='discrete')
    errCov_disc = rmse(predictions_disc, data, manifestNames, getErrCov=TRUE)
    A_cont = logm(A_disc)/ (if(hyperparams$approx == 'unitIntv') 1 else mean_dt)
    b_cont = solve(A_disc-diag(l)) %*% A_cont %*% b_disc
    A_hash = A_cont %x% diag(l) + diag(l) %x% A_cont
    Q_cont = matrix(solve(solve(A_hash) %*% (expm(A_hash*mean_dt) - diag(l*l))) %*% c(errCov_disc), nrow=l, ncol=l)
    
    A_cont[A_cont == 0] = 0.01
    b_cont[b_cont == 0] = 0.01
    A_cont[A_cont == 1] = 0.99
    b_cont[b_cont == 1] = 0.99
    
    #print(matrix(colMeans(data[,paste0(manifestNames, '_T0')], na.rm=TRUE), nrow=l, ncol=1))
    #print(matrix(0, nrow=m, ncol=1, byrow=TRUE))
    #print(matrix(diag(1, m), nrow=m, ncol=l,byrow=TRUE))
    ##(A_cont)
    #print(b_cont)
    
    mod = voelkleModel(n.manifest   = m,                                                                       # total number of indicators for ALL factors. 
                       n.latent     = l,                                                                       # number of latent variables
                       PHI1         = cov(data[,paste0(manifestNames, '_T0')], use = "pairwise.complete.obs"), # var/cov matrix of latent variables at first time point
                       latentM1     = matrix(colMeans(data[,paste0(manifestNames, '_T0')], na.rm=TRUE), nrow=l, ncol=1),   # means of latent variables at first time point
                       manifestM    = matrix(0, nrow=m, ncol=1, byrow=TRUE),                                   # intercepts of manifest variables (usually fixed to zero)
                       LAMBDA       = matrix(diag(1, m), nrow=m, ncol=l,byrow=TRUE),                           # factor loading matrix
                       THETA        = matrix(0,nrow=m, ncol=m),                                                # var/cov matrix of measurement error
                       DRIFT        = A_cont,                                                         # crude approximation of drift matrix
                       CINT         = matrix(b_cont, ncol=l, nrow=1, byrow=TRUE),                     # continuous time intercepts
                       Q            = round(Q_cont,4),                                                         # diffusion matrix	
                       delta_t      = intvls,                                                                  # vector of (possibly different) time intervals
                       data         = data[,(1:(m*Tpoints))],
                       Tpoints      = Tpoints, 
                       measurements = 1:(Tpoints-1)) 
    
    # Fit the model 
    fit = mxTryHard(mod, extraTries = 10)
    summary = summary(fit)$parameters[,c('name', 'matrix', 'Estimate')]
    A_cont = matrix(summary[summary[,'matrix']=='DRIFT', 'Estimate'], ncol=m)
    b_cont = matrix(summary[summary[,'matrix']=='CINT', 'Estimate'], ncol=1)
    estimates = matrix(cbind(A_cont, b_cont), nrow=m, dimnames=list(manifestNames, c(manifestNames, 'intcp')))

    #print(A_cont)
    #print(b_cont)
    #print(summary)   
    
    return(estimates)
  }
}

pred = function(estimates, data, type) {
  manifestNames = rownames(estimates)
  Tpoints = ceiling(ncol(data)/(length(manifestNames)+1))
  predictions = data.frame(matrix(nrow=nrow(data), ncol=0))
  # Predictions dependent on time interval
  if(type=='stanct' || type=='omx' || type=='voelkle') {
    # Derive dicrete drift matrixes for all possible time intervals
    intervals = unique(as.vector(as.matrix(data[,(ncol(data)-Tpoints+2):ncol(data)])))
    estimates_disc = matrix(ncol=length(manifestNames)+1+1, nrow=0)
    for(interval in intervals) {
      estimate_disc = contToDisc(estimates=estimates, interval=interval)
      estimates_disc = rbind(estimates_disc, cbind(interval=interval, estimate_disc))
    }
    # Make predictions
    for (T in 1:(Tpoints-1)) {
      for (var in manifestNames) {
        rel_cols = cbind(data.matrix(data[,c(paste0('dT',T),paste0(manifestNames, '_T', T-1))]), 1)
        prediction = c()
        for (i in 1:nrow(rel_cols)) {
          val = rel_cols[i,-1]
          est = estimates_disc[(estimates_disc[,'interval']==rel_cols[i,1]) & (rownames(estimates_disc)==var), -1]
          prediction = rbind(prediction, val %*% est)
        }
        predictions[paste0(var,'_T',T)] = prediction
      }
    }
  }
  
  # Predictions independent of time interval
  else {
    for (T in 1:(Tpoints-1)) {
      for (var in manifestNames) {
        rel_cols = cbind(data.matrix(data[,paste0(manifestNames, '_T', T-1)]), 1)
        prediction = (rel_cols %*% estimates[var,])
        predictions[paste0(var,'_T',T)] = prediction
      }
    }
  }
  
  return(predictions)
}

rmse = function(predictions, data, manifestNames, getErrCov=FALSE) {
  
  residuals_all = matrix(ncol=0, nrow=nrow(data))
  for (var in manifestNames) {
    residuals_var = c()
    Tpoints = ceiling(ncol(data)/(length(manifestNames)+1))
    for (T in 1:(Tpoints-1)) {
      col = paste0(var,'_T',T)
      residuals = data[col] - predictions[col]
      residuals = residuals[complete.cases(residuals),]
      residuals_var = c(residuals_var, residuals)
    }
    residuals_all = cbind(residuals_all, residuals_var)
  }
  
  if(getErrCov) {
    return(var(residuals_all))
  }
  else {
    return(sqrt(colMeans(residuals_all^2)))
  }
}

evalMod = function(data, manifestNames, type, hyperparams=list(), k=1, valSize=0.2, startVals=NULL) {
  
  m=length(manifestNames)
  
  set.seed(42)
  if(k>1) { folds = cut(seq(1, nrow(data)), breaks=k, labels=FALSE) } # Create equally size folds
  
  # Build the model
  Tpoints = ceiling(ncol(data)/(length(manifestNames)+1))
  mod = build(data, manifestNames, Tpoints, type, hyperparams)
  
  # Perform CV for parameter estimation
  estimates_sum = matrix(0, nrow=length(manifestNames), ncol=length(manifestNames)+1)
  res = data.frame(paramType=character(), fold=integer(), target=character(), sec=double(), rmse_train=double(), rmse_val=double()) 
  
  for(fold in 1:k){
    
    # Split into train & validation set
    if(k==1) { # i.e. if no CV enabled
      set.seed(42)
      data <- data[sample(nrow(data)), ] # Shuffle dataset
      valIndeces = 1:floor(nrow(data)*valSize)
    }
    else { # i.e. if CV enabled
      valIndeces = which(folds==fold, arr.ind=TRUE)
    }
    train = data[-valIndeces, ]
    val   = data[valIndeces, ]

    # Retrieve starting values (only for Voelkle model)
    if(!is.null(startVals)) {
      startVals_fold = as.matrix(startVals[startVals[,'fold']==fold, c(manifestNames, 'intcp')])
      rownames(startVals_fold) = manifestNames
    }
    else { startVals_fold = NULL }
    
    # Fit the model
    start_time = Sys.time()
    estimates = fit(mod, train, manifestNames, type, hyperparams, startVals_fold)
    end_time = Sys.time()
    
    # Evaluate model on training and test data
    pred_train = pred(estimates, train, type)
    rmse_train = rmse(pred_train, train, manifestNames)
    pred_val = pred(estimates, val, type)
    rmse_val = rmse(pred_val, val, manifestNames)
    
    # Save results (For continuous models, save continuous and discrete parameters in separate rows)
    if (type=='omx' || type=='stanct' || type=='voelkle') {
      res = rbind.fill(res, data.frame(paramType='cont', fold=fold, target=rownames(estimates), estimates, sec=as.numeric(end_time-start_time), rmse_train=rmse_train))
      intervals = c(as.matrix(data[,(m*Tpoints+1):ncol(data)]))
      delta_t = mean(intervals[intervals>=1])
      estimates_disc = contToDisc(estimates=estimates, interval=delta_t)
      res = rbind.fill(res, data.frame(paramType='disc', fold=fold, target=rownames(estimates_disc), estimates_disc, sec=as.numeric(end_time-start_time), rmse_train=rmse_train, rmse_val=rmse_val))
    }
    else {
      res = rbind.fill(res, data.frame(paramType='disc', fold=fold, target=rownames(estimates), estimates, sec=as.numeric(end_time-start_time), rmse_train=rmse_train, rmse_val=rmse_val))
    }
  }

  return(res)
}

modSel = function(models, dat, manifestNames, mode='normal', k=1, split=list(train=0.6, val=0.2, test=0.2)) {
  
  # Deactivate linear model if Tpoints > 2
  if(ncol(dat) != ((length(manifestNames)+1)*2)-1) { models = models[names(models)!='lm'] }
  
  res_combined = data.frame()
  for (i in (1:length(models))) {
    modelID = sort(names(models))[i]
    modelType = if(is.null(models[[modelID]]$type)) modelID else models[[modelID]]$type
    startVals = if(modelType=='voelkle' & nrow(res_combined)>0) res_combined[res_combined[,'model']=='discrete_fiml',] else NULL
    cat('Definition & Training of', modelID, 'model begins...\n')
    res = evalMod(dat, manifestNames, modelType, hyperparams=models[[modelID]], k=k, valSize=split$val, startVals=startVals)
    res_combined = rbind.fill(res_combined, cbind(model=modelID, res))
    cat('Evaluation of', modelID, 'model is comple.\n')
  }
  
  printRes(res_combined)
  return(res_combined)
}

trainTestSplit = function(dat, subset='train', fold=NULL, manifestNames, exclWaves=NULL, keepOnlyWaves=NULL, exclNANs=FALSE, k=1, split=list(train=0.6, val=0.2, test=0.2), subj_types=c('anchor', 'partner'), genders=c(1,2)) {
  
  #Preprocess data
  dat = prepData(dat, manifestNames,  timeVar='wave', minTpoints=2, exclWaves=exclWaves, keepOnlyWaves=keepOnlyWaves, exclNANs=exclNANs, subj_types=subj_types, genders=genders)
  
  # Split into train and and test set
  set.seed(42)
  dat <- dat[sample(nrow(dat)), ] # Shuffle dataset
  testIDs = seq(1,nrow(dat)*split$test)
  test = dat[testIDs,]
  trainVal = dat[-testIDs,]
  
  if(is.null(fold) == FALSE) {
    # Split data into k folds
    set.seed(42)
    if(k>1) { folds = cut(seq(1, nrow(trainVal)), breaks=k, labels=FALSE) } # Create equally size folds
    
    # Retrieve data for requested fold (if no k specified, use valSize argument to split into train and validation set)
    valIndeces = if(k>1) which(folds==fold, arr.ind=TRUE) else 1:floor(nrow(trainVal)*split$val)
    val    = trainVal[valIndeces,]
    train  = trainVal[-valIndeces, ]
  }
  
  return(get(subset))
}

printRes = function(results, paramType='disc', sortby='rmse_val') {
  
  results = results[results$paramType==paramType,]
  
  # Aggregate results from CV
  aggregated = aggregate(results[,-(1:4)], by=list(model=results$model, target=results$target), FUN=mean)
  
  # Compute variances across folds
  variances = matrix(ncol=1, nrow=0, dimnames=list(c(), c('var')))
  for (i in 1:nrow(aggregated)) {
    row = aggregated[i,]
    variances = rbind(variances, mean(diag(var(results[results$model==row$model & results$target==row$target, (8:ncol(results))])), na.rm=TRUE))
  }
  aggregated = cbind(aggregated, variances)
  
  # Round values to 3 decimals
  aggregated[,-c(1,2,4,5,ncol(aggregated))] = round(aggregated[,-c(1,2,4,5,ncol(aggregated))], 3)
  aggregated[,c(4,5,ncol(aggregated))]  = round(aggregated[,c(4,5,ncol(aggregated))], 4)
  
  # Output
  options(scipen = 50)
  print(aggregated[order(aggregated$target, aggregated[sortby]),])
}

writeRes = function(resNames=NULL) {
  if(is.null(resNames)) { resNames=ls(envir=.GlobalEnv, pat='res') }
  for(resName in resNames) {
    write.csv(get(resName), paste0('../../results/', resName, '.csv'))
  }
}

readRes = function(resNames=NULL) {
  if(is.null(resNames)) { resNames=sub('.csv', '', list.files(path='../../results/', pattern='res')) }
  for(resName in resNames) {
    assign(resName, read.csv(paste0('../../results/', resName, '.csv')), envir=.GlobalEnv)
  }
}



# Model comparison for Tpoints=2

models = list(naive             = list(type='naive'),
              lm                = list(type='lm'),
              voelkle_precise   = list(type='voelkle',  approx='precise'),
              voelkle_crude     = list(type='voelkle',  approx='crude'),
              voelkle_unit      = list(type='voelkle',  approx='unit'),
              ctsem_cint        = list(type='omx', cint='free', manifestmeans='zero', traitvar=NULL),
              simple            = list(type='simple',   missing='fiml'),
              discrete_listwise = list(type='discrete', missing='listwise'),
              discrete_fiml     = list(type='discrete', missing='fiml'))

manifestNames=c('sat6', 'per1i6')
res1a = modSel(models, dat=prepData(data7, manifestNames, keepOnlyWaves=c(4,5), exclNANs=TRUE), manifestNames=manifestNames, k=5)
res1b = modSel(models, dat=prepData(data7, manifestNames), manifestNames=manifestNames, k=5)

manifestNames=c('sat6', 'attAnx')
res2a = modSel(models, dat=prepData(data7, manifestNames, keepOnlyWaves=c(1,3), exclNANs=TRUE), manifestNames=manifestNames, k=5)
res2b = modSel(models, dat=prepData(data7, manifestNames), manifestNames=manifestNames, k=5)

manifestNames=c('sat6', 'attAvd')
res3a = modSel(models, dat=prepData(data7, manifestNames, keepOnlyWaves=c(1,3), exclNANs=TRUE), manifestNames=manifestNames, k=5)
res3b = modSel(models, dat=prepData(data7, manifestNames), manifestNames=manifestNames, k=5)

manifestNames=c('attAnx', 'attAvd')
res4a = modSel(models, dat=prepData(data7, manifestNames, keepOnlyWaves=c(1,3), exclNANs=TRUE), manifestNames=manifestNames, k=5)
res4b = modSel(models, dat=prepData(data7, manifestNames), manifestNames=manifestNames, k=5)

manifestNames=c('attAnx', 'per1i6')
res5a = modSel(models, dat=prepData(data7, manifestNames, keepOnlyWaves=c(1,5), exclNANs=TRUE), manifestNames=manifestNames, k=5)
res5b = modSel(models, dat=prepData(data7, manifestNames), manifestNames=manifestNames, k=5)

manifestNames=c('attAvd', 'per1i6')
res6a = modSel(models, dat=prepData(data7, manifestNames, keepOnlyWaves=c(1,5), exclNANs=TRUE), manifestNames=manifestNames, k=5)
res6b = modSel(models, dat=prepData(data7, manifestNames), manifestNames=manifestNames, k=5)

writeRes()






# Model comparison for complete dataset

manifestNames=c('sat6', 'per1i6', 'attAnx', 'attAvd', 'pattAnx', 'pattAvd')

# Split data into trainVal and test set while assuring a 80/20 split for both genders
set.seed(42)
dat7 = data7[sample(nrow(data7)), ] # Shuffle dataset
dat7.m          = prepData(data7, manifestNames, genders=1)
dat7.m.trainVal = dat7.m[seq(1,nrow(dat7.m)*0.8),]
dat7.m.test     = dat7.m[-seq(1,nrow(dat7.m)*0.8),]
dat7.f          = prepData(data7, manifestNames, genders=2)
dat7.f.trainVal = dat7.f[seq(1,nrow(dat7.f)*0.8),]
dat7.f.test     = dat7.f[-seq(1,nrow(dat7.f)*0.8),]
dat7.trainVal = rbind(dat7.m.trainVal, dat7.f.trainVal)
dat7.test = rbind(dat7.m.test, dat7.f.test)

# Compare models using trainVal set
models = list(naive             = list(type='naive'),
              lm                = list(type='lm'),
              voelkle_precise   = list(type='voelkle',  approx='precise'),
              ctsem_cint        = list(type='omx', cint='free', manifestmeans='zero', traitvar=NULL),
              simple            = list(type='simple',   missing='fiml'),
              discrete_fiml     = list(type='discrete', missing='fiml'))
res7.sel = modSel(models, dat7.trainVal, manifestNames, k=5)
res7.m.sel = modSel(models, dat7.m.trainVal, manifestNames, k=5)
res7.f.sel = modSel(models, dat7.trainVal, manifestNames, k=5)

# Final fit of best models using entire data set
best_models = models[names(models)=='voelkle_precise' | names(models)=='discrete_fiml']
res7.fin = modSel(best_models, dat7.test, manifestNames)
res7.m.fin = modSel(best_models, dat7.m.test, manifestNames)
res7.f.fin = modSel(best_models, dat7.f.test, manifestNames)

writeRes()






manifestNames=c('attAnx', 'per1i6')
dat = prepData(data7, manifestNames)
mod = ctModel(type='omx', Tpoints=5, n.manifest=length(manifestNames), 
              manifestNames=manifestNames, 
              MANIFESTMEANS = matrix(0, nrow=length(manifestNames), ncol=1),
              MANIFESTVAR = matrix(0, nrow=length(manifestNames), ncol=length(manifestNames)),
              CINT = matrix(paste0("cint_",manifestNames),ncol=1),
              TRAITVAR = NULL, 
              LAMBDA=diag(length(manifestNames)))
set.seed(42)
fit = ctFit(ctmodelobj=mod, dat=dat, dataform="wide")
res = summary(fit)
res$CINT
res$DRIFT

mean(dat$per1i6_T0, na.rm=TRUE)
mean(dat$per1i6_T1, na.rm=TRUE)
mean(dat$per1i6_T2, na.rm=TRUE)
mean(dat$per1i6_T3, na.rm=TRUE)
mean(dat$per1i6_T4, na.rm=TRUE)
mean(dat$per1i6_T5, na.rm=TRUE)

expm(res$DRIFT*3)
expm(res$DRIFT*4)
expm(res$DRIFT*5)
expm(res$DRIFT*6)
expm(res$DRIFT*7)
expm(res$DRIFT*8)
expm(res$DRIFT*9)
expm(res$DRIFT*10)
expm(res$DRIFT*11)
expm(res$DRIFT*12)
expm(res$DRIFT*13)
expm(res$DRIFT*14)
expm(res$DRIFT*15)


start = c(-1,1)
for (i in c(0:15)) {
  A = expm(res$DRIFT*i)
  b = solve(res$DRIFT) %*% (A-diag(1,2)) %*% res$CINT
  print(A%*%start+b) 
}

start = c(-1,1)
for (i in c(0:15)) {
  A = expm(res$DRIFT*i)
  b = solve(res$DRIFT) %*% (A-diag(1,2)) %*% res$CINT
  print(b) 
}


eigen(res$DRIFT)
