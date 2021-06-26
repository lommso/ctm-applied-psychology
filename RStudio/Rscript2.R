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
data.pairfam = read.csv(file = '../../data/samples/data8.csv')





#--- MANUAL IMPLEMENTATION OF CONTINUOUS TIME MODEL --------------------------------------------------------------------------------------------------------------------#

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
  
  return(model2)
}





#--- MODEL EVALUATION PIPELINE --------------------------------------------------------------------------------------------------------------------#

prepData <- function(df_long, manifestNames, timeVar="wave", minTpoints=2, exclWaves=NULL, keepOnlyWaves=NULL, exclNANs=FALSE, subj_types=c('anchor', 'partner'), genders=c(1,2), N_samples=NULL) {
  
  # Keep only relevant subjects (parameter: subj_type)
  df_long = df_long[df_long$subj_type %in% subj_types,]
  df_long = df_long[df_long$sex %in% genders,]
  
  # Exclude/Keep specifed waves from the dataset (parameters: exclWaves, keepWaves)
  exclWaves = if(is.null(keepOnlyWaves)) exclWaves else setdiff(1:11, keepOnlyWaves)
  for(exclWave in exclWaves) {
    df_long = df_long[df_long$wave != exclWave,]
  }
  
  # Remove time points with only NANs for at least one variable (anchors and partners treated separately)
  aggr = aggregate(df_long[,manifestNames], by=list(wave=df_long$wave, subj=df_long$subj_type), FUN=mean, na.rm=TRUE) # Mean variable values for each manifest variable per wave and subject (0=anchor; 1=partner) combination
  emptyTpoints = aggr[rowSums(is.na(aggr[,manifestNames]))>0, c('wave')] # Waves without any data for at least one of the selected subj_type
  for(emptyTpoint in emptyTpoints) {
    df_long = df_long[df_long$wave!=emptyTpoint,]
  }
  
  # Remove rows with NANs, if specified (parameter: 'exclNANs')
  if(exclNANs) {
    df_long = df_long[complete.cases(df_long[,manifestNames]), ]
  }
  
  # Keep only subjects with data for at least minTpoints time points (parameter: minTpoints)
  frequencies = data.frame(table(df_long['id']))
  df_long = merge(x=df_long, y=frequencies, by.x="id", by.y='Var1', all.x = TRUE)
  df_long = df_long[df_long[,'Freq']>=minTpoints,]
  
  # Reduce sample size if required (parameter: 'sample_size')
  ids = unique(df_long$id)
  N_samples = if(is.null(N_samples) || N_samples>length(ids)) length(ids) else N_samples
  ids = ids[sample(N_samples)] # Generate random sample of IDs
  df_long = merge(x=df_long, y=ids, by.x='id', by.y=1, all.y=TRUE)
  
  # Standardize variables
  df_long[, manifestNames] = scale(df_long[, manifestNames])
  
  # Transform & Intervalise long to wide data format
  df_wide = reshape(df_long[,c('id', 'wave', manifestNames)], idvar="id", timevar="wave", sep=paste0("_T"), direction="wide") # Reshape long to wide format
  rownames(df_wide) = df_wide$id # Set ID column as index
  df_wide$id = NULL # Drop ID column
  df_wide = df_wide[, order(as.integer(sub(".*_T", "", colnames(df_wide))))] # Reorder time points
  colnames(df_wide) = paste0(sub("\\_.*", "", colnames(df_wide)), '_T', rep(0:ncol(df_wide), each=length(manifestNames), len=ncol(df_wide))) # Rename columns
  delta_t = diff(sort(unique(df_long[,'wave'])))
  df_intv = cbind(df_wide, matrix(delta_t, nrow=1, dimnames=list(c(),paste0('dT',(1:length(delta_t)))))) # Add time interval columns
  
  # Print dataset characteristics
  cat('Number of subjects:           ', nrow(df_intv), '\n')
  cat('Average number of time points:', round(nrow(df_long)/nrow(df_intv),2), '\n')
  cat('Mean interval:                ', round(mean(delta_t),2), '\n')
  cat('Variance of time intervals:   ', round(var(delta_t),2), '\n')
  cat('Percentage of NANs:           ', round(sum(is.na(df_long[,manifestNames]))/(nrow(df_long)*length(manifestNames))*100,1), '%\n')
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

fit = function(mod, data, manifestNames, type, hyperparams=list()) {
  
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
    retryattempts = if (is.null(hyperparams$retryattempts)) 10 else hyperparams$retryattempts
    if(is.null(hyperparams$carefulFit)) { hyperparams$carefulFit=TRUE }
    set.seed(42)
    fit = ctFit(ctmodelobj=mod, dat=data, dataform="wide", carefulFit=hyperparams$carefulFit, retryattempts=retryattempts)
    res = summary(fit)
    estimates = matrix(cbind(res$DRIFT, res$CINT, res$MANIFESTMEANS), nrow=m, dimnames=list(manifestNames, c(manifestNames, 'intcp', 'mmeans')))
    return(estimates)
  }
  
  else if (type=='stanct') {
    set.seed(42)
    fit = ctStanFit(datalong=data, ctstanmodel=mod, cores=6, verbose=1)
  }
  
  else if (type=='voelkle') {
    m = length(manifestNames)
    l = m #  number of latentvariables
    Tpoints = (ncol(data)+1)/(length(manifestNames)+1)
    intvls = as.vector(as.matrix(data[1,(m*Tpoints+1):ncol(data)]))
    mean_dt = mean(intvls)
    retryattempts = if (is.null(hyperparams$retryattempts)) 10 else hyperparams$retryattempts
    
    # Obtain starting values for parameters from discrete model
    mod_disc = build(data, manifestNames, Tpoints, type='discrete')
    startVals = fit(mod_disc, data, manifestNames, type='discrete', hyperparams=list(missing='fiml'))
    
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
    fit_voelkle = mxTryHard(mod, extraTries = retryattempts, scale=0.1)
    summary = summary(fit_voelkle)$parameters[,c('name', 'matrix', 'Estimate')]
    A_cont = matrix(summary[summary[,'matrix']=='DRIFT', 'Estimate'], ncol=m)
    b_cont = matrix(summary[summary[,'matrix']=='CINT', 'Estimate'], ncol=1)
    estimates = matrix(cbind(A_cont, b_cont), nrow=m, dimnames=list(manifestNames, c(manifestNames, 'intcp')))
    
    return(estimates)
  }
}

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

evaluate = function(mod, dat.train, dat.val, manifestNames, mod_type, hyperparams=list()) {
  
  # Fit the model
  start_time = Sys.time()
  estimates = fit(mod, dat.train, manifestNames, mod_type, hyperparams)
  end_time = Sys.time()
  secs = difftime(end_time, start_time, units = "secs")[[1]]
  
  # Evaluate model on training and validation data
  pred_train = pred(estimates, dat.train, mod_type)
  rmse_train = rmse(pred_train, dat.train, manifestNames)
  if(!is.null(dat.val)) {
    pred_val = pred(estimates, dat.val, mod_type)
    rmse_val = rmse(pred_val, dat.val, manifestNames)
  }
  else { rmse_val = NA }
  
  # Wrap results for continuous time models
  res = data.frame()
  if (mod_type=='omx' || mod_type=='stanct' || mod_type=='voelkle') {
    res = rbind.fill(res, data.frame(paramType='cont', target=rownames(estimates), 
                                     estimates, secs=secs, rmse_train=rmse_train, rmse_val=rmse_val))
    Tpoints = (ncol(dat.train)+1)/(length(manifestNames)+1)
    intervals = c(as.matrix(dat.train[,(length(manifestNames)*Tpoints+1):ncol(dat.train)]))
    delta_t = mean(intervals[intervals>=1])
    estimates_disc = contToDisc(estimates=estimates, interval=delta_t)
    res = rbind.fill(res, data.frame(paramType='disc', target=rownames(estimates_disc), 
                                     secs=secs, estimates_disc, rmse_train=rmse_train, rmse_val=rmse_val))
  }
  
  # Wrap results for discrete time models
  else { 
    res = rbind.fill(res, data.frame(paramType='disc', target=rownames(estimates), 
                                     estimates, secs=secs, rmse_train=rmse_train, rmse_val=rmse_val))
  }

  return(res)
}

cv = function(mod, dat, manifestNames, mod_type, hyperparams=list(), k=NULL) {
  
  set.seed(42)
  folds = cut(seq(1, nrow(dat)), breaks=k, labels=FALSE) # Create equally size folds
  
  res = data.frame() 
  
  for(fold in 1:k){
    
    # Separate train and validation set
    set.seed(42)
    sample = which(folds!=fold, arr.ind=TRUE)
    dat.train = dat[sample, ]
    dat.val   = dat[-sample, ]
    
    # Evaluate model
    res_fold = evaluate(mod, dat.train, dat.val, manifestNames, mod_type, hyperparams=hyperparams)
    res_fold['fold'] = fold
    res = rbind.fill(res, res_fold)
    print(res_fold)
  }
  
  return(res)
}

modEval = function(dat.trainVal, manifestNames, mod_type, hyperparams=list(), k=NULL, val_size=NULL, dat.test=NULL, final_fit=FALSE) {
  
  # Build model
  Tpoints = ceiling(ncol(dat.trainVal)/(length(manifestNames)+1))
  mod = build(dat.trainVal, manifestNames, Tpoints, mod_type, hyperparams)
  
  res = data.frame() 
  
  # Cross validation, if requested through parameter k
  if(!is.null(k)) { 
    res_cv = cv(mod, dat.trainVal, manifestNames, mod_type, hyperparams, k=k)
    res_cv['data_train'] = 'train'
    res = rbind.fill(res, res_cv)
  }
  
  # Ordinary performance evaluation, if no cross validation, but specified val_size
  else if(!is.null(val_size)) {
    dat.subsets = trainTestSplit(dat.trainVal, val_size)
    dat.train = dat.subsets['train']
    dat.val = dat.subsets['test']
    
    res_val = evaluate(mod, dat.train, dat.val, manifestNames, mod_type, hyperparams)
    res_cv['data_train'] = 'train'
    res = rbind.fill(res, res_test)
  }
  
  # Evaluate performance on test data using training and validation data combined
  else if(!is.null(dat.test)) {
    res_test = evaluate(mod, dat.trainVal, dat.test, manifestNames, mod_type, hyperparams)
    res_test['data_train'] = 'trainVal'
    res = rbind.fill(res, res_test)
  }
  
  # Final fit with complete dataset
  if(final_fit || (is.null(k) & is.null(val_size) & is.null(dat.test))) {
    res_final = evaluate(mod, dat.trainVal, NULL, manifestNames, mod_type, hyperparams)
    res_final['data_train'] = 'trainValTest'
    res = rbind.fill(res, res_final)
  }
  
  return(res)
}

modComp = function(models, dat.trainVal, manifestNames, mode='normal', k=NULL, val_size=NULL, data.test=NULL) {
  
  # Deactivate linear model if Tpoints > 2
  if(ncol(dat.trainVal) != ((length(manifestNames)+1)*2)-1) { models = models[names(models)!='lm'] }
  
  res = data.frame()
  for (i in (1:length(models))) {
    mod_id = sort(names(models))[i]
    mod_type = if(is.null(models[[mod_id]]$type)) mod_id else models[[mod_id]]$type
    cat('Model ', mod_id, 'build & evaluation begins\n')
    res_mod = modEval(dat.trainVal, manifestNames, mod_type, hyperparams=models[[mod_id]], k=k, val_size=val_size)
    res = rbind.fill(res, cbind(model=mod_id, res_mod))
    cat('Model ', mod_id, 'build & evaluation complete\n')
  }
  
  printRes(res)
  return(res)
}

trainTestSplit = function(dat, test_size) {
  
  set.seed(42)
  sample = sample.int(n=nrow(dat), size=floor((1-test_size)*nrow(dat)))
  dat.train = dat[sample, ]
  dat.test = dat[-sample, ]
  
  return(list(train=dat.train, test=dat.test))
}

printRes = function(res, paramType='disc', sortby='rmse_val') {

  res = res[res$paramType==paramType,]
  
  # Aggregate results from CV
  non_numeric_cols = c("model", "data_train", "paramType", "fold", "target")
  res.agg = aggregate(res[, !names(res) %in% non_numeric_cols], 
                      by=list(model=res$model, data_train=res$data_train, target=res$target, paramType=res$paramType), FUN=mean)
  
  # Round values to 3 decimals
  res.agg[, !names(res.agg) %in% non_numeric_cols] = round(res.agg[, !names(res.agg) %in% non_numeric_cols], 4)
  
  # Output
  options(scipen = 50)
  print(res.agg[order(res.agg$target, res.agg[sortby]),])
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





#--- PART 1: PREDICTION ACCURACY ACROSS MODELS --------------------------------------------------------------------------------------------------------------------#

models = list(naive             = list(type='naive'),
              lm                = list(type='lm'),
              voelkle           = list(type='voelkle',  approx='precise'),
              ctsem             = list(type='omx', cint='free', manifestmeans='zero', traitvar=NULL),
              simple            = list(type='simple',   missing='fiml'),
              discrete_listwise = list(type='discrete', missing='listwise'),
              discrete_fiml     = list(type='discrete', missing='fiml'))

manifestNames=c('sat6', 'attAnx')
res1a = modComp(models, dat.trainVal=prepData(data.pairfam, manifestNames, keepOnlyWaves=c(1,3), exclNANs=TRUE), manifestNames=manifestNames, k=5)
res1b = modComp(models, dat.trainVal=prepData(data.pairfam, manifestNames, keepOnlyWaves=c(1,5,7,9,11)), manifestNames=manifestNames, k=5)

manifestNames=c('sat6', 'attAvd')
res2a = modComp(models, dat.trainVal=prepData(data.pairfam, manifestNames, keepOnlyWaves=c(1,3), exclNANs=TRUE), manifestNames=manifestNames, k=5)
res2b = modComp(models, dat.trainVal=prepData(data.pairfam, manifestNames, keepOnlyWaves=c(1,5,7,9,11)), manifestNames=manifestNames, k=5)

manifestNames=c('attAnx', 'attAvd')
res3a = modComp(models, dat.trainVal=prepData(data.pairfam, manifestNames, keepOnlyWaves=c(1,3), exclNANs=TRUE), manifestNames=manifestNames, k=5)
res3b = modComp(models, dat.trainVal=prepData(data.pairfam, manifestNames, keepOnlyWaves=c(1,5,7,9,11)), manifestNames=manifestNames, k=5)


manifestNames=c('sat6', 'per1i6', 'attAnx', 'attAvd')
res4 = modComp(models, dat.trainVal=prepData(data.pairfam, manifestNames, keepOnlyWaves=c(1,5,7,9,11)), manifestNames=manifestNames, k=5)

writeRes()
readRes()
printRes(res3a)



#--- PART 2: COMPUTATIONAL COMPLEXITY OF THE MODELS --------------------------------------------------------------------------------------------------------------------#

performance_test = function(hps, models, data) {
  res = data.frame()
  for (hp in hps) {
    for (N in hp$N) {
      dat = prepData(data, hp$manifestNames, keepOnlyWaves=hp$waves, exclNANs=TRUE, N_samples=N, minTpoints=hp$minTpoints)
      Tpoints = ceiling(ncol(dat)/(length(hp$manifestNames)+1))
      if(is.null(hp$duplicates)) {hp$duplicates = 1}
      for (dups in hp$duplicates) {
        dat.d = do.call("rbind", replicate(dups, dat, simplify = FALSE)) # Appends dataframe multiple times to itself
        for (model in names(models)) {
          start_time = Sys.time()
          for (i in c(1:hp$iterations)) {
            mod = build(dat.d, hp$manifestNames, Tpoints, models[[model]]['type'], models[[model]])
            fit(mod, dat.d, hp$manifestNames, models[[model]]['type'], models[[model]])
          }
          end_time = Sys.time()
          res = rbind(res, list(manifests=paste(unlist(hp$manifestNames), collapse=', '), 
                                N=N, duplicates=dups, waves=paste(unlist(hp$waves), collapse=', '), 
                                model=model, sec=difftime(end_time, start_time, units = "secs")[[1]]/hp$iterations))
        }
      }
    }
  }
  return(res)
}

# Experiment 1: Varying sample size for a simple model and for a complex model
models = list(discrete = list(type='discrete', missing='fiml'),
              voelkle  = list(type='voelkle',  approx='precise', retryattempts=0),
              ctsem    = list(type='omx', cint='free', manifestmeans='zero', traitvar=NULL, retryattempts=0))
hyperparams = list(a = list(manifestNames=c('sat6', 'attAnx'), 
                       waves=c(1,3),
                       N = c(5032),
                       duplicates = c(1,5,10,20,30,50),
                       minTpoints=2,
                       iterations = 25))
res5 = performance_test(hyperparams, models, data.pairfam)
res5

# Experiment 2: Varying number of time points
models = list(discrete = list(type='discrete', missing='fiml'),
              voelkle  = list(type='voelkle',  approx='precise', retryattempts=10),
              ctsem    = list(type='omx', cint='free', manifestmeans='zero', traitvar=NULL, retryattempts=10))
hyperparams = list(a = list(manifestNames=c('sat6', 'attAnx'), 
                            waves=c(1,3),
                            N = c(1000),
                            minTpoints=2,
                            iterations = 50),
                   b = list(manifestNames=c('sat6', 'attAnx'), 
                            waves=c(1,3,5),
                            N = c(1000),
                            minTpoints=3,
                            iterations = 50),
                   c = list(manifestNames=c('sat6', 'attAnx'), 
                            waves=c(1,3,5,7),
                            N = c(1000),
                            minTpoints=4,
                            iterations = 50),
                   d = list(manifestNames=c('sat6', 'attAnx'), 
                            waves=c(1,3,5,7,9),
                            N = c(1000),
                            minTpoints=5,
                            iterations = 50),
                   e = list(manifestNames=c('sat6', 'attAnx'), 
                            waves=c(1,3,5,7,9,11),
                            N = c(1000),
                            minTpoints=6,
                            iterations = 50))
res6 = performance_test(hyperparams, models, data.pairfam)
res6

# Experiment 3: Varying number of parameters
models = list(discrete = list(type='discrete', missing='fiml'),
              voelkle  = list(type='voelkle',  approx='precise', retryattempts=10),
              ctsem    = list(type='omx', cint='free', manifestmeans='zero', traitvar=NULL, retryattempts=10))
hyperparams = list(a = list(manifestNames=c('sat6', 'attAnx'), 
                            waves=c(1,5),
                            N = c(1964),
                            minTpoints=2,
                            iterations = 100),
                   b = list(manifestNames=c('sat6', 'attAnx', 'attAvd'), 
                            waves=c(1,5),
                            N = c(1964),
                            minTpoints=2,
                            iterations = 35),
                   c = list(manifestNames=c('sat6', 'attAnx', 'attAvd', 'per1i6'), 
                            waves=c(1,5),
                            N = c(1964),
                            minTpoints=2,
                            iterations = 15),
                   d = list(manifestNames=c('sat6', 'attAnx', 'attAvd', 'per1i6', 'pattAnx'), 
                            waves=c(1,5),
                            N = c(1964),
                            minTpoints=2,
                            iterations = 7),
                   e = list(manifestNames=c('sat6', 'attAnx', 'attAvd', 'per1i6', 'pattAnx', 'pattAvd'), 
                            waves=c(1,5),
                            N = c(1964),
                            minTpoints=2,
                            iterations = 3))
res7 = performance_test(hyperparams, models, data.pairfam)
res7




#--- PART 3: COMPARABILITY OF PARAMETERS ACROSS INTERVALS --------------------------------------------------------------------------------------------------------------------#

models = list(voelkle       = list(type='voelkle',  approx='precise'),
              discrete_fiml = list(type='discrete', missing='fiml'))

# Generate a list with all possible combinations of waves
wave_instances = list()
for (i in seq(from=1, to=11, by=2)) {
  j=i+2
  while (j<=11) {
    wave_instances[[length(wave_instances)+1]] = c(i,j)
    k=j+2
    while (k<=11) {
      wave_instances[[length(wave_instances)+1]] = c(i,j,k)
      l=k+2
      while (l<=11) {
        wave_instances[[length(wave_instances)+1]] = c(i,j,k,l)
        m=l+2
        while (m<=11) {
          wave_instances[[length(wave_instances)+1]] = c(i,j,k,l,m)
          m=m+2
        }
        l=l+2
      }
      k=k+2
    }
    j=j+2
  }
}
manifestNames=c('sat6', 'attAnx')
res8a=data.frame()
for(waves in wave_instances) {
  res = modComp(models, dat.trainVal=prepData(data.pairfam, manifestNames, keepOnlyWaves=waves), manifestNames=manifestNames)
  res['waves']=paste(unlist(waves), collapse=', ')
  res8a = rbind(res8a, res)
}
res8a[res8a['paramType']=='disc',]

writeRes()



#--- PART 4: EXAMINATION OF PSYCHOLOGICAL HYPOTHESES --------------------------------------------------------------------------------------------------------------------#

models = list(naive             = list(type='naive'),
              voelkle           = list(type='voelkle',  approx='precise'),
              ctsem             = list(type='omx', cint='free', manifestmeans='zero', traitvar=NULL),
              autoregressive    = list(type='simple',   missing='fiml'),
              discrete_fiml     = list(type='discrete', missing='fiml'))

manifestNames=c('sat6', 'per1i6', 'attAnx', 'attAvd', 'pattAnx', 'pattAvd') # All variables considered now

# Create 80%/20% split for male subjects
dat9.m = prepData(data.pairfam, manifestNames, genders=1)
dat9.m.trainVal = trainTestSplit(dat9.m, test_size=0.2)$train
dat9.m.test = trainTestSplit(dat9.m, test_size=0.2)$test

# Create 80%/20% split for female subjects
dat9.f = prepData(data.pairfam, manifestNames, genders=2)
dat9.f.trainVal = trainTestSplit(dat9.f, test_size=0.2)$train
dat9.f.test = trainTestSplit(dat9.f, test_size=0.2)$test

# Combine both subsets
dat9 = rbind(dat9.m, dat9.f)
dat9.trainVal = rbind(dat9.m.trainVal, dat9.f.trainVal)
dat9.test = rbind(dat9.m.test, dat9.f.test)

# Test what model works best for that particular use case
res9a = modComp(models, dat9.trainVal, manifestNames, k=5)
res9a2 = modComp(models, dat9.trainVal, manifestNames, k=5)
res9a3 = modComp(models, dat9.trainVal, manifestNames, k=5)

# Evaluate prediction accuracy for chosen models
best_model = models$ctsem
res9b = modEval(dat9.trainVal, manifestNames, mod_type=best_model$type, hyperparams=best_model, dat.test=dat9.test)

# Obtain final results
res9c = modEval(dat9, manifestNames, mod_type=best_model$type, hyperparams=best_model, final_fit=TRUE)

# Analyze results separately for both genders
res9b.m2 = modEval(dat9.m, manifestNames, mod_type=best_model$type, hyperparams=best_model, dat.test=dat9.m.test)
res9b.f = modEval(dat9.f.trainVal, manifestNames, mod_type=best_model$type, hyperparams=best_model, dat.test=dat9.f.test)


i=2
drift = as.matrix(res9b[res9b$paramType=='cont',c(manifestNames)])
startVals = c(0,0,-1,-1,-1,-1)
expm(drift*i) %*% startVals


writeRes()

