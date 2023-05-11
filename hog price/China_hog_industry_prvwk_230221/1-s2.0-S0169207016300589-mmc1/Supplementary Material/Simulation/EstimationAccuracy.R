####################################################################################################
### R-SCRIPT SIMULATION STUDY: ESTIMATION ACCURACY LOW-DIMENSIONAL SPARSE DESIGN q=4, T=500, r=1 ###
####################################################################################################
rm(list=ls())

# LIBRARIES
library(huge)
library(mvtnorm)
library(MASS)
library(glmnet)
library(glasso)
library(lars)
library(urca)
library(R.utils)
library(pls)

# SET WORKING DIRECTORY AND LOAD FUNCTIONS
setwd("W:/SparseCointegration/Simulation")
source('../Functions.R')

# PARAMETER SETTINGS
q<-4 # number of time series
p<-2 # order of VAR in levels
n<-500 # sample size
burnin<-100 # burn-in length
varcov<-diag(1,q) # error covariance matrix
for(i.var in 1:q){
  for(j.var in 1:q){
    varcov[i.var,j.var]<-0.0^(abs(i.var-j.var)) 
  }
}
beta0<-matrix(c(1,0,0,0),ncol=1) # cointegrating vector
r<-ncol(beta0) # cointegration rank
a.values<-seq(from=-0.2,to=-0.8,by=-0.2) # adjustment coefficients
GAMMA1<-diag(0.1,q) # short run effects

# STORE RESULTS FOR DIFFERENT ADJUSTMENT COEFFICIENTS
angle_ML<-matrix(NA,ncol=length(a.values),nrow=1)
angle_PML<-matrix(NA,ncol=length(a.values),nrow=1)
ITERATION_PML<-array(NA,c(1,length(a.values),1))
TIME_PML<-array(NA,c(1,length(a.values),1))

# PARALLELL COMPUTING: LOAD REQUIRED PACKAGES
library(foreach)
library(doSNOW)
ncl<-parallel:::detectCores()-1 # check number of cores
R<-500 # number of simulations
simestcl <- makeCluster(ncl) 
registerDoSNOW(simestcl)


AccuracyLD_R1 <- foreach(isim = 1:R, .packages = c('huge','mvtnorm', 'MASS','glmnet','glasso','lars','urca','R.utils','pls')) %dopar% {
  for (ia in 1:length(a.values)){
    set.seed(isim*ia)
    
    # MATRIX ADJUSTMENT COEFFICIENTS
    a<-a.values[ia]
    alpha0<-a*matrix(c(1,1,0,0),ncol=1)
    
    ###########################
    # Data Generating Process #
    ###########################
    Y<-matrix(data=NA,ncol=q,nrow=n+burnin)
    Y[1:2,]<-0
    for (i in 3:(n+burnin)){
      Y[i,]<- matrix(data=Y[i-1,],ncol=1) + alpha0%*%t(beta0)%*%matrix(data=Y[i-1,],ncol=1) + GAMMA1%*%matrix(data=(Y[i-1,]-Y[i-2,]),ncol=1) + matrix(data=rmvnorm(1,mean=rep(0,q),sigma=varcov),ncol=1)
    }
    Y<-Y[-(1:burnin),]
    
    #################
    # ML estimation #
    #################   
    dsY<-embed(diff(Y),p)
    Z<-Y[p:(n-1),]
    Y<-dsY[,1:q]
    X<-dsY[,-(1:q)]

    BETA_DELTA<-ginv(t(X)%*%X)%*%t(X)%*%Y
    BETA_LEVEL<-ginv(t(X)%*%X)%*%t(X)%*%Z
    
    RESID_DELTA<-Y-X%*%BETA_DELTA
    RESID_LEVEL<-Z-X%*%BETA_LEVEL
    
    # Canonical Correlation Analysis  
    cancor_fit=cancor(RESID_DELTA,RESID_LEVEL,xcenter=F,ycenter=F)
    BETAHAT_JOH<-cancor_fit$ycoef[,1:r] 
    angle_ML[1,ia]<-principal_angles(BETAHAT_JOH,a=beta0) # Estimation accuracy cointegration space
    
    ##################
    # PML estimation #
    ################## 
    ptm <- proc.time()  # Start timinig
    FIT3<-try(SparseCointegration_Lasso(Y=Y,X=X,Z=Z,beta.init=NULL,alpha.init=NULL,p=p,r=r,max.iter=10,conv=10^-2,lambda.gamma=matrix(seq(from=0.01,to=0.001,length=5),nrow=1),rho.glasso=seq(from=1,to=0.1,length=5),lambda_beta=matrix(seq(from=0.1,to=0.001,length=100),nrow=1)))
    FIT.time<-proc.time() - ptm   # Stop timining # 
    
    ITERATION_PML[1,ia,1]<-FIT3$it # number of iterations
    TIME_PML[1,ia,1]<-FIT.time[1] # computing time

    BETAHAT_PML<-FIT3$BETAhat
    angle_PML[1,ia]<-principal_angles(BETAHAT_PML,a=beta0) # Estimation accuracy cointegration space
  }
  list("angle_PML"=angle_PML,"angle_ML"=angle_ML,"ITERATION_PML"=ITERATION_PML,"TIME_PML"=TIME_PML)  # return values
  
}
stopCluster(simestcl)

angleML<-matrix(NA,ncol=length(a.values),nrow=R)
colnames(angleML)<-paste("a=",a.values)
anglePML<-matrix(NA,ncol=length(a.values),nrow=R)
colnames(angleML)<-colnames(angleML)
iterationsPML<-matrix(NA,ncol=length(a.values),nrow=R)
colnames(iterationsPML)<-colnames(angleML)
timePML<-matrix(NA,ncol=length(a.values),nrow=R)
colnames(timePML)<-colnames(angleML)

for(i in 1:R){
  angleML[i,]<-AccuracyLD_R1[[i]]$angle_ML
  anglePML[i,]<-AccuracyLD_R1[[i]]$angle_PML
  iterationsPML[i,]<-AccuracyLD_R1[[i]]$ITERATION_PML
  timePML[i,]<-AccuracyLD_R1[[i]]$TIME_PML
}
round(apply(angleML,2,mean),3) # Estimation accuracy ML estimator
round(apply(anglePML,2,mean),3) # Estimation accuracy PML estimator
round(apply(iterationsPML,2,mean),2) # number of iterations until convergence PML estimator
round(apply(timePML,2,mean),2) # computation time PML esitmator
