#################################################################################################
### R-SCRIPT SIMULATION STUDY RANK DETERMINATION: LOW-DIMENSIONAL SPARSE DESIGN q=4, T=500, r=1 #
#################################################################################################
rm(list=ls())

# LIBRARIES
library(mvtnorm)
library(MASS)
library(pls)
library(parcor)
library(glmnet)

# SET WORKING DIRECTORY AND LOAD FUNCTIONS
setwd("W:/SparseCointegration/Simulation")
source('..Functions.R')

# PARAMETER SETTINGS
Nboot<-399
Nsim<-500
p<-4
k<-2 # order of VAR in levels
n<-500 # sample size
burnin<-100
varcov<-diag(1,p)
beta<-matrix(c(1,0,0,0),ncol=1)
r<-ncol(beta)
alpha<--0.4*matrix(c(1,1,0,0),ncol=1)
GAMMA1<-diag(0.1,p)

# Determination of cointegration rank
Johansen_rank=matrix(data=0,ncol=p+1,nrow=Nsim)
Bootstrap_rank=matrix(data=0,ncol=p+1,nrow=Nsim)
Bartlett_rank=matrix(data=0,ncol=p+1,nrow=Nsim)
RSC=matrix(data=0,ncol=p+1,nrow=Nsim)

# Critical values Johansen approach
cv<-c(48.28, 31.52 , 17.95   ,8.18) # Critical values Johansen trace test
for (nsim in 1:Nsim){
  cat("Start",nsim,"\n")
  
  ###########################
  # Data Generating Process #
  ###########################
  X<-matrix(data=NA,ncol=p,nrow=n+burnin)
  X[1:2,]<-0
  for (i in 3:(n+burnin)){
    X[i,]<- matrix(data=X[i-1,],ncol=1) + alpha%*%t(beta)%*%matrix(data=X[i-1,],ncol=1) + GAMMA1%*%matrix(data=(X[i-1,]-X[i-2,]),ncol=1) + matrix(data=rmvnorm(1,mean=rep(0,p),sigma=varcov),ncol=1)
  }
  X<-X[-(1:burnin),]
  
  #################
  # ML ESTIMATION #
  #################
  DELTA_X<-diff(X,differences=1)
  LEVEL_X<-X[-nrow(X),]
  data<-embed(DELTA_X,dimension=k)
  
  DELTA_X_AUX<-data[,(1:p)]
  LEVEL_X_AUX<-cbind(LEVEL_X[-(1:(k-1)),])
  DESIGN_X_AUX<-cbind(data[,-(1:p)])
  
  BETA_DELTA<-solve(t(DESIGN_X_AUX)%*%DESIGN_X_AUX)%*%t(DESIGN_X_AUX)%*%DELTA_X_AUX
  BETA_LEVEL<-solve(t(DESIGN_X_AUX)%*%DESIGN_X_AUX)%*%t(DESIGN_X_AUX)%*%LEVEL_X_AUX
  
  RESID_DELTA<-DELTA_X_AUX-DESIGN_X_AUX%*%BETA_DELTA
  RESID_LEVEL<-LEVEL_X_AUX-DESIGN_X_AUX%*%BETA_LEVEL
  
  # Canonical Correlation Analysis
  cancor=cancor(RESID_DELTA,RESID_LEVEL)
  eigenvalues=cancor$cor^2
  
  # JOHSANSEN TRACE STATISTIC
  for (r in 0:(p-1)){
    test_product<--n*log(1-eigenvalues[(r+1):p])
    test_r<-sum(test_product)
    if (test_r<cv[r+1] & sum(Johansen_rank[nsim,])==0) {Johansen_rank[nsim,r+1]=1} 
    else {r=r+1}
  }
  if (sum(Johansen_rank[nsim,])==0){Johansen_rank[nsim,p+1]=1}
  
  
  # BARTLETT CORRECTED TRACE STATISTIC
  S11<-(t(RESID_LEVEL)%*%RESID_LEVEL)/n
  S10<-(t(RESID_LEVEL)%*%RESID_DELTA)/n
  S01<-(t(RESID_DELTA)%*%RESID_LEVEL)/n
  S00<-(t(RESID_DELTA)%*%RESID_DELTA)/n
  S11_eigenvalues<-eigen(S11)$values
  S11_eigenvectors<-eigen(S11)$vectors 
  S11_root<-S11_eigenvectors%*%diag(1/sqrt(S11_eigenvalues))%*%t(S11_eigenvectors)
  matrix<-S11_root%*%S10%*%solve(S00)%*%S01%*%S11_root
  decomp<-eigen(matrix)
  V<-S11_root%*%decomp$vectors 
  
  
  Z0t<-DELTA_X_AUX
  Z1t<-LEVEL_X_AUX
  Z2t<-DESIGN_X_AUX
  M02<-(t(Z0t)%*%Z2t)/n
  M22<-(t(Z2t)%*%Z2t)/n
  M12<-(t(Z1t)%*%Z2t)/n
  
  for (r in 0:(p-1)){
    test_product<--n*log(1-eigenvalues[(r+1):p])
    test_r<-sum(test_product)
    Bartlettcorection<-Bartlettcorrection(r=r,n=n,p=p,k=k,cancor=cancor,S01=S01,V=V,S11=S11,M02=M02,M22=M22,M12=M12,S00=S00)$correctionterm
    testbartlett_r<-test_r/Bartlettcorection
    if (testbartlett_r<cv[r+1] & sum(Bartlett_rank[nsim,])==0) {Bartlett_rank[nsim,r+1]=1} 
    else {r=r+1}
  }
  if (sum(Bartlett_rank[nsim,])==0){Bartlett_rank[nsim,p+1]=1}
  
  # BOOTSTRAP RANK DETERMINATION
  test_stat_Johansen<-matrix(NA,ncol=1,nrow=p)
  for (ip in 1:p){
    test_stat_Johansen[ip,1]<-sum(-n*log(1-eigenvalues[((ip-1)+1):p]))
  }
  for (rank in 0:(p-1)){
    bootrun=BOOTSTRAP(r=rank,V=V,p=p,M02=M02,M22=M22,M12=M12,LEVEL_X_AUX=LEVEL_X_AUX,DESIGN_X_AUX=DESIGN_X_AUX,DELTA_X_AUX=DELTA_X_AUX,Nboot=Nboot,test_stat_Johansen)
    if (bootrun$Bootstrap_flag>0){Bootstrap_error[nsim,1]=1}
    if (bootrun$fraction>0.05){Bootstrap_rank[nsim,(rank+1)]=1}
    if (sum(Bootstrap_rank[nsim,])==1){break}
  }
  if (sum(Bootstrap_rank[nsim,])==0){Bootstrap_rank[nsim,p+1]=1}
  
  # RANK SELECTION CRITERION
  Y<- DELTA_X_AUX
  X<- DESIGN_X_AUX
  Z<- LEVEL_X_AUX
  RSC_rhat<-determine_rank(Y=Y,X=X,Z=Z,r.init=p,p=k,beta.init=beta,alpha.init=alpha)
  RSC[nsim,RSC_rhat$rhat+1]<-1  
}

RANK_RESULTS<-matrix(NA,ncol=5,nrow=4)
RANK_RESULTS[1,]<-apply(Johansen_rank[,],2,mean)*100 
RANK_RESULTS[2,]<-apply(Bartlett_rank[,],2,mean)*100
RANK_RESULTS[3,]<-apply(Bootstrap_rank[,],2,mean)*100
RANK_RESULTS[4,]<-apply(RSC[,],2,mean)*100 
colnames(RANK_RESULTS)<-paste("r=",rep(0:p))
rownames(RANK_RESULTS)<-c("Johansen","Bartlett","Bootstrap","RSC")
round(RANK_RESULTS,2)


