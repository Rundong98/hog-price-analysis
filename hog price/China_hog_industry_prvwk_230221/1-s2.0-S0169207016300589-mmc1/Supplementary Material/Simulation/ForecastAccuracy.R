####################################################################################################
### R-SCRIPT SIMULATION STUDY: FORECAST ACCURACY LOW-DIMENSIONAL SPARSE DESIGN q=4, T=500, r=1 ###
####################################################################################################
rm(list=ls())

# LIBRARIES
library(mvtnorm)
library(glmnet)
library(R.utils)
library(pls)
library(MASS)
library(parcor)
library(urca)
library(glasso)
library(lars)

# SET WORKING DIRECTORY AND LOAD FUNCTIONS
setwd("W:/SparseCointegration/Simulation")
source('../Functions.R')


# PARAMETER SETTINGS
q<-4 # number of time series
p<-2 # order of VAR in levels
nlength<-500 # sample size
burnin<-100 # burn-in length
varcov<-diag(1,q) # error covariance matrix
beta0<-matrix(c(1,0,0,0),ncol=1) # cointegrating vector
r<-ncol(beta0) # cointegration rank
a<--0.4  # adjustment coefficient
GAMMA1<-0.1*diag(1,q) # short-run effects

hseq<-c(1,3,6,12) #DIFFERENT FORECAST HORIZONS h
Sseq<-c(4,8,12) #DIFFERENT ROLLING WINDOW SIZES S=4*YRS

MMAFETABLEML<-matrix(NA,ncol=4,nrow=3)
MMAFETABLEPML<-matrix(NA,ncol=4,nrow=3)
colnames(MMAFETABLEPML)<-paste("h=",hseq)
rownames(MMAFETABLEPML)<-paste("S=",Sseq*12)
colnames(MMAFETABLEML)<-colnames(MMAFETABLEPML)
rownames(MMAFETABLEML)<-rownames(MMAFETABLEPML)

# FORECAST FOR DIFFERENT FORECAST HORIZONS h
for(h in c(1,3,6,12)){
  cat("start",h,"\n")
  # FORECAST FOR DIFFERENT ROLLING WINDOW SIZES S=12*YRS
  for(YRS in c(4,8,12)){
    cat("start",YRS,"\n")
    h<-h
    YRS<-YRS
    i.fend<-nlength-YRS*12-h-1
    
    # PARALLELL COMPUTING
    library(foreach)
    library(doSNOW)
    ncl<-parallel:::detectCores()-1 # check number of cores
    R<-500# !! only for two simulation runs, change this to R=500 
    simestcl <- makeCluster(ncl) 
    registerDoSNOW(simestcl)
    
    Forecast <- foreach(isim = 1:R, .packages = c('huge','mvtnorm', 'MASS','glmnet','glasso','lars','urca','R.utils','pls')) %dopar% {
      set.seed(isim)
      alpha0<-a*matrix(c(1,1,0,0),ncol=1)
      

      ###########################
      # Data Generating Process #
      ###########################
      Y<-matrix(data=NA,ncol=q,nrow=nlength+burnin)
      Y[1:2,]<-0
      for (i in 3:(nlength+burnin)){
        Y[i,]<- matrix(data=Y[i-1,],ncol=1) + alpha0%*%t(beta0)%*%matrix(data=Y[i-1,],ncol=1) + GAMMA1%*%matrix(data=(Y[i-1,]-Y[i-2,]),ncol=1) + matrix(data=rmvnorm(1,mean=rep(0,q),sigma=varcov),ncol=1)
      }
      Y<-Y[-(1:burnin),]
      
      # STORE FORECASTS
      FORECAST<-array(NA,c(q,3,i.fend))
      
      # BEGIN ROLLING WINDOW FORECAST
      for(i.f in 1:i.fend){
        cat("start",i.f,"\n")
        # TRAINING DATA
        trainingDATA<-Y[(1+i.f-1):(YRS*12+2+(h-1)+i.f-1),1:q]
        trainingDATADIFF<-embed(diff(trainingDATA),p+h-1)
        trainingY<-trainingDATADIFF[,1:q]
        trainingZ<-trainingDATA[-(1:(p+h-1)),]
        trainingX<-trainingDATADIFF[,(ncol(trainingDATADIFF)-q+1):ncol(trainingDATADIFF)]
        
        # TEST DATA
        TESTYLEVEL<-as.matrix(Y[(YRS*12+2+(h-1)+i.f-1)+1,]) 
        
        FORECAST[,1,i.f]<-c(TESTYLEVEL)-trainingDATA[nrow(trainingDATA),] # true value
        LEVELYtest<-trainingDATA[nrow(trainingDATA),]
        DELTALAGYtest<-trainingDATA[nrow(trainingDATA),]-trainingDATA[nrow(trainingDATA)-1,]
        
        
        # PML
         FIT3<-try(SparseCointegration_Lasso(Y=trainingY,X=trainingX,Z=trainingZ,beta.init=NULL,alpha.init=NULL,p=p,r=r,max.iter=10,conv=10^-2,lambda.gamma=matrix(seq(from=0.1,to=0.01,length=5),nrow=1),lambda_beta=matrix(seq(from=1,to=0.1,length=10),nrow=1)))
         FORECAST[,2,i.f]<-FIT3$ALPHAhat%*%t(FIT3$BETAhat)%*%LEVELYtest + FIT3$ZBETA%*%DELTALAGYtest
        
        # ML
        RESID_DELTA<-trainingY-trainingX%*%solve(t(trainingX)%*%trainingX)%*%t(trainingX)%*%trainingY
        RESID_LEVEL<-trainingZ-trainingX%*%solve(t(trainingX)%*%trainingX)%*%t(trainingX)%*%trainingZ
        n<-nrow(trainingY)
        M02<-(t(trainingY)%*%trainingX)/n
        M22<-(t(trainingX)%*%trainingX)/n
        M12<-(t(trainingZ)%*%trainingX)/n
        M01<-(t(trainingY)%*%trainingZ)/n
        M11<-(t(trainingZ)%*%trainingZ)/n
        S01<-M01-M02%*%solve(M22)%*%t(M12)
        S11<-M11-M12%*%solve(M22)%*%t(M12)
  
        BETAJOH<-matrix(cancor(RESID_DELTA,RESID_LEVEL)$ycoef[,1:r],ncol=r)
        ALPHAJOH<-S01%*%BETAJOH%*%solve(t(BETAJOH)%*%S11%*%BETAJOH)
        GAMMAJOH<-M02%*%solve(M22)-ALPHAJOH%*%t(BETAJOH)%*%M12%*%solve(M22)
        FORECAST[,3,i.f]<-ALPHAJOH%*%t(BETAJOH)%*%LEVELYtest + GAMMAJOH%*%DELTALAGYtest
        
      }
      
      # CALCULATE ABSOLUTE SCALED FORECAST ERROR
      FORECASTERROR<-array(NA,c(q,2,i.fend))
      colnames(FORECASTERROR)<-c("PML","ML")
      for(i.t in 1:i.fend){
        FORECASTERROR[,,i.t]<-abs(FORECAST[,-1,i.t]-FORECAST[,1,i.t])
        for(i.var in 1:q){
          FORECASTERROR[i.var,,i.t]<-FORECASTERROR[i.var,,i.t]/apply(diff(as.matrix(Y)),2,sd)[i.var]
        }
      }
      FEADJ<-apply(FORECASTERROR,c(2,3),mean) #mean over all variables
      MMAFE<-round(apply(FEADJ,1,mean),2) #mean over all time points
      
      list("MMAFE"=MMAFE,"FORECAST"=FORECAST)  # return values
      }

    stopCluster(simestcl)
    
    MMAFE_ALL<-matrix(NA,nrow=R,ncol=2)
    colnames(MMAFE_ALL)<-c("PML","ML")
    for(i in 1:R){
      MMAFE_ALL[i,]<-Forecast[[i]]$MMAFE
    }
    apply(MMAFE_ALL,2,mean)
    print(apply(MMAFE_ALL,2,mean))
    i.horizon<-which(hseq==h)
    i.window<-which(Sseq==YRS)
    MMAFETABLEPML[i.window,i.horizon]<-apply(MMAFE_ALL,2,mean)[1]
    MMAFETABLEML[i.window,i.horizon]<-apply(MMAFE_ALL,2,mean)[2]
  }
}

MMAFETABLEPML # MMAFE PML estimator
MMAFETABLEML # MMAFE ML estimator