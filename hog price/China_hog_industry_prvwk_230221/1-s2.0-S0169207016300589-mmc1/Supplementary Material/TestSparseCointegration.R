### R-SCRIPT TEST CODE SPARSE COINTEGRATION ###

rm(list=ls())
library(mvtnorm)
library(glmnet)
library(R.utils)
library(pls)
library(MASS)
library(parcor)
library(urca)
library(glasso)
library(lars)
library(dplyr)
# ADJUST  WORKING DIRECTORY
setwd("D:/Research/hog price/China_hog_industry_prvwk_230221/1-s2.0-S0169207016300589-mmc1/Supplementary Material")
source('Functions.R')

data = read.csv("China_Hog_PriceAnalysis.csv")
data = data %>% select(-c(X,Date))

# DATA SETTINGS
set.seed(1)
q<-81 # number of time series
p<-2 # order of VAR in levels
n<-287 # time series length
burnin<-100 # burn-in period
varcov<-diag(1,q) # error covariance matrix
beta0<-matrix(c(1,0,0,0),ncol=1) # cointegrating vector
r<-ncol(beta0) # cointegration rank
alpha0<--0.8*matrix(c(1,1,0,0),ncol=1) #adjustment coefficients
GAMMA1<-0.1*diag(1,q) # short-run effects

# Data Generating Process #
    
Y<-matrix(data=NA,ncol=q,nrow=n+burnin)
Y[1:2,]<-0
for (i in 3:(n+burnin)){
 Y[i,]<- matrix(data=Y[i-1,],ncol=1) + alpha0%*%t(beta0)%*%matrix(data=Y[i-1,],ncol=1) + GAMMA1%*%matrix(data=(Y[i-1,]-Y[i-2,]),ncol=1) + matrix(data=rmvnorm(1,mean=rep(0,q),sigma=varcov),ncol=1)
}
Y<-Y[-(1:burnin),]
    
# Johansen's ML 
Y<-data %>% as.matrix()
dsY<-embed(diff(Y),p)
Z<-Y[p:(n-1),]
Y<-dsY[,1:q]
X<-dsY[,-(1:q)]
BETA_DELTA<-ginv(t(X)%*%X)%*%t(X)%*%Y
BETA_LEVEL<-ginv(t(X)%*%X)%*%t(X)%*%Z
RESID_DELTA<-Y-X%*%BETA_DELTA
RESID_LEVEL<-Z-X%*%BETA_LEVEL
cancor_fit=cancor(RESID_DELTA,RESID_LEVEL,xcenter=F,ycenter=F)
BETAHAT_JOH<-cancor_fit$ycoef[,1:r] 
angle_JOH<-principal_angles(BETAHAT_JOH,a=beta0)

# Sparse Cointegration
ptm <- proc.time()  # Start timinig
FITPML<-try(SparseCointegration_Lasso(Y=Y,X=X,Z=Z,beta.init=NULL,alpha.init=NULL,p=p,r=r,max.iter=10,conv=10^-2,lambda.gamma=matrix(seq(from=0.01,to=0.001,length=5),nrow=1),rho.glasso=seq(from=1,to=0.1,length=5),lambda_beta=matrix(seq(from=0.1,to=0.001,length=100),nrow=1)))
FIT.time<-proc.time() - ptm   # Stop timining # 
BETAHAT_PML<-FITPML$BETAhat
angle_PML<-principal_angles(BETAHAT_PML,a=beta0) 
FIT.time[1] # computation time
FITPML$it # number of iterations

re = determine_rank(Y,X,Z,r.init=NULL,p,max.iter.lasso=3,conv.lasso=10^-2,max.iter.r=5,beta.init=NULL,alpha.init=NULL,rho.glasso=0.1,lambda.gamma=matrix(seq(from=0.1,to=0.001,length=10),nrow=1),lambda_beta=matrix(seq(from=2,to=0.001,length=100),nrow=1),glmnetthresh=1e-04)
