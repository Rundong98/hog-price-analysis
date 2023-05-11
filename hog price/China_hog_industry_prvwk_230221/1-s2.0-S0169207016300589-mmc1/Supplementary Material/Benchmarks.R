###########################################################################
# BENCHMARKS USED IN CONSUMPTION FORECAST SECTION OF PAPER                #
# Wilms I. and CROUX C. (2016), "Forecasting using sparse cointegration." #
###########################################################################

SparseVAR_Lasso<-function(p,Y,X,Z,r,alpha.init=NULL,beta.init,max.iter=25,conv=10^-3,lambda.gamma=matrix(seq(from=0.1,to=0.001,length=10),nrow=1),lambda_beta=matrix(seq(from=2,to=0.001,length=100),nrow=1),rho.glasso=0.5,cutoff=0.8){
  ### FUNCTION TO PERFORM PENALIZED MAXIMUM LIKELIHOOD OF VAR ###
  
  ## INPUT
  # p: number of lagged differences
  # Y: Response Time Series
  # X: Time Series in Differences
  # Z: Time Series in Levels
  # r: cointegration rank
  # alpha.init: initial value for adjustment coefficients
  # beta.init: initial value for cointegrating vector
  # max.iter: maximum number of iterations
  # conv: convergence parameter
  # lambda.gamma: tuning paramter short-run effects
  # lambda_beta: tuning paramter cointegrating vector
  # rho.glasso: tuning parameter inverse error covariance matrix
  # cutoff: cutoff value time series cross-validation approach
  
  ## OUTPUT
  # ZBETA: estimate of short-run effects
  # OMEGA: estimate of inverse covariance matrix
  
  ### START CODE
  
  # Preliminaries
  q <- dim(Y)[2]
  n <- dim(Y)[1]
  
  # Initial value
  if (is.null(alpha.init)){
    Pi.init<-matrix(0,ncol=q,nrow=q)
  } else {
    Pi.init<-alpha.init%*%t(beta.init)
  }
  
  # Convergence parameters: initialization
  it<-1
  diff.obj<-10*conv
  Omega.init=diag(1,q)
  value.obj<-matrix(NA,ncol=1,nrow=max.iter+1)
  RESID<-Y-X%*%matrix(rep(diag(1,q),p-1),ncol=q,byrow=T)
  value.obj[1,]<-(1/n)*sum(diag(RESID%*%Omega.init%*%t(RESID)))-log(det(Omega.init))
  
  
  while( (it<max.iter) &  (diff.obj>conv) )
  {
    
    # Obtain gamma
    FIT1<-NTS.GAMMA.LassoReg(Y=Y,X=X,Z=Z,Pi=Pi.init,p=p,Omega=Omega.init,lambda.gamma=lambda.gamma)
    Resid<-Y-X%*%FIT1$ZBETA
    covResid<-cov(Resid)
    # Obtain omega
    GLASSOfit<-huge(covResid,lambda=rho.glasso,method="glasso",cov.output=T,verbose=F)
    OMEGA<-GLASSOfit$icov[[1]]

    
    # Check convergence
    value.obj[1+it,]<- (1/n)* sum(diag((Resid)%*%OMEGA%*%t(Resid))) - log(det(OMEGA))
    diff.obj<-abs(value.obj[1+it,]-value.obj[it,])/abs(value.obj[it,])
    Omega.init<-OMEGA
    it<-it+1
  }
  
  out<-list(it=it,ZBETA=FIT1$ZBETA,OMEGA=OMEGA)
}

DFM<-function(Y,X,rank=NULL,p=1,Ylag){
  # Function to estimate a Dynamic Factor Model for stationary time series
  
  ### INPUT
  # Y: Response 
  # X: Design variables
  # p: order of VAR model
  # Ylag: lagged variables
  
  ### OUTPUT
  # Betahat: estimated regression coefficients 
  # PCloading: loadings of principal component analysis
  # rank: number of factor
  
  # PCA
  decomp<-eigen(cov(X))
  if(is.null(rank)){
    MER<-decomp$values[-length(decomp$values)]/decomp$values[-1] #maximum eigenvalue ratio
    rank<-which.max(MER)
  }
  Factor<-X%*%decomp$vector[,1:rank]
  

  # Estimate VAR Model with Factor
  Response <- c(Y) # Response
  q<-ncol(Y)
  Y.lag<-as.matrix(Ylag)
  Xmatrix.kron <- cbind(Y.lag,Factor)
  # Fit by OLS
  BETA.OLS<-ginv(t(Xmatrix.kron)%*%Xmatrix.kron)%*%t(Xmatrix.kron)%*%Response 
  
  
  out<-list(Betahat=BETA.OLS,PCloading=decomp$vector[,1:rank],rank=rank)
}

LBVAR<-function(Y,X,p,delta=rep(0,M),lambda=sqrt(0.05),epsilon=10^-5){
  
  ## References:  Banbura, Giannone, Reichlin (2010)
  #              "Large Bayesian Vector Auto Regressions"
  #               Kadiyala and Karlsson (1997)
  #               "Numerical Methods for estimation and inference in Bayesian VAR-models"
  
  ## INPUT: 
  # Y=vector of endogenous variables
  # p=lag length
  # delta: contains M elements, prior on own lags for lag 1, default value see Kadiyala (1997)
  # lambda: hyperparamter in prior covariance matrix associated to own lags, default value see Kadiyala (1997)
  # epsilon: uninformative prior for the intercept, default value see Kadiyala (1997)
  
  ## OUTPUT: 
  # sigma.matrix: covariance matrix of error terms
  # post.coef: post mean matrix of VAR coefficients
  # post.cov.B: post covariance matrix of coefficients
  # post.cov.df: degrees of freedom post Wishart distribution of covmatrix
  
  
  
  # Parameter and Variable setting
  M<-ncol(Y) #Number of endogenous variables
  kx<-ncol(X)
  K<-kx #Number of parameters in each equation
  vec<-c(Y);y=matrix(data=vec,ncol=1) #Stack Y variables
  X=cbind(X,rep(1,nrow(X)))
  t=nrow(Y)
  
  #STEP 1: Covariance matrix of error terms
  sigma=matrix(data=0,ncol=M,nrow=M) #Store results from AR(p)-model
  for (i in 1:M){
    fit=ar.ols(Y[,i],aic=F,order.max=p)
    sigma[i,i]=fit$var.pred 
  }
  sigma.diag=matrix(data=diag(sigma),ncol=M,nrow=1);
  
  
  #STEP 2: Creating  dummy variables
  #Y variable
  Yd=matrix(data=0,ncol=M, nrow=M+M*(p-1)+M+1)
  
  #Priors autoregressive coefficients
  Y.BLOCK1=matrix(data=0,ncol=M,nrow=M)              
  diag(Y.BLOCK1)=delta*sqrt(sigma.diag)/lambda       
  Y.BLOCK2=matrix(data=0,ncol=M,nrow=M*(p-1))                 
  
  #Priors covariance matrix
  Y.BLOCK3=matrix(data=0,ncol=M,nrow=M)
  diag(Y.BLOCK3)=sqrt(sigma.diag)
  
  #Uninformative prior intercept
  Y.BLOCK4=matrix(data=0,nrow=1,ncol=M)
  
  #Matrix of dummy variable Y
  Yd=rbind(Y.BLOCK1,Y.BLOCK2,Y.BLOCK3,Y.BLOCK4)
  
  
  #X variable
  Xd=matrix(data=0,ncol=p*M+1,nrow=p*M +M + 1)
  
  #Priors autoregressive coefficients
  Jp=matrix(data=0,nrow=p,ncol=p)
  diag(Jp)=rep(1:p)
  sigma.lambda=sqrt(sigma)/lambda
  kron.prod=kronecker(Jp,sigma.lambda)
  X.BLOCK1=cbind(kron.prod,matrix(data=0,ncol=1,nrow=M*p))
  
  #Priors covariance matrix
  X.BLOCK2=matrix(data=0,ncol=p*M+1,nrow=M)
  
  #Uninformative prior intercept
  X.BLOCK3=matrix(data=0,ncol=p*M+1,nrow=1)
  X.BLOCK3[1,(p*M+1)]=epsilon
  
  #Matrix of dummy variable X
  Xd=rbind(X.BLOCK1,X.BLOCK2,X.BLOCK3)
  
  #New variables augmented regression model
  Y.star=rbind(Y,Yd)
  X.star=rbind(X,Xd)
  
  
  
  #STEP 3: Posterior mean of B
  B.postmean=ginv(t(X.star)%*%X.star)%*%t(X.star)%*%Y.star
  
  
  #STEP 4: Posterior mean of covmatrix error terms
  S.postmean=t(Y.star - X.star%*%B.postmean)%*%(Y.star - X.star%*%B.postmean)
  
  
  #STEP 5: Posterior covariance of B
  B0=ginv(t(Xd)%*%Xd)%*%t(Xd)%*%Yd
  sigma0=t(Yd-Xd%*%B0)%*%(Yd-Xd%*%B0)
  B.postcov=kronecker(S.postmean,ginv(t(X.star)%*%X.star))
  
  #STEP6: Degrees of freedom Wishart distribution
  W.df=nrow(Yd) + 2 + t - K
  
  # Transform VAR coefficients
  intercept=B.postmean[nrow(B.postmean),] #Vector of intercepts
  
  B=array(0,c(M,M,p)) #Array of p matrices associated to each lag
  for (i in 1:p){
    beginrow=1+(i-1)*M
    endrow=beginrow+ (M-1)
    B[,,i]=t(B.postmean[(beginrow:endrow),])
  }
  output=list(sigma.matrix=sigma,post.intercept=intercept,post.coef=B,post.cov.B=B.postcov,post.cov.df=W.df,post.mean.S=S.postmean)
  return(output)
}


BayesianRRP<-function(Bposterior,r=NULL){
  ### Bayesian reduced rank regression ###
  
  ### INPUT
  # Bposterior: MxN posterior mean of B
  # r: rank
  
  ### OUTPUT
  # Bpostr: rank r approximation of posterior mean
  # r: rank
  
  ### START CODE
  SVDdecomp<-svd(Bposterior) # Singular value decomposition
  MER<-SVDdecomp$d[-length(SVDdecomp$d)]/SVDdecomp$d[-1] #maximum eigenvalue ratio
  if(is.null(r)==T){
    r<-which.max(MER)
  } 
  Bpostr<-SVDdecomp$u[,1:r]%*%diag(SVDdecomp$d[1:r],r)%*%t(SVDdecomp$v[,1:r])
  
  # RETURN
  out<-list(Bpostr=Bpostr,rank=r)
}

DFMLEVELS<-function(DX,X,rank){
  ### BARIGOZZI ET AL 2016, Dynamic Factor Model for non-stationary time series
  
  # INPUT
  # DX: time series in differences
  # X: time series in levels
  # rank: rank
  
  # OUTPUT
  # LAMBDA: estimated factor loadings
  # ALPHA: estimated alpha
  # BETA: estimated beta
  # GAMMA: estimated short-run dynamics
  # FACTOR: estimatedfactors
  # rank: number of factors
  
  # CALCULATE LOADINGS
  decomp<-eigen(cov(DX))
  if(is.null(rank)){
    MER<-decomp$values[-length(decomp$values)]/decomp$values[-1] #maximum eigenvalue ratio
    rank<-which.max(MER)
  }
  LAMBDA<-decomp$vector[,1:rank]
  # FACTORS
  FACTOR<-X%*%LAMBDA/ncol(X)
  FACTOR.DIFF<-diff(FACTOR)
  
  # COINTEGRATION MODEL FOR FACTORS
  DATA.FACTOR<-embed(FACTOR.DIFF,dimension=2)
  FACTOR.Y<-DATA.FACTOR[,1:ncol(FACTOR)]
  FACTOR.X<-DATA.FACTOR[,-(1:ncol(FACTOR))]
  FACTOR.Z<-FACTOR[2:(nrow(FACTOR)-1),]
  
  # OBTAIN ESTIMATES VIA CCA ANALYSIS
  RESID0<-FACTOR.Y-FACTOR.X%*%ginv(t(FACTOR.X)%*%FACTOR.X)%*%t(FACTOR.X)%*%FACTOR.Y
  RESID1<-FACTOR.Z-FACTOR.X%*%ginv(t(FACTOR.X)%*%FACTOR.X)%*%t(FACTOR.X)%*%FACTOR.Z
  
  M02<-(t(FACTOR.Y)%*%FACTOR.X)/nrow(RESID0)
  M22<-(t(FACTOR.X)%*%FACTOR.X)/nrow(RESID0)
  M12<-(t(FACTOR.Z)%*%FACTOR.X)/nrow(RESID0)
  M01<-(t(FACTOR.Y)%*%FACTOR.Z)/nrow(RESID0)
  M11<-(t(FACTOR.Z)%*%FACTOR.Z)/nrow(RESID0)
  S01<-M01-M02%*%ginv(M22)%*%t(M12)
  S11<-M11-M12%*%ginv(M22)%*%t(M12)
  
  cancor_fit=cancor(RESID0,RESID1,xcenter=F,ycenter=F)
  BETA<-matrix(cancor_fit$ycoef[,1],ncol=1)
  ALPHA<-S01%*%BETA%*%ginv(t(BETA)%*%S11%*%BETA)
  GAMMA<-M02%*%ginv(M22)-ALPHA%*%t(BETA)%*%M12%*%ginv(M22)
  
  
  out<-list("LAMBDA"=LAMBDA,"ALPHA"=ALPHA,"BETA"=BETA,"GAMMA"=GAMMA,"FACTOR"=FACTOR,"rank"=rank)
}

