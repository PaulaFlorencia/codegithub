################################################################################################
################################################################################################
### This document contains all the functions that are neccesary for the                      ###
### utilisation of Tests_p1.R                                                                ###   
################################################################################################
################################################################################################

library(evd)

################################################################################################
############################ p12_t and p13_t computation #######################################

######  Kernels  ###############################################################################

dEpan <- function(x){
  ## Function of Epanechnikov density distribution
  k <- (3/4)* (1-x^2)
  k[-1>x] <- 0
  k[x>1] <- 0
  return (k)
} 

# Normal distribution
Kern_gauss <- function(x)exp(-x^2/2)

#### p12_t Non parametric computation ##############################################################

p12_NonPar <- function(X.mat,Z.mat,tt,t_eval,h,kern= dEpan){
  
  ## Computation of p12_t using kernels   
  ##
  ##    \hat{p}12_t = Wtj %*% GmZ_tj 
  ##
  ## for each time step t 
  ##
  ## Input :
  ## - Counterfactual an factual trajectories X_t and Z_t
  ## - tt, vector tt of size J containing the Z's time steps
  ## - evaluation vector t_eval (in practice t_eval = tt)
  ## - kernel bandwith h 
  ##
  ## Output :
  ## - vector of length J containing \hat{p}12_t
  ##
  ## Used in : IC, CV_error, P12_P13_estimation , matcovp12p13_t
  ##
  ## Requires : dEpan (or Kern_gauss)
  
  X.mat <- as.matrix(X.mat)
  Z.mat <- as.matrix(Z.mat)
  
  GmZ.mat <- matGZ_func(X.mat,Z.mat)
  N <- dim(GmZ.mat)[2]
  J <- dim(GmZ.mat)[1]
  Kij <- outer(t_eval,tt,function(zz,z) kern((zz - z) / h))
  W <- Kij / rowSums(Kij)
  p12.mat <- matrix(NA, ncol = N, nrow = length(t_eval))
  for(i in 1:N){
    GZ_2<-GmZ.mat[,i]
    p12.mat[,i] <- W %*% GZ_2
  }
  return( list("p12.mat"=p12.mat,"GmZ.mat"=GmZ.mat) ) 
}

#### p13_t and p1r_t Non parametric computation ##############################################################

p13_NonPar <- function(X.mat,Z.mat,tt,t_eval,h,kern= dEpan){

  ## Computation of p12_t using kernels   
  ##
  ##    \hat{p}12_t = Wtj %*% (GmZ_tj)^2 
  ##
  ## for each time step t 
  ##
  ## Input :
  ## - Counterfactual an factual trajectories X_t and Z_t
  ## - tt, vector tt of size J containing the Z's time steps
  ## - evaluation vector t_eval (in practice t_eval = tt)
  ## - kernel bandwith h 
  ##
  ## Output :
  ## - vector of length J containing \hat{p}1r_t
  ##
  ## Used in : P12_P13_estimation , matcovp12p13_t
  ##
  ## Requires : dEpan (or Kern_gauss)
  
  X.mat <- as.matrix(X.mat)
  Z.mat <- as.matrix(Z.mat)
  
  GmZ.mat <- matGZ_func(X.mat,Z.mat)
  N <- dim(GmZ.mat)[2]
  J <- dim(GmZ.mat)[1]
  Kij <- outer(t_eval,tt,function(zz,z) kern((zz - z) / h))
  W <- Kij / rowSums(Kij)
  p13.mat <- matrix(NA, ncol = N, nrow = length(t_eval))
  for(i in 1:N){
    GZ_2<-(GmZ.mat[,i])^2
    p13.mat[,i] <- W %*% GZ_2
  }
  return(p13.mat)
}

p1r_NonPar <- function(X.mat,Z.mat,r.vec,tt,t_eval,h,kern= dEpan){
  
  ## Computation of p12_t using kernels   
  ##
  ##    \hat{p}1r_t = Wtj %*% (GmZ_tj)^(r-1) 
  ##
  ## for each time step t 
  ##
  ## Input :
  ## - Counterfactual an factual trajectories X_t and Z_t
  ## - tt, vector tt of size J containing the Z's time steps
  ## - evaluation vector t_eval (in practice t_eval = tt)
  ## - kernel bandwith h 
  ##
  ## Output :
  ## - vector of length J containing \hat{p}1r_t
  ##
  ## Used in : Not used in other functions
  ##
  ## Requires : dEpan (or Kern_gauss)  

  X.mat <- as.matrix(X.mat)
  Z.mat <- as.matrix(Z.mat)
  r.vec<-as.vector(r.vec)
  
  GmZ.mat <- matGZ_func(X.mat,Z.mat)
  N <- dim(GmZ.mat)[2]
  J <- dim(GmZ.mat)[1]
  Kij <- outer(t_eval,tt,function(zz,z) kern((zz - z) / h))
  W <- Kij / rowSums(Kij)
  p1r.mat <- matrix(NA, ncol = N, nrow = length(t_eval))
  for(i in 1:N){
    GZ_2<-(GmZ.mat[,i])^(r.vec[i]-1)
    p1r.mat[,i] <- W %*% GZ_2
  }
  return(p1r.mat)
}

#### Others fonctions ###############3########################################################################

# fonction that creates "real" factual run
Z_real_comp <- function(size,scale_Z_real,shape_Z_real) {
  mean=0; coef=1 # parameters of gev
  Z_real<-rep(0,size)
  for (i in 1:size/4){
    Z_real[i]<-Z_real[i]+rgev(1,loc=mean, scale_Z_real, shape_Z_real)
  }
  for (i in 251 : size){
    Z_real[i]<-Z_real[i]+rgev(1,loc=mean+(coef*tt[i]-(size/4))/size, scale_Z_real, shape_Z_real)
  }
  return(Z_real)
}

################################### Confidence interval ###########################################

#### Auxiliary fonctions 
dEpan_2<- function(x){
  ## Function of Epanechnikov squared density distribution
  k <- ((3/4)* (1-x^2))^2
  k[-1>x] <- 0
  k[x>1] <- 0
  return (k)
}

kern_dens <- function(tt,t_eval,h,kern=dEpan) {
  # f(t) density distibution
  fh <- numeric(length(t_eval))
  for (i in seq_along(tt))
    for (j in seq_along(t_eval))
      fh[i] <- fh[i] + kern((t_eval[i] - tt[j])/h) 
    fh <- fh / ( h * length(tt))
    return(fh)
}

#### p12_t's Confidence interval computation 
IC <- function(X,Z,tt,t_eval,h,kern=dEpan) {
  
  ## Computation of p12_t's confidence interval   
  ##
  ## Input :
  ## - Counterfactual an factual vectors X and Z
  ## - tt, vector tt of length(Z) containing the Z's time steps
  ## - evaluation vector t_eval (in practice t_eval = tt)
  ## - kernel bandwith value h 
  ##
  ## Output :
  ## - two lists containining upper and lower bonds of the estimator 
  ##
  ## Used in : Not used in other functions
  ##
  ## Requires : p12_NonPar, dEpan_2, dEpan,  
  
  n <- length(Z)
  G_emp <- ecdf(X)
  G_Z <- G_emp(Z)
  
  p12_hat<-p12_NonPar(X,Z,tt,t_eval,h,kern)$p12.mat
  K22<-integrate(dEpan_2,-1,1)$value
  f_t<-kern_dens(tt,tt,h)
  sigma_GZ <- as.numeric(var(G_Z-p12_hat))
  VarA <- (sigma_GZ * K22 * f_t )/(n*h)
  
  Khj <- outer(t_eval, tt, function(zz,z) dEpan((zz - z) / h))
  Khi <- outer(t_eval, tt, function(zz,z) dEpan((zz - z) / h))
  E <- mean(outer(G_Z,G_Z, pmin) - outer(G_Z,G_Z,"*"))
  VarB_num <- as.numeric(rowSums(Khj) %*% rowSums(Khi))
  VarB_denom <- length(tt)*(rowSums(Khj))^2
  VarB <- VarB_num*E/VarB_denom
  
  sigma_m <- VarA + VarB
  s <- (sigma_m)^1/2
  error <- qnorm(0.975)*s/sqrt(length(tt))
  left <- p12_hat-error
  right <- p12_hat+error
  return(list(low=left,high=right))
}

################################ Optimal bandwidth ################################################

#### Cross-validation Error

# as imput we have a sequence of h and as output d'error associated with each one
CV_error <- function(X,Z,tt,t_eval,h_seq,kern= dEpan){
  
  ## This routine uses cross validation to chose the optimal bandwith of the kernel   
  ## from two vectors.
  ##
  ## Input :
  ## - Counterfactual an factual vectors X and Z
  ## - tt, vector tt of length(Z) containing the Z's time steps
  ## - evaluation vector t_eval (in practice t_eval = tt)
  ## - h_seq, vector containing the set from which the optimal bandwith will be chosen 
  ## - kernel bandwith value h 
  ##
  ## Output : list of 3 elements
  ## - optimal bandwith
  ## - optimal error
  ## - vector of errors
  ##
  ## Used in : Not used in other functions yet 
  ##
  ## Requires : p12_Nonpar  
  
  n <- length(tt)
  G_emp <- ecdf(X)
  G_Z <- G_emp(Z)
  CV_err_h = rep(NA,length(h_seq))
  for(j in 1:length(h_seq)){
    h_using = h_seq[j]
    CV_err=rep(NA,n)
    for(i in 1:n){
      tt_test <- tt[i]; Z_test<-Z[i]; GZ_test <-G_Z[i] # validation set
      tt_tr <- tt[-i]; Z_tr<-Z[-i] # training set
      GZ_predict <-p12_NonPar(X,Z_tr,tt_tr,t_eval,h_using,kern)$p12.mat 
      CV_err[i]<-(GZ_test - GZ_predict[i])^2 
    }
    CV_err_h[j]<-(mean(CV_err))^0.5
  }
  h_opt <- h_seq[which(CV_err_h == min(CV_err_h))]
  opt_err <-min(CV_err_h)
  return(list("optimal_bandwith"= h_opt,"optimal_error"= opt_err ,"list_errors"= CV_err_h,"list_h"=h_seq) )#vector of errors
}

CVmat_error <- function(matX,matZ,tt,t_eval,mat_h_seq,kern= dEpan){
  
  ## Similar to CV_error(), but it allows Z and can be matrix. mat_hseq must be a matrix of the same number 
  ## of columns as matZ as matX. il allows to compare to find for each column a different optimal bandwith
  ## from a different set of possible bandwidths
  ##
  ## Used in : Optimal_h_cmip()
  
  matX <- as.matrix(matX)
  matZ <- as.matrix(matZ)
  GmZ.mat <- matGZ_func(matX,matZ)
  N <- dim(matZ)[2]
  J <- dim(matZ)[1]
  hopt.vec = rep(NA,N)
  opterr.vec = rep(NA,N)
  for(n in 1:N){
    G_Z<-GmZ.mat[,n]
    Z <- matz[,n]
    X <- matx[,n]
    h_seq = mat_h_seq[,n]
    CV_err_h = rep(NA,length(h_seq))
    for(j in 1:length(h_seq)){
      h_using = h_seq[j]
      CV_err=rep(NA,J)
      for(i in 1:J){
        tt_test <- tt[i]; Z_test<-Z[i]; GZ_test <-G_Z[i] # validation set
        tt_tr <- tt[-i]; Z_tr<-Z[-i] # training set
        GZ_predict <-p12_NonPar(X,Z_tr,tt_tr,t_eval,h_using,kern)$p12.mat 
        CV_err[i]<-(GZ_test - GZ_predict[i])^2 
      }
      CV_err_h[j]<-(mean(CV_err))^0.5
    }
    hopt.vec[n] <- h_seq[which(CV_err_h == min(CV_err_h))]
    opterr.vec[n] <-min(CV_err_h)
  }
  return(list("optimal_bandwith"= hopt.vec,"optimal_error"= opterr.vec) )
}

Optimal_h_cmip<- function(matx,matz,tt,t_eval,kern=dEpan){
  
  ## for each columns of our matrix matx and matz, it applies two times the fonction CVmat_error()
  ## and finds the optimal bandwidth of each column.
  ##
  ## Uses : CVmat_error()
  ##
  ## Used in :
  
  N <- dim(matz)[2]
  J <- dim(matz)[1]
  hseq_first <- c(7,10,20,30,40,50,60,70,80,90,100, 110, 120, 130,140,150,160,170,180,190,200)
  mat_h_first <- cbind(hseq_first,hseq_first)
  h_opt_first <- CVmat_error(matx,matz,tt,tt,mat_h_first,kern)$optimal_bandwith
  
  mat_h <- matrix(numeric(10*J), nrow = 10, ncol = N) # empty matrix
  for (i in 1:N){
    h_opt_i <- h_opt_first[i]
    hseqfinal <- c(h_opt_i-5,h_opt_i-4,h_opt_i-3,h_opt_i-2,h_opt_i-1,h_opt_i,h_opt_i+1,h_opt_i+2,h_opt_i+3,h_opt_i+4)
    mat_h[,i] <- hseqfinal
  }
  testerror <-CVmat_error(matx,matz,tt,tt,mat_h,kern)
}


#### General fonction #############################################################################

# Estimate p12_t and p13_t for each trajectory {(X)_t,(Z)_t}
P12_P13_estimation <- function(X.mat,Z.mat,tt,t_eval,h,kern=dEpan){
  
  ## This routine computes both \hat{p}12_t and  \hat{p}13_t 
  ##
  ## Input :
  ## - Counterfactual an factual matrix  matX and matZ , where each column is a trayectory Xi or Zi
  ## - tt, vector tt of length(Z) containing the Z's time steps
  ## - evaluation vector t_eval (in practice t_eval = tt)
  ## - kernel bandwith value h 
  ##
  ## Output : list of 2 elements
  ## - matrix matp12 of \hat{p}12_t , where the column i is associated to {Xi,Zi}
  ## - same for matrix matp13 
  ##
  ## Used in : FastTestforp1r_gumbel 
  ##
  ## Requires : p12_Nonpar, p13_NonPar 
  
  X.mat <- as.matrix(X.mat)
  Z.mat <- as.matrix(Z.mat) 
  dimnsZ=dim(Z.mat)

  if (dimnsZ[1]!=length(tt)){
    stop("Z trajectory and tt must have the same length")
  }
  p12.mat<-p12_NonPar(X.mat,Z.mat,tt,tt,h)$p12.mat
  GmZ.mat<-p12_NonPar(X.mat,Z.mat,tt,tt,h)$GmZ.mat
  p13.mat<-p13_NonPar(X.mat,Z.mat,tt,tt,h)
  return(list("matp12"=p12.mat,"matp13"=p13.mat,"matGmZ"=GmZ.mat))
}

################################################################################################
############################ lambda_t and k_t computation ######################################
###############################  (multiple méthods) ############################################
################################################################################################

library(stats4); library(gmm); library(stats); library(np); library(EWGoF)

#### functions for E(G(Z)^j)
funcLaplace <- function(x,m,lam,k,a){ 
  ##  Utilitary function used in the integrate() statement in function 
  ## laplaceWeibull() defined below
  (1/a) * exp( -(m*lam/a^(1/k))*(-log(x))^(1/k) ) * x^(1/a - 1)
}

laplaceWeibull <- function(j,lambda,k,lowerbnd=10^(-6),upperbnd=1,fac=1,tol=10^(-5)){
  
  ## This function computes E(G(Z)^j) 
  ## for any given integer j, where W is a Weibull(lambda,k) variable.
  ## 
  ## The computation is based on integration rather than partial power series, 
  ##
  ## Input :
  ## - single values of lambda and k (in practice issued from weibullGMM_NonStationaire
  ##   or other estimation method) 
  ## - j value , in practice j = r-1
  ## - other optional parameters controll the way the numerical integration is conducted 
  ##
  ##
  ## Output :
  ## - the numerical evaluation of E(exp(-j*W)) when W~Weibull(lambda,k)
  ##
  ## Used in function : many functions in this file
  ## Requires : funcLaplace()
  ## 

  cat("lam=",lambda,",k=",k,"\n")
  vala=fac*(j*lambda)^k  
  upperbndmodif=upperbnd^vala   
  lowerbndmodif=upperbndmodif*10^(-5) 
  I <- integrate(f=funcLaplace,
                 lower=lowerbndmodif,
                 upper=upperbndmodif,
                 subdivisions=1000L,
                 rel.tol=tol,
                 m=j,lam=lambda,k=k,a=vala,
                 stop.on.error = FALSE)
  resultat <- I$value
  return(resultat)
}

#### Jacobian computation
jacobianFunctiong12 <- function(lam.vec,k.vec,debugg=FALSE){
  
  ## This function computes the J (2x2) jacobian matrices of function
  ## g : (lambda_tj,k_tj) -> 
  ##     ( g1(lambda_tj,k_tj), g2(lambda,k) ) = ( E(G(Z_tj)) , E(G^2(Z_tj)) )
  ##
  ## at the values (lambda_t,k_t) given as inputs.
  ##
  ## Input :
  ## - vectors of same sizes containing values of lambda and k
  ##   (in practice, these values are estimates of lambda and k issued from weibullGMMestim())
  ##
  ## Output :
  ## - a list containing the J 2x2 jacobian matrices described above 
  ##
  ## Used : weibullFsolve() (as J12), matcovtheta_t(), varp1rfar_t()
  ##
  ## Requires : dgjoverdlambdafunc(), dgjoverdkfunc()
  
  lvec <- length(lam.vec)
  listejacobiennes <- list()
  for (i in 1:lvec){
    lambda <- lam.vec[i] ; k <- k.vec[i]
    dg1surdlambda <- dgjoverdlambdafunc(1,lambda,k)
    if (debugg){ cat("i=",i,": dg1dlam,") }
    dg1surdk <- dgjoverdkfunc(1,lambda,k)
    if (debugg){ cat("dg1dk,")}
    dg2surdlambda <- dgjoverdlambdafunc(2,lambda,k)
    if (debugg){ cat("dg2dlam,")}
    dg2surdk <- dgjoverdkfunc(2,lambda,k)
    if (debugg){ cat("dg2dk\n")}
    
    listejacobiennes[[i]] <- matrix(c(dg1surdlambda,dg1surdk,dg2surdlambda,dg2surdk),
                                    2,2,byrow=TRUE)
  }
  return(list(listejacobiennes))
}

foncdgjoverdlambda <- function(u,j,lam,k,a){
  ## Utilitary function used inside dgjoverdlambdafunc() 
  (-j/a^((1/k)+1)) * (-log(u))^(1/k) * exp( -(j*lam/a^(1/k))*(-log(u))^(1/k) ) * u^(1/a - 1) 
} 

dgjoverdlambdafunc <- function(j,lambda,k,lowerbnd=10^(-6),fac=0.5){
  
  ## This utilitary function computes the partial derivative lambda
  ## of the expectation E( exp(-j*W) ) where W is Weibull(lambda,k).
  ##
  ## Used in : jacobianFunctiong12(), jacobianFunctiongrminus1()
  ##
  ## Requires : foncdgjoverdlambda()
  
  vala=fac*(j*lambda)^k
  I <- integrate(f=foncdgjoverdlambda,lower=lowerbnd,upper=1,subdivisions=1000L,
                 j=j,lam=lambda,k=k,a=vala,
                 stop.on.error = FALSE)
  return(I$value) 
}

foncdgjoverdk <- function(u,j,lam,k,a){
  ## Utilitary function used inside dgjoverdkfunc() 
  (-lam/k^2) * log( (1/a)*(-log(u)) ) * foncdgjoverdlambda(u,j,lam,k,a)
} 

dgjoverdkfunc <- function(j,lambda,k,lowerbnd=10^(-6),fac=1){
  
  ## This utilitary function computes the partial derivative k
  ## of the expectation E( exp(-j*W) ) where W is Weibull(lambda,k).
  ##
  ## Used in : jacobianFunctiong12(), jacobianFunctiongrminus1()
  ##
  ## Requires : foncdgjoverdk()
  
  vala=fac*(j*lambda)^k
  I <- integrate(f=foncdgjoverdk,lower=lowerbnd,upper=1,subdivisions=1000L,
                 j=j,lam=lambda,k=k,a=vala,stop.on.error = FALSE)
  return(I$value) 
}

################################################################################################
############################ fsolve method #####################################################

functionp12p13 <- function(theta,vecx){
  
  ## This routine computes both \hat{p}12_t and  \hat{p}13_t 
  ##
  ## Input :
  ## - Counterfactual an factual matrix  matX and matZ , where each column is a trayectory Xi or Zi
  ## - tt, vector tt of length(Z) containing the Z's time steps
  ## - evaluation vector t_eval (in practice t_eval = tt)
  ## - kernel bandwith value h 
  ##
  ## Output : list of 2 elements
  ## - matrix matp12 of \hat{p}12_t , where the column i is associated to {Xi,Zi}
  ## - same for matrix matp13 
  ##
  ## Used in : FastTestforp1r_gumbel 
  ##
  ## Requires : p12_Nonpar, p13_NonPar
  
  lambdaval=theta[1] ; kval = theta[2]
  p12 = vecx[1] ; p13 = vecx[2]
  p12val <- laplaceWeibull(j=1,lambda=lambdaval,k=kval,lowerbnd=10^(-8))
  p13val <- laplaceWeibull(j=2,lambda=lambdaval,k=kval,lowerbnd=10^(-8))
  M <- c( p12 - p12val ,  p13 - p13val )
  return(M)
}

fg <- function(theta){
  ## Utilitary function used inside weibullFsolve()
  ## it allows us to use functionp12p13() just giving it theta as input 
  functionp12p13(theta,vecx=phat)}


J12 <- function(x){
  ## Utilitary function used inside weibullFsolve()
  ## it allows us to use functionp12p13() just giving a vector as imput 
  jacobianFunctiong12(x[1],x[2])}

weibullFsolve <- function(matp12,matp13,truevalues=NULL){
  
  ## This function computes the estimates of lambda_t and k_t 
  ## from the output of the the estimations \hat{p}12, \hat{p}13
  ## solving the equation system with Newton-Raphson method
  ##
  ## Input :
  ## - matrix matp12 and matp13, contaning \hat{p}12 and  \hat{p}13 estimations
  ##   where each column is associated with trayectory Xi or Zi.
  ## - truevalues : an optional matrix, which can contain appropriate 
  ##   starting values for the estimation process. 
  ##   For instance, in a simulation setting, this optional matrix could
  ##   contain the true underlying values of the Weibull parameters lambda and k. 
  ##   (this matrix should contain the lambda values in the first row, and 
  ##    the k values in the second row ; the number of columns must be same as that of matGm)
  ##
  ## Output :
  ## - a list of 2 elements, the first one ($lambdahat) containing the estimate
  ##   of the scale lambda parameter for each column of matp12 and matp13 (ie of matX and matZ),  
  ##   and the second one ($khat) containing the estimate of the shape k shape parameter.
  ##
  ## Requires : the fsolve() routine and the utilitary function functionp12p13() 

  matp12 <- as.matrix(matp12)
  matp13 <- as.matrix(matp13)
  Nligne <- dim(matp12)[1]
  Ncol <- dim(matp12)[2] 
  lambdahat.mat <- matrix(nrow=Nligne,ncol=Ncol)
  khat.mat <- matrix(nrow=Nligne,ncol=Ncol)
  if ( is.null(truevalues) ){
    startvalueslambd <- matrix(1,nrow=Nligne,ncol=Ncol)
    startvaluesk <- matrix(1,nrow=Nligne,ncol=Ncol)
    # starting values, if not precised, are set to 1 and 1 for lambda and k
  } else {
    startvalues <- truevalues
  }
  for (j in 1:Ncol){
    p12hat = matp12[,j]
    p13hat = matp13[,j]
    
    lambdahat.vec = lambdahat.mat[,j]
    khat.vec = khat.mat[,j]
    
    startvalueslambd.vec <- startvalueslambd[,j]
    startvaluesk.vec <- startvaluesk[,j]
    
    for (i in 1:Nligne){
      cat("\n",i,"eme ligne :")
      p12hat_t=p12hat[i]
      p13hat_t=p13hat[i]
      phat=c(p12hat_t,p13hat_t)
      
      startvalueslambd_t=startvalueslambd.vec[i]
      startvaluesk_t=startvaluesk.vec[i]
      startval_t=as.numeric(c(startvalueslambd_t,startvaluesk_t))
      
      fg <- function(theta){functionp12p13(theta,vecx=phat)}
      
      Fsolveweibull <- fsolve(fg,x0=startval_t,J=J12
      )
      
      lambdahat.vec[i] <- Fsolveweibull$x[1]
      khat.vec[i] <- FsolveMweibull$x[2] 
    }
    lambdahat.mat[,j] <- lambdahat.vec
    khat.mat[,j] <- khat.vec
  }   
  return( list("lambdahat"=lambdahat.mat,"khat"=khat.mat) ) 
}

################################################################################################
############################ fsolve  modified method ###########################################

fsolve_modif <- function (f, x0, J = NULL, maxiter = 100, tol = .Machine$double.eps^(0.5), 
          ...) 
{
  ## Utilitary function used inside weibullFsolve_modif()
  ## It is a modified version of the fuction fsolve() that uses broyden_modif()
  ## as numerical method, instead of broyden()
  
  if (!is.numeric(x0)) 
    stop("Argument 'x0' must be a numeric vector.")
  x0 <- c(x0)
  fun <- match.fun(f)
  f <- function(x) fun(x, ...)
  n <- length(x0)
  m <- length(f(x0))
  if (n == 1) 
    stop("Function 'fsolve' not applicable for univariate root finding.")
  if (!is.null(J)) {
    Jun <- match.fun(J)
    J <- function(x) J(x, ...)
  }
  else {
    J <- function(x) jacobian(f, x)
  }
  if (m == n) {
    sol = broyden_modif(f, x0, J0 = J(x0), maxiter = maxiter, 
                  tol = tol)
    xs <- sol$zero
    fs <- f(xs)
  }
  else {
    sol <- gaussNewton(x0, f, Jfun = J, maxiter = maxiter, 
                       tol = tol)
    xs <- sol$xs
    fs <- sol$fs
    if (fs > tol) 
      warning("Minimum appears not to be a zero -- change starting point.")
  }
  return(list(x = xs, fval = fs))
}

broyden_modif<-function (Ffun, x0,J0 = NULL, ..., maxiter = 100, tol = .Machine$double.eps^(1/2)) 
{
  ## Numerical method used in fsolve_modif()
  ## It is a modified version of the fuction broyden() in wich we constrait the parameters
  ## lambda and k from taking negativ values. 
  ## the constraints added are: 
  ## lambda_xnew <- (lambda_x0)/2 when lamba_xnew<0 & k_xnew <- (k_x0)/2 when k_xnew<0
  ##
  ## The paramatric jacobian matrix is not used, instead a numerical one is estimated

  if (!is.numeric(x0)) 
    stop("Argument 'x0' must be a numeric (row or column) vector.")
  fun <- match.fun(Ffun)
  F <- function(x) fun(x, ...)
  y0 <- F(x0)
  if (length(x0) != length(y0)) 
    stop("Function 'F' must be 'square', i.e. from R^n to R^n .")
  if (length(x0) == 1) 
    stop("Function 'F' must not be a univariate function.")
  if (is.null(J0)) {
    A0 <- jacobian(F, x0)
  }
  else {
    A0 <- J0
  }
  B0 <- inv(A0)
  if (any(is.infinite(B0))) 
    B0 <- diag(length(x0))
  xnew <- x0 - B0 %*% y0
  #
  #if (xnew[1]<0){xnew <- x0/2}
  if (xnew[1]<0){xnew[1] <- x0[1]/2}
  if (xnew[2]<0){xnew[2] <- x0[2]/2}
  #
  ynew <- F(xnew)
  k <- 1
  while (k < maxiter) {
    s <- xnew - x0
    d <- ynew - y0
    if (norm(s, "F") < tol || norm(as.matrix(ynew), "F") < 
        tol) 
      break
    B0 <- B0 + (s - B0 %*% d) %*% t(s) %*% B0/c(t(s) %*% 
                                                  B0 %*% d)
    x0 <- xnew
    a <- xnew - B0 %*% ynew
    if(isTRUE(a[1] > 0))
    {xnew <- a}
    else{
      xnew <- x0/2
      #xnew[1]<-x0/2
      #xnew[2] <- xnew[2] - B0[2] %*% ynew[2]
    }
    y0 <- ynew
    ynew <- F(xnew)
    k <- k + 1
  }
  if (k >= maxiter) 
    warning(paste("Not converged: Max number of iterations reached."))
  fnew <- sqrt(sum(ynew^2))
  return(list(zero = c(xnew), fnorm = fnew, niter = k))
}

weibullFsolve_modif <- function(matp12,matp13,truevalues=NULL){
  
  ## Similar to weibullFsolve(). This function computes the estimates of
  ## lambda_t and k_t from the output of the the estimations \hat{p}12, \hat{p}13
  ## solving the equation system with a modified Newton-Raphson method
  ## that prevents lambda and k from taking negative values
  ##
  ## Input :
  ## - matrix matp12 and matp13, contaning \hat{p}12 and  \hat{p}13 estimations
  ##   where each column is associated with trayectory Xi or Zi.
  ## - truevalues : an optional matrix, which can contain appropriate 
  ##   starting values for the estimation process. 
  ##   For instance, in a simulation setting, this optional matrix could
  ##   contain the true underlying values of the Weibull parameters lambda and k. 
  ##   (this matrix should contain the lambda values in the first row, and 
  ##    the k values in the second row ; the number of columns must be same as that of matGm)
  ##
  ## Output :
  ## - a list of 2 elements, the first one ($lambdahat) containing the estimate
  ##   of the scale lambda parameter for each column of matp12 and matp13 (ie of matX and matZ),  
  ##   and the second one ($khat) containing the estimate of the shape k shape parameter.
  ##
  ## Requires : the fsolve_modif() routine and the utilitary function functionp12p13() 
  
  matp12 <- as.matrix(matp12)
  matp13 <- as.matrix(matp13)
  Nligne <- dim(matp12)[1]
  Ncol <- dim(matp12)[2] 
  lambdahat.mat <- matrix(nrow=Nligne,ncol=Ncol)
  khat.mat <- matrix(nrow=Nligne,ncol=Ncol)
  if ( is.null(truevalues) ){
    startvalueslambd <- matrix(1,nrow=Nligne,ncol=Ncol)
    startvaluesk <- matrix(1,nrow=Nligne,ncol=Ncol)
    # starting values, if not precised, are set to 1 and 1 for lambda and k
  } else {
    startvalues <- truevalues
  }
  for (j in 1:Ncol){
    p12hat = matp12[,j]
    p13hat = matp13[,j]
    
    lambdahat.vec = lambdahat.mat[,j]
    khat.vec = khat.mat[,j]
    
    startvalueslambd.vec <- startvalueslambd[,j]
    startvaluesk.vec <- startvaluesk[,j]
    
    for (i in 1:Nligne){
      cat("\n",i,"eme ligne :")
      p12hat_t=p12hat[i]
      p13hat_t=p13hat[i]
      phat=c(p12hat_t,p13hat_t)
      
      startvalueslambd_t=startvalueslambd.vec[i]
      startvaluesk_t=startvaluesk.vec[i]
      startval_t=as.numeric(c(startvalueslambd_t,startvaluesk_t))
      
      Fsolveweibull_modif <- fsolve_modif(functionp12p13,x0=startval_t,vecx=phat
                            )
      
      lambdahat.vec[i] <- Fsolveweibull_modif$x[1]
      khat.vec[i] <- Fsolveweibull_modif$x[2] 
    }
    lambdahat.mat[,j] <- lambdahat.vec
    khat.mat[,j] <- khat.vec
  }   
  return( list("lambdahat"=lambdahat.mat,"khat"=khat.mat) ) 
}

################################################################################################
############################ Optisation approcah ###############################################

#### Using quasi-Newton methods (BFGS and L-BFGS-B)

library(pracma)

fgoptim<- function (theta,vecx){
  ## Utilitary function used in weibullOptim() that allow us to solve the system
  ## as the optimisation of the crossproduct of the equations.
  crossprod(functionp12p13(theta,vecx))
  }

# fg2 <- function (theta){fgoptim(theta,vecx=phat)}
# Gradfg <- t(JX)%*%(functionp12p13(theta,phat)[[1]]) + crossprod(fg,JX) 

weibullOptim<- function(matp12,matp13,truevalues=NULL){
  
  ## This function computes the estimates of lambda_t and k_t 
  ## from the output of the the estimations \hat{p}12, \hat{p}13
  ## solving the equation system as a optimisation problem 
  ##
  ## Input :
  ## - matrix matp12 and matp13, contaning \hat{p}12 and  \hat{p}13 estimations
  ##   where each column is associated with trayectory Xi or Zi.
  ## - truevalues : an optional matrix, which can contain appropriate 
  ##   starting values for the estimation process. 
  ##   For instance, in a simulation setting, this optional matrix could
  ##   contain the true underlying values of the Weibull parameters lambda and k. 
  ##   (this matrix should contain the lambda values in the first row, and 
  ##    the k values in the second row ; the number of columns must be same as that of matGm)
  ##
  ## Output :
  ## - a list of 2 elements, the first one ($lambdahat) containing the estimate
  ##   of the scale lambda parameter for each column of matp12 and matp13 (ie of matX and matZ),  
  ##   and the second one ($khat) containing the estimate of the shape k shape parameter.
  ##
  ## Requires : the optim() routine (we use  L-BFGS-B or BFGS method)
  ##            and optiwith and the utilitary function fgoptim() 
  
  matp12 <- as.matrix(matp12)
  matp13 <- as.matrix(matp13)
  Nligne <- dim(matp12)[1]
  Ncol <- dim(matp12)[2] 
  lambdahat.mat <- matrix(nrow=Nligne,ncol=Ncol)
  khat.mat <- matrix(nrow=Nligne,ncol=Ncol)
  if ( is.null(truevalues) ){
    startvalueslambd <- matrix(1,nrow=Nligne,ncol=Ncol)
    startvaluesk <- matrix(1,nrow=Nligne,ncol=Ncol)
  } else {
    startvalues <- truevalues
  }
  for (j in 1:Ncol){
    p12hat = matp12[,j]
    p13hat = matp13[,j]
    
    lambdahat.vec = lambdahat.mat[,j]
    khat.vec = khat.mat[,j]
    
    startvalueslambd.vec <- startvalueslambd[,j]
    startvaluesk.vec <- startvaluesk[,j]
    
    for (i in 1:Nligne){
      cat("\n",i,"eme ligne :")
      p12hat_t=p12hat[i]
      p13hat_t=p13hat[i]
      phat=c(p12hat_t,p13hat_t)
      
      startvalueslambd_t=startvalueslambd.vec[i]
      startvaluesk_t=startvaluesk.vec[i]
      startval_t=as.numeric(c(startvalueslambd_t,startvaluesk_t))
      Optimweibull <- optim(par=startval_t,fn=fgoptim,method="L-BFGS-B",vecx=phat,lower = 10^(-5)
      )
      #Optimweibull <- optim(par=startval_t,fn=fgoptim,method="BFGS",vecx=phat
      #)
      
      lambdahat.vec[i] <- Optimweibull$par[1]
      khat.vec[i] <- Optimweibull$par[2] 
    }
    lambdahat.mat[,j] <- lambdahat.vec
    khat.mat[,j] <- khat.vec
  }   
  return( list("lambdahat"=lambdahat.mat,"khat"=khat.mat) ) 
}

################################################################################################
############################ Modified GMM méthod ###############################################

library(gmm)

matGZ_func<- function(matX,matZ){
  
  ## This function of computes G(Ztj) for a given trajectory X and Z, where the imput variables are two matrix
  ## matX, matZ, conteining each one multiple trajectories X and Z
  ## This allow us to compute G(Ztj) multiple times for differents trajectories. 
  ## matX and matZ must have the same number of columns
  ##
  ## Input :
  ## - matrices matX and matZ of trajectories X and Z, where each column is a trayectory
  ##
  ## Output :
  ## - matrix of dimentions ( dim(matX)[1] x J ) containong G(Ztj) values 
  ##
  ## Used in : FastTestforp1r_gumbel()
  
  matX <- as.matrix(matX)
  matZ <- as.matrix(matZ) 
  dimnsZ=dim(matZ)
  matGm <- matrix(nrow=dimnsZ[1],ncol=dimnsZ[2])
  for (j in 1:dimnsZ[2]){
    X <- matX[,j]
    Z <- matZ[,j]
    G_empirique<-ecdf(X) 
    matGm[,j] <- G_empirique(Z)
  }
  return(matGm)
}

function_gmm_noyaux<-function(theta,vecx,index,tt.vec,t_eval.vec,bandwidth,kern=dEpan){
  
  ## Utilitary function required for the M-estimation of (lambda_t,k_t) 
  ## in function weibullGMM_NonStationaire() and weibullGMM_NonStationaire_startval_1() 
  ##
  ## Input :
  ## - theta : matrix of dimesion 2 x J containing the tentative values of lambda_t and k_t
  ## - vecx  : a vector of size n supposed to contain the values \hat{G}_m(Z_tj)
  ## - index : value i that tell us the time step tt[i] in wich we are doing the estimation
  ## - tt: vector tt of length(Z) containing the Z's time steps
  ## - t_eval: evaluation vector (in practice t_eval = tt)
  ## - bandwidth: kernel bandwith value 
  ##
  ## Output :
  ## - a J x 2 matrix which first and second columns respectively contain the values 
  ##      (  J * Kij_tj/sum(Kij_tj)*\hat{G}_m(Z_tj) - p_12t(lambda_tj,k_tj)  ) and 
  ##        (  J * Kij_tj/sum(Kij_tj)*(\hat{G}_m(Z_tj)^2) - p_13t(lambda_tj,k_tj)  )  
  ##
  ## Used in : weibullGMM_NonStationaire(),  weibullGMM_NonStationaire_startval_1()
  ## Requires : laplaceWeibull() 
  
  point_eval=t_eval.vec[index]
  lambdaval=theta[1]; kval=theta[2]
  
  p12val <- laplaceWeibull(j=1,lambda=lambdaval,k=kval,lowerbnd=10^(-8))
  p13val <- laplaceWeibull(j=2,lambda=lambdaval,k=kval,lowerbnd=10^(-8))
  
  Kij_ti <- outer(point_eval,tt.vec,function(zz,z)dEpan((zz-z)/bandwidth))
  Kij_ti <- t(Kij_ti)
  W <- Kij_ti/sum(Kij_ti)
  
  p12_nonMoyenne <- W*vecx
  p12_nonMoyenne <- as.vector (p12_nonMoyenne)
  p13_nonMoyenne <- W*(vecx)^2
  p13_nonMoyenne <- as.vector (p13_nonMoyenne)
  
  M <- cbind( length(vecx)*p12_nonMoyenne - p12val ,  length(vecx)*p13_nonMoyenne - p13val )
  return(M)
}

####  Initialisation 1 : Using the previus calculated lambda and k as start values ##################

weibullGMM_NonStationaire <- function(matGm, tt, t_eval, h, kern=dEpan, truevalues=NULL){ # ... also possible
  
  ## This function computes the estimates of multiple lambda_t and k_t from the output of the 
  ## function matGZ_func(), via M-estimation
  ##
  ## Input :
  ## - matGm : a matrix issued from matGZ_func(), containing the values \hat{G}m(Z_tj) , where each column is associated 
  ##   to a different trayectory X_t and Z_t.
  ## - tt: vector tt of length(Z) containing the Z's time steps
  ## - t_eval: evaluation vector (in practice t_eval = tt)
  ## - truevalues : an optional vector, which can contain appropriate starting values for the GMM estimation process of 
  ##   lambda_1 and k_1, i.e, the estimations for the first time step.
  ##   if the starting values ar not given, they will be chosen as [1,1]
  ##   for the next times steps t=2, .... the starting value at t will be the optimal found in t-1.
  ## - kernel bandwith value h 
  ##
  ## Output :a list of 2 elements
  ##   the first one ($lambdahat) containing the GMM estimate of the scale lambda_t parameter for each column of matGm
  ##   (ie of matX and matZ) and for each time step t.  
  ##   the first one ($khat) containing the GMM estimate of the scale k_t parameter for each column of matGm
  ##   (ie of matX and matZ) and for each time step t.  
  ##
  ## Requires : the gmm() routine from the gmm package and the utilitary function function_gmm_noyaux()
  ## 
  ## In gmm(), the "nlminb" optimization routine is used.
  matGm <- as.matrix(matGm)
  Nligne <- dim(matGm)[1] 
  Ncol <- dim(matGm)[2]
  lambdahat.mat <- matrix(nrow=Nligne,ncol=Ncol)
  khat.mat <- matrix(nrow=Nligne,ncol=Ncol)
  
  if ( is.null(truevalues) ){
    startvalueslambd <- 1
    startvaluesk <- 1
  } else {
    startvalueslambd <- truevalues[1] # truevalues is a vector c(1,1)
    startvaluesk <- truevalues[2]
  }
  for (j in 1:Ncol){
    Gnhat = matGm[,j]
    
    lambdahat.vec = lambdahat.mat[,j]
    khat.vec = khat.mat[,j]
    
    startvalueslambd_t <- startvalueslambd 
    startvaluesk_t<- startvaluesk
    
    for (i in 1:Nligne){
      
      startval_t=as.numeric(c(startvalueslambd_t,startvaluesk_t))
      
      fg_noyaux <- function(theta,vecx){function_gmm_noyaux(theta,vecx,index=i,tt.vec=tt,t_eval.vec=t_eval,bandwidth=h)}
      
      EGMMweibull_NonStationary <- gmm(g=fg_noyaux,
                                       x=Gnhat,
                                       t0=startval_t,
                                       optfct="nlminb",
                                       lower=c(10^(-8),10^(-8)),upper=c(Inf,Inf),
                                       onlyCoefficients = TRUE
      )
      
      lambdahat.vec[i] <- EGMMweibull_NonStationary$coefficients[1]
      khat.vec[i] <- EGMMweibull_NonStationary$coefficients[2] 
      
      startvalueslambd_t <- lambdahat.vec[i]
      startvaluesk_t <- khat.vec[i]
      
    }
    lambdahat.mat[,j] <- lambdahat.vec
    khat.mat[,j] <- khat.vec
  }   
  return( list("lambdahat"=lambdahat.mat,"khat"=khat.mat) ) 
  
}

####  Initialisation 2 : start values are the same at each iteration ##################
# If truevalues== NULL the algorthim starts at (1,1).

weibullGMM_NonStationaire_startval_1 <- function(matGm, tt, t_eval, h, kern=dEpan, truevalues=NULL){ # ... also possible
  
  ## This function computes the estimates of multiple lambda_t and k_t from the output of the 
  ## function matGZ_func(), via M-estimation
  ##
  ## Input :
  ## - matGm : a matrix issued from matGZ_func(), containing the values \hat{G}m(Z_tj) , where each column is associated 
  ##   to a different trayectory X_t and Z_t.
  ## - tt: vector tt of length(Z) containing the Z's time steps
  ## - t_eval: evaluation vector (in practice t_eval = tt)
  ## - truevalues : an optional matrix of dimsion 2xJ, which can contain appropriate starting values for the GMM estimation
  ##   process of lambda_t and k_t. It is a matrix so we can chose different starting valueas at each time step t 
  ##   if the starting values ar not given, they will be chosen as [1,1]
  ## - kernel bandwith value h 
  ##
  ## Output :a list of 2 elements
  ##   the first one ($lambdahat) containing the GMM estimate of the scale lambda_t parameter for each column of matGm
  ##   (ie of matX and matZ) and for each time step t.  
  ##   the first one ($khat) containing the GMM estimate of the scale k_t parameter for each column of matGm
  ##   (ie of matX and matZ) and for each time step t.  
  ##
  ## Requires : the gmm() routine from the gmm package and the utilitary function function_gmm_noyaux()
  ## 
  ## In gmm(), the "nlminb" optimization routine is used.
  matGm <- as.matrix(matGm)
  Nligne <- dim(matGm)[1] 
  Ncol <- dim(matGm)[2]
  lambdahat.mat <- matrix(nrow=Nligne,ncol=Ncol)
  khat.mat <- matrix(nrow=Nligne,ncol=Ncol)
  
  if ( is.null(truevalues) ){
    startvalueslambd <- matrix(1,nrow=Nligne,ncol=Ncol)
    startvaluesk <- matrix(1,nrow=Nligne,ncol=Ncol)
  } else {
    startvalueslambd <- truevalues[[1]] # truevalues doit être une liste avec deux matrices nxN
    startvaluesk <- truevalues[[2]]
  }
  for (j in 1:Ncol){
    Gnhat = matGm[,j]
    
    lambdahat.vec = lambdahat.mat[,j]
    khat.vec = khat.mat[,j]
    
    startvalueslambd.vec <- startvalueslambd[,j]
    startvaluesk.vec <- startvaluesk[,j]
    
    for (i in 1:Nligne){
      
      startvalueslambd_t=startvalueslambd.vec[i]
      startvaluesk_t=startvaluesk.vec[i]
      startval_t=as.numeric(c(startvalueslambd_t,startvaluesk_t))
      
      fg_noyaux <- function(theta,vecx){function_gmm_noyaux(theta,vecx,index=i,tt.vec=tt,t_eval.vec=t_eval,bandwidth=h)}
      
      EGMMweibull_NonStationary <- gmm(g=fg_noyaux,
                                       x=Gnhat,
                                       t0=startval_t,
                                       optfct="nlminb",
                                       lower=c(10^(-8),10^(-8)),upper=c(Inf,Inf),
                                       onlyCoefficients = TRUE
      )
      
      lambdahat.vec[i] <- EGMMweibull_NonStationary$coefficients[1]
      khat.vec[i] <- EGMMweibull_NonStationary$coefficients[2] 
    }
    lambdahat.mat[,j] <- lambdahat.vec
    khat.mat[,j] <- khat.vec
  }   
  return( list("lambdahat"=lambdahat.mat,"khat"=khat.mat) ) 
  
}

################################################################################################
############################ p1r_t and far_t computation #######################################
################################################################################################

p1rfarW_temps<- function(lam.mat,k.mat,r.mat,lowerbnd=10^(-5)){
  
  ## This function computes the parametric estimates of p_1r.t and far(r).t from each column of 
  ## the matrices of parametes, i.e, lambda_t , k_t  and r_t (in practice we will compute p1r.t 
  ## for the same r at each time step, so we will use r_t=rep(r,) )
  ##
  ## Accepting matrices as imput allow us to repeat the procedure multiple times. computtaions 
  ## for each column are independent.
  ##
  ## Input :
  ## - matrices of the same dimentions containing values of lambda_t, k_t, r
  ##   (in practice, these values are the estimated lambda_t and k_t issued from weibullGMM_NonStationaire()
  ##   or weibullGMM_NonStationaire_startval_1()
  ## - an optional technical parameter, lowerbnd (used in laplaceWeibull() )
  ## 
  ## Output :
  ## - if the values (lambda_t,k_t) contain estimates of lambda_t and k_t, 
  ##   then the ouptut is the list of associated couples of values ( \hat{p}^(W)_1r.t , \hat{far}^(W)(r).t ) 
  ## - if the values (lambda_t,k_t) contain true values of lambda and k
  ##   (theoretical), then the output is the list of associated couples of values ( p_1r.t , far(r).t ) 
  ##
  ## Used in: FastTestforp1r_gumbel(), matcovp12p13_t(), varp1rfar_t(), CI_p1rfar()
  ##
  ## Requires : laplaceWeibull() 
  
  lam.mat<-as.matrix(lam.mat)
  k.mat<-as.matrix(k.mat)
  n = dim(lam.mat)[1]
  N = dim(lam.mat)[2]
  p1r.mat = matrix(0,n,N)
  for (j in 1:N){#  for each trajectory
    lam.vec = lam.mat[,j]
    k.vec = k.mat[,j]
    r.vec = r.mat[,j]
    p1r.vec = p1r.mat[,j]
    for (i in 1:n){ # for each t in the trajectory
      p1r_t <- laplaceWeibull(r.vec[i]-1,lam.vec[i],k.vec[i],lowerbnd=lowerbnd)
      p1r.vec[i] = p1r_t
    }
    p1r.mat[,j] <- p1r.vec
  }  
  far.mat <- 1 - 1 / (r.mat * p1r.mat)
  return( list("p1r"=p1r.mat,"far"=far.mat) )
}


#### Function that applies optim/gmm method to multiple trajectories {(X)t,(Z)t}
#### for Gumbel trayectories, where Xt is stationary and Zt oth linear trend on muz

FastTestforp1r_gumbel <- function(tt,h,r,I,J,N,sigX.vec,sigZ.vec,muX.vec,muZ.vec,methode="optim"){
  
  ## This function creates multiple Gumbel trayectories {X_t,Z_t}  and for a given r , 
  ## and estimes its p1r_t  using optimisation or M-moments(with learned starting values) method.
  ## 
  ## All the trayectories X_t created come from the same density distribution, same for Z_t
  ##
  ## Input :
  ## - tt, vector tt of length(Z) containing the Z's time steps
  ## - kernel bandwith value h 
  ## - desired value of r for the study 
  ## - I: length of trajextory X_t
  ## - J: length of trajectory Z_t
  ## - N: numer od trajectories created
  ## - GEV's paramteres: 
  ##   sigX: scale parameter for X_t trajectory
  ##   sigZ.vec: scale parameter for Z_t trajectory, it is a vector because Z_t can be non-stationnary 
  ##   muX.vec: position parameter of X_t, it's a vector so it can be chosen to be variable in time
  ##   muZ.vec: position parameter of Z_t, it's a vector so it can be chosen to be variable in time
  ##   method: paramter to chose the estimation method, it can be "gmm" or "optim"
  ## 
  ## Output :
  ## - list of two objects, a matrix ($matp1r) containing the the p1r_t estimated trajectories,
  ##   and a vector ($p1r_mean) mean betwwen the different trajectories.
  ## 
  ## Used in: nowhere, just a function that allow us to test the different methods 
  ##
  ## Requires : P12_P13_estimation(), weibullOptim(), matGZ_func(), weibullGMM_NonStationaire(), p1rfarW_temps()
  
  theta_theo <- 1 / exp(muZ.vec)
  p1r_t_theo <- 1 / (1 + (r-1)*theta_theo)
  matX <- matrix(nrow=I,ncol=N)
  matZ <- matrix(nrow=J,ncol=N)
  matp <- matrix(nrow=I,ncol=N)
  for (j in 1:N){
    matX[,j] <- rgev(I, loc = muX.vec, scale = sigX.vec, shape = 0)
    matZ[,j] <- rgev(J, loc = muZ.vec, scale = sigZ.vec, shape = 0)
  }
  if (methode=="optim"){
    matp<-P12_P13_estimation(matX,matZ,tt,tt,h,kern=dEpan)
    thetahat<-weibullOptim(matp$matp12,matp$matp13,truevalues=NULL)
  }
  if (methode=="gmm"){
    GZ<-matGZ_func(matX,matZ)
    thetahat<-weibullGMM_NonStationaire(GZ, tt, tt, h, kern=dEpan, truevalues=NULL)
  }
  p1rfar<-p1rfarW_temps(thetahat[[1]],thetahat[[2]],matrix(r,ncol=N,nrow=J))
  p1r_mean <- rowMeans(p1rfar$p1r)
  plot(tt,p1r_t_theo,type="l",col="black",xlab="time",ylab=expression(p["1,10,t"]), main=expression(Evolution ~  over  ~ time  ~ of ~ p["1,10"]))
  lines(tt,p1r_mean,col="red")
  #for (i in 1:N){
   # lines(tt,p1rfar$p1r[,i],col="gray")
  #}
  return(list("matp1r"=p1rfar$p1r,"p1r_mean"=p1r_mean))
}

################################################################################################
############################ Trajectory simulations ############################################
################################################################################################

#### Simulate W-class trajectories with non-sationnary Z 

simulWclass_nonStationnary_general <-function(I,N,ksiX,sigX=1,muX=0,ksiZ=NULL,muZ=NULL,sigZ=NULL,unknown="location", graph=TRUE){
  
  ## This function simulates X_t and Z_t trayectories fulfilling Wclass condition. all of X_t's parameters must be know as well
  ## as two from three of Z_t parameters, the third parameter will be computated in order to fufill the Wclass condition.
  ## The function allow us to chose the unknow parameter and give us the possibility of visualizing the distribution density at 
  ## different times. 
  ##
  ## X_t length must be given directly  as a parameter of the function and Z_t length will be taken from the length of its shape 
  ## parameter ksiZ. If ksi is unknown,the length of sigmaZ will be the one taken in consideration.
  ## it's paramrataers inputs.
  ##
  ## Input:
  ## - m: X_t's length
  ## - n: the number of Xt and Zt trajectories that we create
  ## - ksix: X_t's shape parameter
  ## - sigX: X_t's scale parameter. if not given, its default value is  1
  ## - muX: X_t's position parameter. if not given, its default value is 0
  ## - ksiZ: Z_t's shape parameter
  ## - muZ: Z_t's position parameter
  ## - sigZ: Z_t's scale parameter
  ## - unknown: "location", "scale" or "shape"
  ## - graph: TRUE or FALSE
  ##
  ## Output :
  ## - list of multiple objects. matrix matZ, vector lambda_t, vector k_t, vector muZ, vector sigZ,
  ##   vector ksiZ, W ( W[,j] <- -log(Gm(matZ[,j])) ) 
  ## 
  ## Used in: nowhere , it is used to create tests
  ##
  ## Requires : just gev package
  
  
  if (ksiX!=0){
    if (unknown == "location"){
      # We take our Z's sample size from the length of shape input
      J <- length(ksiZ)
      # We remind that if our shape input is just a value, the sample size will be 1
      if (J == 1){
        warning("Z's shape parameters length give us the size of this strajectory")
      }
      # Verification. Shape parameters must be of the same sign
      for (i in 1:J){
        if (sign(ksiZ[i])!=sign(ksiX)){
          stop("X and Z can not be W-class, shape parameters must have the same sign")
        }
      }
      # We must verify that only Z's location parameter is unknown
      if (!is.null(muZ)){
        stop("Location parameter of Z must be unknown")
      }
      
      if (is.null(sigZ)){
        stop("Scale parameter of Z must be known")
      }
      if (is.null(ksiZ)){
        stop("Shape parameter of Z must be known")
      }
      # Verification that X's parameters are just values
      if (length(muX)!=1 | length(sigX)!=1 | length(ksiX)!=1){
        stop("X's parameters must be of just 1 value")
      }
      # Z's parameters must be of the same length
      if (length(sigZ) == 1) {
        sigZ <- rep(sigZ,J)
      }
      if (length(sigZ)!=J){
        stop("Z's parameters are not from the same length")
      }
      # support value
      support <- muX - sigX/ksiX
      # We impose the condition of equal support
      muZ <- support + sigZ/ksiZ
    }
    if (unknown == "scale"){
      # We take our sample size from the length of shape inputs
      J <- length(ksiZ)
      # We remind that if our shape inputs are just a value, the sample size will be 1
      if (J == 1){
        warning("Z's shape parameters length will give us the size of this strajectory")
      }
      # Verification. Shape parameters must be of the same sign
      for (i in 1:J){
        if (sign(ksiZ[i])!=sign(ksiX)){
          stop("X and Z can not be W-class, shape parameters must have the same sign")
        }
      }
      # We must verify that only Z's scale parameter is unknown
      if (!is.null(sigZ)){
        stop("Scale parameter of Z must be unknown")
      }
      if (is.null(muZ)){
        stop("Location parameter of Z must be known")
      }
      if (is.null(ksiZ)){
        stop("Shape parameter of Z must be known")
      }
      # Verification that X's parameters are just values
      if (length(muX)!=1 | length(sigX)!=1 | length(ksiX)!=1){
        stop("X's parameters must be of just 1 value")
      }
      # Z's parameters must be of the same length
      if (length(muZ) == 1) {
        muZ <- rep(muZ,J)
      }
      if (length(muZ)!=J){
        stop("Z's parameters are not from the same length")
      }
      # support value
      support <- muX - sigX/ksiX
      # We impose the condition of equal support
      sigZ <- (support - muZ)*(-ksiZ)
    }
    if (unknown == "shape"){
      # We take our sample size from the length of scale inputs
      J <- length(sigZ)
      # We remind that if our shape inputs are just a value, the sample size will be 1
      if (J == 1){
        warning("Z's scale parameters length will give us the size of this strajectory")
      }
      # We must verify that only Z's location parameter is unknown
      if (!is.null(ksiZ)){
        stop("Shape parameter of Z must be unknown")
      }
      if (is.null(muZ)){
        stop("Location parameter of Z must be known")
      }
      if (is.null(sigZ)){
        stop("Scale parameter of Z must be known")
      }
      # Verification that X's parameters are just values
      if (length(muX)!=1 | length(sigX)!=1 | length(ksiX)!=1){
        stop("X's parameters must be of just 1 value")
      }
      # Z's parameters must be of the same length
      if (length(muZ) == 1) {
        muZ <- rep(muZ,J)
      }
      if (length(muZ)!=J){
        stop("Z's parameters are not from the same length")
      }
      
      # support value
      support <- muX - sigX/ksiX
      # We impose the condition of equal support
      ksiZ <- sigZ/(muZ - support)
      # Verification. Shape parameters must be of the same sign
    }
    
      for (i in 1:J){
        ksiZ_i<- as.numeric(ksiZ[i])
        if ((sign(ksiZ_i)==sign(ksiX))==FALSE){
          stop("X and Z can not be W-class, shape parameters must have the same sign")
        }
      }
      
    matX<-matrix(nrow=I,ncol=N)
    matZ<-matrix(nrow=J,ncol=N)
    
    for(j in 1:N){
      matX[,j] <- rgev(I, loc=muX, scale= sigX, shape=ksiX)
      for(i in 1:J)
        matZ[i,j] <-rgev(1, loc=muZ[i], scale=sigZ[i], shape=ksiZ[i])}
    
    lam_t <- ((sigZ/ksiZ)*(ksiX/ksiZ))^(-1/ksiX)
    k_t <- ksiX/ksiZ
  }
  
  if (ksiX==0){
    J <- length(sigZ)
    # We remind that if our shape inputs are just a value, the sample size will be 1
    # matX<-matrix(nrow=m,ncol=N)
    # matZ<-matrix(nrow=n,ncol=N)
    if (J == 1){
      warning("Z's scale parameters length will give us the size of this strajectory")
      }
    # In this case we must known all the parameters
    if (is.null(sigZ)){
      stop("For Gumbel trajectories all Z parameters must be unknown")
      }
    if (is.null(muZ)){
      stop("For Gumbel trajectories all Z parameters must be unknown")
      }
    if (is.null(ksiZ)){
      stop("For Gumbel trajectories all Z parameters must be unknown")
      }
    # ksiZ will be re defined as 0 as precaution
    ksiZ<-rep(0,J)
  
    # Verification that X's parameters are just values
    if (length(muX)!=1 | length(sigX)!=1 | length(ksiX)!=1){
      stop("X's parameters must be of just 1 value")
      }
    # Z's parameters must be of the same length
    if (length(muZ) == 1) {
      muZ <- rep(muZ,J)
      }
    if (length(sigZ) == 1) {
      sigZ <- rep(sigZ,J)
      }
    if (length(muZ)!=J){
      stop("Z's parameters are not from the same length")
      }
    if (length(sigZ)!=J){
      stop("Z's parameters are not from the same length")
    }
    
    k_t <- sigX/sigZ
    lam_t <- exp((muX-muZ)/sigX)
    
  }
  
  matX<-matrix(nrow=I,ncol=N)
  matZ<-matrix(nrow=J,ncol=N)
  
  for(j in 1:N){
    matX[,j] <- rgev(I, loc=muX, scale= sigX, shape=ksiX)
    for(i in 1:J){
      matZ[i,j] <-rgev(1, loc=muZ[i], scale=sigZ[i], shape=ksiZ[i])}
    
  }
  if (graph==TRUE){
    
    J <-length(ksiZ)

    muZ_t <- c(muZ[1],muZ[J/4],muZ[J/2], muZ[3*J/4], muZ[J])
    sigZ_t <- c(sigZ[1],sigZ[J/4],sigZ[J/2], sigZ[3*J/4],sigZ[J])
    ksiZ_t <- c(ksiZ[1],ksiZ[J/4],ksiZ[J/2], ksiZ[3*J/4], ksiZ[J])

    plot(density(matX[,1]),col="red",main= expression(paste(X[t]," and ",Z[t]," distribution evolution")))

    Z1 = rgev(J, loc = muZ_t[1], scale = sigZ_t[1], shape = ksiZ_t[1])
    lines(density(Z1),col="gray")
    Z2 = rgev(J, loc = muZ_t[2], scale = sigZ_t[2], shape = ksiZ_t[2])
    lines(density(Z2),col="green")
    Z3 = rgev(J, loc = muZ_t[3], scale = sigZ_t[3], shape = ksiZ_t[3])
    lines(density(Z3),col="blue")
    Z4 = rgev(J, loc = muZ_t[4], scale = sigZ_t[4], shape = ksiZ_t[4])
    lines(density(Z4),col="purple")
    legend(26.0,1.5,legend=c(expression(X[t]),expression(Z[t=1]),expression(Z[t=50]),expression(Z[t=100]),expression(Z[t=200])),col=c("red","gray","green","blue","purple"),lty=1:2,cex=0.7)
    
  }
  
  W=matrix(nrow=J, ncol=N)
  for(j in 1:N){
      Gm <- ecdf(matX[,j])
      W[,j] <- -log(Gm(matZ[,j]))
      
    }
  # Calculation of lambda an k
  return (list("matX"=matX,
               "matZ"=matZ,
               "lam"=lam_t,
               "k"=k_t,
               "mu_Z"=muZ,
               "sig_Z"=sigZ,
               "ksi_Z"=ksiZ,
               "W"=W))
  }

############################ Less general fonctions  ############################################

# plot for two Gumbel, X stationary, Z non-sationary with linear growth in location parameter)
plotd_time <- function(size,mux,muz,sigx,sigz){
  
  ## This function plots d'evolution in time of the density distributions of two gumbel trajectories X_t, Z_t
  ## where for instance Z_t is non-stationnary and X_t stationnaty, in order to visualise de resemblance beetwen
  ## the trajectories and the posibility of fulfilling the Wclass condition
  ##
  ## Input :
  ## - size: legth of the trajectories, both have the same length
  ## - GEV's paramteres: 
  ##   sigX: scale parameter value for X_t trajectory
  ##   sigZ.vec: scale parameter vector for Z_t trajectory
  ##   muX: position parameter value of X_t
  ##   muZ: position parameter vector of Z_t
  ##
  ## Output :
  ## - plot of the density distributions
  ## 
  ## Used in: nowhere, just a function that allow to visualise gumbel trajectories
  
  X = rgev(size, loc = mux, scale = sigx, shape = 0)
  muz_t= c(muz[1],mu[length(muz)/4],mu[length(muz)/2], mu[3*length(muz)/4], mu[length(muz)])
  plot(density(X),col="black",)
  for (i in 1:length(muz_t)){
    Z = rgev(size, loc = muz_t[i], scale = sigz, shape = 0)
    lines(density(Z),col="gray")
  }
}

################################################################################################
################################### CI for lambda_t and k_t ####################################
################################################################################################

##################### Variance - covariamce matrix of {p12_t, p13_t} ###########################
#####################                                                ###########################

#### VarCovar matrix ###########################################################################

matcovp12p13_t <- function(X.vec, Z.vec,tt,t_eval,h){
  
  ## This routine computes the variance-covariance matrix of  
  ##
  ##    N = ( (\hat{p}_12.t,\hat{p}_13.t) - (p_12.t,p_13.t) )
  ##
  ## for each time step t 
  ##
  ## Input :
  ## - vectors of size I and J containing the trajectories X and Z
  ## - vector tt of size J containing the Z's time steps
  ## - evaluation vector t_eval (in practice t_eval = tt)
  ## - kernel bandwith h 
  ##
  ## Output :
  ## - Array of length J containg 2x2 asymptotic covariance matrices 
  ##
  ## Used in : ... nowhere yet
  ##
  ## Requires : matcovNA_t(), matcovNB_t()
  
  J <- length(Z.vec)
  I <- length(X.vec)
  
  G_emp <- ecdf(X.vec)
  GmZ <- G_emp(Z.vec)
  
  cat("start lambda_t and k_t computation", "\n")
  theta_t <-weibullGMM_NonStationaire (GmZ, tt, t_eval, h, kern=dEpan, truevalues=NULL)
  lambda_t <- theta_t[[1]]
  k_t <- theta_t[[2]]
  cat("end lambda_t and k_t computation", "\n")
  
  cat("start p12_hat.vec and p13_hat.vec computation", "\n")
  p12_hat.vec <- p12_NonPar(X.vec,Z.vec,tt,t_eval,h)$p12.mat
  p13_hat.vec<- p13_NonPar(X.vec,Z.vec,tt,t_eval,h)
  cat("end p12_hat.vec and p13_hat.vec computation", "\n")
  
  list_weighed_matcov_N_t<- array(NA,c(2,2,J))
  
  list_unweighed_matcov_NB <- (1/I) * 2* matcovNB(p12_hat.vec, p13_hat.vec,lambda_t,k_t)
  list_unweighed_matcov_NA_t<- 2 * matcovNA_t(p12_hat.vec, p13_hat.vec,lambda_t,k_t,J)

  A1ji <- list_unweighed_matcov_NB[1,1,]
  B1ji <- list_unweighed_matcov_NB[2,2,]
  D1ji <- list_unweighed_matcov_NB[1,2,]
  C1ji <- list_unweighed_matcov_NB[2,1,]
  
  A2ji <- list_unweighed_matcov_NA_t[1,1,]
  B2ji <- list_unweighed_matcov_NA_t[2,2,]
  D2ji <- list_unweighed_matcov_NA_t[1,2,]
  C2ji <- list_unweighed_matcov_NA_t[2,1,]
  
  Kh <- outer(t_eval, tt, function(zz,z) dEpan((zz - z) / h))
  
  for (index_t in 1:J){
    
    denom<-(sum(Kh[index_t,]))^2
    
    list_Wji_t <- c()
    list_khj_t<-c()
    
    for (j in 1:dim(Kh)[2]){
      Khj<-Kh[index_t,j]
      
      list_khj_t <- c(list_khj_t, Khj^2)
 
      for (i in 1:j){
        Khi <- Kh[index_t,i]
        
        KhjKhi <- Khj * Khi
        
        list_Wji_t <- c(list_Wji_t, KhjKhi)
      }
    }
    list_weighed_matcov_N_t[1,1,index_t] <- 1/(J^2) *( sum(list_Wji_t * A1ji) + sum(list_khj_t*A2ji)) / denom
    list_weighed_matcov_N_t[1,2,index_t] <- 1/(J^2) *( sum(list_Wji_t * D1ji) + sum(list_khj_t*D2ji) ) /denom
    list_weighed_matcov_N_t[2,1,index_t] <- 1/(J^2) *( sum(list_Wji_t * C1ji) + sum(list_khj_t*C2ji) ) /denom
    list_weighed_matcov_N_t[2,2,index_t] <- 1/(J^2) *( sum(list_Wji_t * B1ji) + sum(list_khj_t*B2ji)) /denom
    cat("time", index_t, "weighed","\n")
  }
  return(list_weighed_matcov_N_t)
}

##### varcovarNB #############################################################################

matcovNB <- function(p12_hat_t, p13_hat_t,lambda_t,k_t){
  
  ## This routine computes the unweighed variance-covariance matrix of NA
  ##
  ##    matcov_NB_unweighed (t) =  ( A1_jit & C1_ijt \\ C1_jit & B1_jit )
  ##
  ## for each time step t 
  ##
  ## Input :
  ## - \hat{p}_12.t and \hat{p}_13.t vectors from kernel estimation
  ## - two vectors of length J with Weibull's estimated parameters \hat{lambda}_t,\hat{k}_t
  ## - J: length of trajectory Z
  ##
  ## Output :
  ## - Array of length J containg 2x2 unweighed covariance matrice of NB 
  ##
  ## Used in : matcovp12p13_t
  ##
  ## Requires : matcovNB_aij_t(), matcovNB_bji(), matcovNB_cji()
  
  A1_ji <- matcovNB_aji(p12_hat_t,lambda_t,k_t)
  B1_ji <- matcovNB_bji(p13_hat_t,lambda_t,k_t)
  C1_ji <- matcovNB_cji(p12_hat_t,p13_hat_t,lambda_t,k_t)
  
  Ncombinations<-length(A1_ji)
  
  list_unweighed_matcov_NB <- array(NA,c(2,2,Ncombinations))
  
  list_unweighed_matcov_NB[1,1,] <- A1_ji
  list_unweighed_matcov_NB[1,2,] <- C1_ji 
  list_unweighed_matcov_NB[2,1,] <- C1_ji 
  list_unweighed_matcov_NB[2,2,] <- B1_ji 
  
  return (list_unweighed_matcov_NB) 
  }

matcovNB_aji <- function(p12_hat_t,lambda_t,k_t){
  
  ## This function computes the unweighted term Aij of NB_n's asymptotic variance-covariance matrix
  ##
  ##    Aij = E( min ( G(Z_tj),G(Z_ti) ) - G(Z_tj )* G(Z_ti) ) =  E ( M_2ji + ( p12_hat.j * p12_hat.i ) )
  ##
  ## Input :
  ## - \hat{p}_12.t and \hat{p}_13.t vectors from kernel estimation
  ## - two vectors of length J containing Weibull's estimated parameters \hat{lambda}_t,\hat{k}_t
  ##
  ## Output :
  ## - Vector containing A1_ji terms
  ##
  ## Used in : matcovNB
  ##
  ## Requires : Mrfuncji()
  
  cat("start: computation of A1ji","\n")
  
  list_unweighed_matcov_NB_aji <- c()
  
  for (j in 1:length(p12_hat_t)){
    cat("start: components for j= ",j,"\n")
    lambda.j <- lambda_t[j]
    k.j <- k_t[j]
    p12_hat.j <- p12_hat_t[j]
    
    for (i in 1:j){
      cat("start: components for i= ",i,"\n")
      lambda.i <- lambda_t[i]
      k.i <- k_t[i]
      p12_hat.i <- p12_hat_t[i]
      
      matcov_NB_aji <-Mrfuncji(2,lambda.j,k.j,lambda.i,k.i) - ( p12_hat.j * p12_hat.i)
      
      list_unweighed_matcov_NB_aji <- c(list_unweighed_matcov_NB_aji, matcov_NB_aji)
      cat("done: unweighed list of components for ",j,", ",i,"\n")
    }
  }
  return (list_unweighed_matcov_NB_aji)
}

matcovNB_cji <- function(p12_hat_t,p13_hat_t,lambda_t,k_t){

  ## This fucntion computes the unweighted term Cij of NB_n's asymptotic variance-covariance matrix
  ##
  ##    Cji =  E( ( 2* G(Z_tj) * min(G(Z_tj),G(Z_ti)) ) - G(Z_tj)^2 *G(Z_ti) )  =  2* E( E_ji - p13_hat.j  * p12_hat.i )
  ##
  ## Input :
  ## - \hat{p}_12.t and \hat{p}_13.t vectors from kernel estimation
  ## - two vectors of length J containing Weibull's estimated parameters \hat{lambda}_t,\hat{k}_t
  ##
  ## Output :
  ## - Vector containing C1_ji terms
  ##
  ## Used in : matcovNB
  ##
  ## Requires : calculEGzjminGzjGzi_partieA1() calculEGzjminGzjGzi_partieB()
  
  cat("start: computation of C1ji","\n")
  
  list_unweighed_matcov_NB_cji <- c()
  
  for (j in 1:length(p12_hat_t)){
    cat("start: components for j= ",j,"\n")
    lambda.j <- lambda_t[j]
    k.j <- k_t[j]
    p12_hat.j <- p12_hat_t[j]
    p13_hat.j <- p13_hat_t[j]
    
    EGzjminGzjGzi_partieA1 <-calculEGzjminGziGzj_partieA1(lambda.j,k.j,lowerbnd=10^(-6),fac=0.5)
  
    for (i in 1:j){
      cat("start: components for i= ",i,"\n")
      lambda.i <- lambda_t[i]
      k.i <- k_t[i]
      p12_hat.i <- p12_hat_t[i]
      p13_hat.i <- p13_hat_t[i]
      
      EGzjminGzjGzi_partieA1bis <-calculEGzjminGziGzj_partieA1(lambda.i,k.i,lowerbnd=10^(-6),fac=0.5)#
      
      EGzjminGzjGzi_partieB <- calculEGzjminGzjGzi_partieB(lambda.j,k.j,lambda.i,k.i,lowerbnd=10^(-5),fac=0.5,tol=10^(-5))
      
      # Adding all parts
      EGzjminGzjGzi <- EGzjminGzjGzi_partieB + p13_hat.j - EGzjminGzjGzi_partieA1
      
      EGzjminGzjGzibis <- EGzjminGzjGzi_partieB + p13_hat.i - EGzjminGzjGzi_partieA1bis
      
      matcov_NB_cji <- 2*(EGzjminGzjGzi+EGzjminGzjGzibis) - 2*((p13_hat.j * p12_hat.i))
    
      list_unweighed_matcov_NB_cji <- c(list_unweighed_matcov_NB_cji, matcov_NB_cji)
      cat("done: unweighed list of components for ",j,", ",i,"\n")
    }
  }
  return (list_unweighed_matcov_NB_cji)
}

matcovNB_bji <- function(p13_hat_t,lambda_t,k_t){

  ## This fucntion computes the unweighted term Cij of NB_n's asymptotic variance-covariance matrix
  ##
  ##    Bij = 4* E( (G(Z_tj) * G(Z_ti) * [ min( G(Z_tj),G(Z_ti) ) - G(Z_tj) * G(Z_ti) ] ) = 4( M_3ji - p13_hat.j * p13_hat.i )
  ##
  ## Input :
  ## - \hat{p}_12.t and \hat{p}_13.t vectors from kernel estimation
  ## - two vectors of length t containing Weibull's estimated parameters \hat{lambda}_t,\hat{k}_t
  ##
  ## Output :
  ## - Vector containing B1_ji terms
  ##
  ## Used in : matcovNB
  ##
  ## Requires : Mrfuncji()
  
  cat("start: computation of B1ji","\n")
  
  list_unweighed_matcov_NB_bji <- c()
  
  for (j in 1:length(p13_hat_t)){
    cat("start: components for j= ",j,"\n")
    lambda.j <- lambda_t[j]
    k.j <- k_t[j]
    p13_hat.j <- p13_hat_t[j]
    
    for (i in 1:j){
      cat("start: components for i= ",i,"\n")
      lambda.i <- lambda_t[i]
      k.i <- k_t[i]
      p13_hat.i <- p13_hat_t[i]
      
      NB_bji <- 4*(Mrfuncji(3,lambda.j,k.j,lambda.i,k.i) - (p13_hat.j * p13_hat.i))
      
      list_unweighed_matcov_NB_bji <- c(list_unweighed_matcov_NB_bji, NB_bji)
      cat("done: unweighed list for j= ", j," i= ",i,"\n")
    }
  }
  return (list_unweighed_matcov_NB_bji)
}


funcLaplaceterji <- function(v,m,lam.j,k.j,lam.i,k.i,a.j,a.i,lowerbnd=10^(-5),tol=10^(-5)){
  # Utilitary function inside Mrfuncji
  nv=length(v)
  mprime.j = (m-1)*lam.j / a.j^(1/k.j)
  mprime.i = (m-1)*lam.i / a.i^(1/k.i)
  facteur1.j = exp( -mprime.j*(-log(v))^(1/k.j) ) 
  facteur1.i = exp( -mprime.i*(-log(v))^(1/k.i) )
  facteur2.j=rep(0,nv)
  facteur2.i=rep(0,nv)
  for (i in 1:nv){
    facteur2.j[i] = laplaceWeibull( j=m, lambda=lam.j, k=k.j, lowerbnd=(10^-6)*v[i],
                                    upperbnd= v[i]^(1/a.i), tol=tol )  # borninf is guaranteed to be both close to 0, AND < v[i] 
  }
  for (i in 1:nv){
    facteur2.i[i] = laplaceWeibull( j=m, lambda=lam.i, k=k.i, lowerbnd=(10^-6)*v[i],
                                    upperbnd= v[i]^(1/a.j), tol=tol )  # borninf is guaranteed to be both close to 0, AND < v[i] 
  }
  facteur3.j=v^(1/a.j - 1)
  facteur3.i=v^(1/a.i - 1)
  return((1/a.j)*facteur1.j*facteur2.i*facteur3.j + (1/a.i)*facteur1.i*facteur2.j*facteur3.i)
}

Mrfuncji <- function(r,lambda.j,k.j,lambda.i,k.i,lowerbnd=10^(-5),fac=0.5,tol=10^(-5)){
  
  ## This function computes the value 
  ##   
  ##    M_rji = E( G(Z_tj)^(r-2)*G(Z_ti)^(r-2) * min(G(Z_tj),G(Z_ti)) )
  ##
  ## which is part of Aji and Bji, component of the covariance matrix of NB
  ## 
  ## Inputs are single values of r, lambda_tj, k_tj , lambda_ti, k_ti
  ##
  ## Used in : matcovNB_aji, matcovNB_bji
  ##
  ## Requires : funcLaplaceterji() 
  
  vala.j=fac*((r-1)*lambda.j)^k.j
  vala.i=fac*((r-1)*lambda.i)^k.i
  I <- integrate(f=funcLaplaceterji,
                 lower=lowerbnd,upper=1,
                 subdivisions=1000L,
                 rel.tol=10^(-5),
                 m=r-1,
                 lam.j=lambda.j,
                 k.j=k.j,
                 lam.i=lambda.i,
                 k.i=k.i,
                 a.j=vala.j,
                 a.i=vala.i,
                 lowerbnd=lowerbnd,
                 stop.on.error = FALSE)
  return(I$value)
}


funcEGzjminGzjGzi_partieB <- function(v,lam.j,k.j,lam.i,k.i,a.j,a.i,lowerbnd=10^(-5),tol=10^(-5)){
  
  ## Utilitary function used inside function calculEGzjminGzjGzi_partieB()
  ##
  ##
  ## Requires : laplaceWeibull()
  
  nv=length(v)
  mprime.j = lam.j / a.j^(1/k.j)
  facteur1.j = exp( -mprime.j*(-log(v))^(1/k.j) ) 
  facteur2.i=rep(0,nv)
  for (i in 1:nv){
    facteur2.i[i] = laplaceWeibull( j=1, lambda=lam.i, k=k.i, lowerbnd=(10^-6)*v[i],
                                    upperbnd= v[i]^(1/a.j), tol=tol )  # borninf is guaranteed to be both close to 0, AND < v[i] 
  }
  facteur3.j=v^(1/a.j - 1)
  return((1/a.j)*facteur1.j*facteur2.i*facteur3.j)
}

calculEGzjminGzjGzi_partieB <- function(lambda.j,k.j,lambda.i,k.i,lowerbnd=10^(-5),fac=0.5,tol=10^(-5)){
  
  ## This function computes a part of the term Cij
  ##
  ## Input :
  ## - \hat{p}_12.t and \hat{p}_13.t vectors from kernel estimation
  ## - two vectors of length t containing Weibull's estimated parameters \hat{lambda}_t,\hat{k}_t
  ##
  ## Used in : matcovNB_cji
  
  vala.j=fac*(lambda.j)^k.j
  vala.i=fac*(lambda.i)^k.i
  I <- integrate(f=funcEGzjminGzjGzi_partieB,
                 lower=lowerbnd,upper=1,
                 subdivisions=1000L,
                 rel.tol=10^(-5),
                 lam.j=lambda.j,
                 k.j=k.j,
                 lam.i=lambda.i,
                 k.i=k.i,
                 a.j=vala.j,
                 a.i=vala.i,
                 lowerbnd=lowerbnd,
                 stop.on.error = FALSE)
  return(I$value)
}

calculEGzjminGziGzj_partieA1 <- function(lambda.j,k.j,lowerbnd=10^(-6),fac=0.5){
  
  ## This function computes a part of the term C1_ji
  ##
  ## Input :
  ## - lambda_tj and k_tj
  ##
  ## Used in : matcovNB_cji
  
  vala=fac*(lambda.j)^k.j
  I <- integrate(f=foncpartiedeA1EGjminGjGi,lower=lowerbnd,upper=1,subdivisions=1000L,
                 lam=lambda.j,k=k.j,a=vala,
                 stop.on.error = FALSE)
  return(I$value) 
}

foncpartiedeA1EGjminGjGi <- function(x,lam,k,a){
  # Utilitary function used inside function calculEGzjminGziGzj_partieA1()
  (1/a) * exp( -(2*lam/a^(1/k)) * (-log(x))^(1/k) ) * x^( 2/a - 1)
}

##### varcovarNA_t #############################################################################

matcovNA_t <- function(p12_hat_t, p13_hat_t,lambda_t,k_t,J){
  
  ## For a given time step t ,this routine computes the unweighed variance-covariance matrix of NA
  ##
  ##    matcov_NA_unweighed (t) =  ( A2_jit & D2_jit \\ C2_jit & B2_jit )
  ##
  ## Input :
  ## - index_t : year t for which we calculate varcov matrix
  ## - \hat{p}_12.t and \hat{p}_13.t vectors from kernel estimation
  ##
  ## Output :
  ## - Array of length J containg 2x2 unweighed covariance matrices of NA
  ##
  ## Used in : matcovp12p13_t
  ##
  ## Requires : matcovNA_A2_jit(), matcovNA_B2_jit(), matcovNA_C2_jit()
  
  cat("start: computation of matcov NA","\n")
  
  p12_hat_t<- as.vector(p12_hat_t)
  p13_hat_t<- as.vector(p13_hat_t)
  
  p1rfar_r4 <- p1rfarW_temps(lambda_t,k_t,matrix(4,ncol=1,nrow=J)) 
  p14W_t <- p1rfar_r4$p1r
  
  p1rfar_r5 <- p1rfarW_temps(lambda_t,k_t,matrix(5,ncol=1,nrow=J)) 
  p15W_t <- p1rfar_r5$p1r
  
  A2_jit <- matcovNA_A2_jit(p12_hat_t, p13_hat_t)
  B2_jit <- matcovNA_B2_jit(p13_hat_t, p15W_t)
  C2_jit <- matcovNA_C2_jit(p12_hat_t, p13_hat_t,p14W_t)
  
  list_unweighed_NA_t<- array(NA,c(2,2,J))
  
  list_unweighed_NA_t[1,1,] <- as.numeric(A2_jit)
  list_unweighed_NA_t[1,2,] <- as.numeric(C2_jit)
  list_unweighed_NA_t[2,1,] <- as.numeric(C2_jit)
  list_unweighed_NA_t[2,2,] <- as.numeric(B2_jit)
  
  return (list_unweighed_NA_t) 
}

matcovNA_A2_jit <- function(p12_hat_t,p13_hat_t){

  ## This fucntion computes the unweighted list A_2ij of NA's variance-covariance matrix
  ##
  ##    A2_ji = p13_hat.ti-(p12_hat.ti)^2
  ##
  ## Input :
  ## - index_t : year t for which are calculating the A2_ji terms (each time step t has associated a different vector A2_ji)
  ## - \hat{p}_12.t vector from kernel estimation
  ##
  ## Output :
  ## - Vector containing A2_ji terms
  ##
  ## Used in : matcovNA_t
  ##
  ## Requires : nothing
  
  cat("start computation of A2ji","\n")
  list_unweighed_matcov_NA_aji <- p13_hat_t - (p12_hat_t)^2
  cat("end computation of A2ji","\n")
  return (list_unweighed_matcov_NA_aji)
}


matcovNA_B2_jit <- function(p13_hat_t,p15W_t){

  ## This fucntion computes the unweighted list B2_ij of NA's variance-covariance matrix
  ##
  ##    B2_ji = p135W.ti - (p13_hat.ti)^2
  ##
  ## Input :
  ## - index_t : year t for which are calculating the A2_ji terms (each time step t has associated a different vector B2_ji)
  ## - \hat{p}_13.t vector from kernel estimation
  ##
  ## Output :
  ## - Vector containing A_2ji terms
  ##
  ## Used in : matcovNA_t
  ##
  ## Requires : nothing
  
  cat("start: computation of B2ji","\n")
  list_unweighed_matcov_NA_bji <- p15W_t -  (p13_hat_t)^2
  cat("end computation of B2ji","\n")
  return (list_unweighed_matcov_NA_bji)
}

matcovNA_C2_jit <- function(p12_hat_t,p13_hat_t,p14W_t){

  ## This fucntion computes the unweighted list C2_ji of NA's variance-covariance matrix
  ##
  ##    C2_ji = p14W_ti - (p12_hat.ti * p13_hat.ti)
  ##
  ## Input :
  ## - index_t : year t for which are calculating the C2_ji terms (each time step t has associated a different vector B2_ji)
  ## - \hat{p}_12.t and \hat{p}_13.t vectors from kernel estimation
  ##
  ## Output :
  ## - Vector containing A2_ji terms
  ##
  ## Used in : matcovNA_t
  ##
  ## Requires : nothing
  
  cat("start: computation of C2ji","\n")
  list_unweighed_matcov_NA_cji <- p14W_t - (p12_hat_t * p13_hat_t)
  cat("end computation of C2ji","\n")
  return (list_unweighed_matcov_NA_cji)
}

###########################   Variance of p1rW_t and  far_t}   ################################
###########################                                    ################################

#### Var {lambda_t, k_t}  #####################################################################

matcovtheta_t <- function(matcovN_t, X.vec, Z.vec, GmZ,tt,t_eval,h){
  # output : array num[1:2, 1:2, 1:J]
  J <- length(Z.vec)
  theta_t <-weibullGMM_NonStationaire (GmZ, tt, t_eval, h, kern=dEpan, truevalues=NULL)
  lam_t <- theta_t[[1]]
  k_t <- theta_t[[2]]
  Jacov <- jacobianFunctiong12(lam_t,k_t)
  
  list_matcovp <- array(NA,c(2,2,J))
  
  for (i in 1:J){
    list_matcovp[,,i]<- inv(Jacov[[1]][[i]]) %*% matcovN_t[,,i] %*% inv(t(Jacov[[1]][[i]]))
  }
  # Warning messages: 
  #   1: In inv(inv(Jacov[[1]][[i]]) %*% matcovN.vec[, , i] %*% t(Jacov[[1]][[i]])) :
  #   Matrix appears to be singular.
  # 2: In inv(inv(Jacov[[1]][[i]]) %*% matcovN.vec[, , i] %*% t(Jacov[[1]][[i]])) :
  #   Matrix appears to be singular.
  return (list_matcovp)
}

#### Var{p1r_t} and Var{far_t}  ###############################################################
varp1rfar_t <- function(r,matcovN_t, X.vec, Z.vec, GmZ,tt,t_eval,h){
  J <- length(Z.vec)
  theta_t <-weibullGMM_NonStationaire (GmZ, tt, t_eval, h, kern=dEpan, truevalues=NULL)
  lam_t <- theta_t[[1]]
  k_t <- theta_t[[2]]
  p1rW_t <- p1rfarW_temps(as.matrix(lam_t),as.matrix(k_t),as.matrix(rep(r,J)))
  
  Jacov12 <- jacobianFunctiong12(lam_t,k_t)
  Jacovrminus1 <- jacobianFunctiongrminus1(lam_t,k_t,rep(r,J))
  
  list_variancep1r_t <- rep(0,J)
  list_variancefar_t <- rep(0,J)
  for (i in 1:J){
    Jacov12_inv <- solve(Jacov12[[1]][[i]])
    Jacov12T_inv <- solve(t(Jacov12[[1]])[[i]])
    list_variancep1r_t[i] <- Jacovrminus1[[i]]%*% Jacov12_inv %*% matcovN_t[,,i] %*% Jacov12T_inv %*% t(Jacovrminus1[[i]])
    #list_variancefar_t[i] <- list_variancep1r_t[i] / (sqrt(J)*r*((p1rW_t$p1r[i])^2))
  }
  return (list("varp1r_t"=list_variancep1r_t,"p1r_t"=p1rW_t$p1r))
}

## Auxiliary function
jacobianFunctiongrminus1 <- function(lam.vec,k.vec,r.vec){

  lvec <- length(lam.vec)
  listejacobiennes <- list()
  for (i in 1:lvec){
    lambda <- lam.vec[i] ; k <- k.vec[i] ; r <- r.vec[i]
    dgoverdlambda <- dgjoverdlambdafunc(r-1,lambda,k)
    dgoverdk <- dgjoverdkfunc(r-1,lambda,k)
    listejacobiennes[[i]] <- matrix(c(dgoverdlambda,dgoverdk),nrow=1,ncol=2)
  }
  return(listejacobiennes)
}

##################   Confidence intervals for p1rW_t and  far_t   ############################
##################                                                 ############################

# temporary
CI_p1rfar <- function(p1r.vec,varp1r.vec,r,t_eval,J,alpha=0.5){
  
  stdp1rW <- sqrt(varp1r.vec)
  zalpha <- qnorm(1-alpha/2)
  lowerbndp1r_t <- p1r.vec*exp(-(zalpha*stdp1rW)/sqrt(J*p1r.vec^2))
  upperbndp1r_t <- p1r.vec*exp(+(zalpha*stdp1rW)/sqrt(J*p1r.vec^2))
  
  plot(t_eval,upperbndp1r_t, type="l",col="gray")
  lines(t_eval,p1r.vec, col="darkgreen")
  lines(t_eval,lowerbndp1r_t,col="gray")
  
  return (list("lowp1r_t"=lowerbndp1r_t,"uppperp1r_t"=upperbndp1r_t))
}

# General fonction
calcul_ICp1r <- function(r,X.vec, Z.vec,tt,t_eval,h,alpha=0.05){

  ## This fucntion computes the upper and lower confidence bounds for p1r
  ##
  ## Input :
  ## - record size r (numeric value)
  ## - counterfactual and factual X.vec and Z.vec (both vectors)
  ## - tt vector containing times steps of Z (vector of the same length of Z.vec)
  ## - t_val , vector of evaluation times
  ## - h: bad+ndwith (numeric value)
  ## - alpha: confidence level (numeric value)
  ##
  ## Output :
  ## - 3 vectors: upper bound, lower bound, p1rt
  ## - 2 plots: the first show us the 3 vectors for a y scale allowing us to differenciate them
  ##            the sencond one has an y scale that allows us to compare p1r to 1/r
  ##
  ## Requires : weibullGMM_NonStationaire(), matcovNB(), matcovNA_t(), varp1rfar_t(), 
  ##
  ## REMARK: WHEN CHANGING r, WE MUST CHANGE LABELS OF AXIS AND TITLE TO THE CORRECT r
  
  J <- length(Z.vec)
  I <- length(X.vec)
  
  G_emp <- ecdf(X.vec)
  GmZ <- G_emp(Z.vec)
  
  cat("start: computation of theta","\n")
  theta_t <-weibullGMM_NonStationaire (GmZ, tt, t_eval, h, kern=dEpan, truevalues=NULL)
  lambda_t <- theta_t[[1]]
  k_t <- theta_t[[2]]
  cat("end: computation of theta","\n")
  
  p12_hat.vec <- p12_NonPar(X.vec,Z.vec,tt,t_eval,h)
  p13_hat.vec<- p13_NonPar(X.vec,Z.vec,tt,t_eval,h)
  
  list_weighed_matcov_N_t<- array(NA,c(2,2,J))
  
  cat("start: computation of unweighted matcovNB","\n")
  list_unweighed_matcov_NB <- (1/I) * matcovNB(p12_hat.vec, p13_hat.vec,lambda_t,k_t)
  cat("end: computation of matcovNB","\n")
  
  Kh <- outer(t_eval, tt, function(zz,z) dEpan((zz - z) / h))
  
  for (index_t in 1:J){
    cat("start: computation of unweighted matcovNA for t=" ,index_t,"\n")
    list_unweighed_matcov_NA_t<-matcovNA_t(index_t, p12_hat.vec, p13_hat.vec,lambda_t,k_t,J)

    list_unweighed_matcov_N_t<- list_unweighed_matcov_NB + list_unweighed_matcov_NA_t
    
    cat("end: unweighed matcov computed for t= ",index_t,"\n")
    
    Aji <- list_unweighed_matcov_N_t[1,1,]
    Bji <- list_unweighed_matcov_N_t[2,2,]
    Dji <- list_unweighed_matcov_N_t[1,2,]
    Cji <- list_unweighed_matcov_N_t[2,1,]
    
    denom<-(sum(Kh[index_t,]))^2
    
    list_Wji_t <- c()
    cat("start: weighting matcovN","\n")
    for (j in 1:dim(Kh)[2]){
      Khj<-Kh[index_t,j]
      
      for (i in 1:j){
        Khi <- Kh[index_t,i]
        
        KhjKhi <- Khj * Khi
        
        list_Wji_t <- c(list_Wji_t, KhjKhi)
      }
    }
    list_weighed_matcov_N_t[1,1,index_t] <- 1/(J^2) *sum(list_Wji_t * Aji) / denom
    list_weighed_matcov_N_t[1,2,index_t] <- 1/(J^2) *sum(list_Wji_t * Dji) /denom
    list_weighed_matcov_N_t[2,1,index_t] <- 1/(J^2) *sum(list_Wji_t * Cji) /denom
    list_weighed_matcov_N_t[2,2,index_t] <- 1/(J^2) *sum(list_Wji_t * Bji) /denom
    cat("end: weighting matcovN ","\n")
  }
  cat("start: CI computation","\n")
  matcov_p1rfar<-varp1rfar_t(r,list_weighed_matcov_N_t, X.vec, Z.vec, GmZ,tt,t_eval,h)
  p1r<-matcov_p1rfar$p1r_t
  stdp1rW <- sqrt(matcov_p1rfar$varp1r_t)
  zalpha <- qnorm(1-alpha/2)
  lowerbndp1r_t <- p1r - zalpha*stdp1rW
  upperbndp1r_t <-p1r + zalpha*stdp1rW
  cat("end: CI computation","\n")
  
  #PLOT 1
  cat("start: plot1","\n")
  plot(t_eval,upperbndp1r_t, type="l",col="gray",xlab="time",ylab=expression(p["1,10,t"]),
       main=expression(Evolution ~  over  ~ time  ~ of ~ p["1,10"]))
  lines(t_eval,lowerbndp1r_t,col="gray")
  lines(t_eval,p1r,col="blue")
  
  # PLOT 2
  cat("start: plot2","\n")
  plot(t_eval,lowerbndp1r_t, type="l",ylim=c((1/r)-0.1,1),xlab="time",ylab=expression(p["1,10,t"]),
       main=expression(Evolution ~  over  ~ time  ~ of ~ p["1,10"]))
  lines(t_eval,lowerbndp1r_t)
  lines(t_eval,p1r,col="blue")
  lines(t_eval,rep(1/r,J),col="red")
  
  return (list("lowp1r_t"=lowerbndp1r_t,"uppperp1r_t"=upperbndp1r_t,"p1r_t"=p1r))
}

matcovNA_alone <- function(X,Z,tt,t_eval,h,lambda_h,k_h,graphiques=TRUE){
  
  # weighed variance of NA
  
  Z<-as.numeric(Z)
  X <- as.numeric(X)
  J <- length(Z)
  p12_hat.vec <- as.vector(p12_NonPar(X,Z,tt,t_eval,h)$p12.mat)
  p13_hat.vec<- p13_NonPar(X,Z,tt,tt,h)
  G_emp <- ecdf(X);matGm<- G_emp(Z)
  tetha<-weibullGMM_NonStationaire(matGm, tt, t_eval, h, kern=dEpan, truevalues=NULL)
  lambda_h <- tetha$lambdahat
  k_h <- tetha$khat
  p14W_t <- p1rfarW_temps(lambda_h,k_h,matrix(4,ncol=1,nrow=J))$p1r
  p15W_t <- p1rfarW_temps(lambda_h,k_h,matrix(5,ncol=1,nrow=J))$p1r
  A2<-matcovNA_A2_jit(p12_hat.vec,p13_hat.vec)
  B2<-matcovNA_B2_jit(p13_hat.vec,p15W_t)
  C2<-matcovNA_C2_jit(p12_hat.vec,p13_hat.vec,p14W_t)
  Kh <- outer(t_eval, tt, function(zz,z) dEpan((zz - z) / h))
  list_unweighed_NA_t<- array(NA,c(2,2,J))
  list_unweighed_NA_t[1,1,] <- as.numeric(A2)
  list_unweighed_NA_t[1,2,] <- as.numeric(C2)
  list_unweighed_NA_t[2,1,] <- as.numeric(C2)
  list_unweighed_NA_t[2,2,] <- as.numeric(B2)
  list_weighed_matcov_N_t<- array(NA,c(2,2,J)) #array
  for (index_t in 1:J){
    denom<-sum(Kh[index_t,])
    list_khj_t<-c()
    for (j in 1:dim(Kh)[2]){
      Khj<-Kh[index_t,j]
      list_khj_t <- c(list_khj_t, (Khj/denom)^2)
    }
    list_weighed_matcov_N_t[1,1,index_t] <- 1/(J^2) * sum(list_khj_t*2*A2) 
    list_weighed_matcov_N_t[1,2,index_t] <-  1/(J^2) *sum(list_khj_t*2*C2)
    list_weighed_matcov_N_t[2,1,index_t] <- 1/(J^2) *sum(list_khj_t*2*C2) 
    list_weighed_matcov_N_t[2,2,index_t] <- 1/(J^2) *sum(list_khj_t*2*B2) 
  }
  if (graphiques==TRUE){
    plot(tt,list_weighed_matcov_N_t[1,1,])
    plot(tt,list_weighed_matcov_N_t[1,2,])
    plot(tt,list_weighed_matcov_N_t[2,2,])
  }
  return(list_weighed_matcov_N_t)
}

matcovNB_alone <- function(I,X,Z,tt,t_eval,h,lambda_h,k_h,graphiques=TRUE){
  
  ## weighed variance of NB
  
  X <- as.numeric(X)
  Z <- as.numeric(Z)
  J<-length(Z)
  I <- length(X)
  uw_matcovNB_a<-matcovNB_aji(p12_hat.vec,lambda_h,k_h)
  uw_matcovNB_b <- matcovNB_bji(p13_hat.vec,lambda_h,k_h)
  uw_matcovNB_c <- matcovNB_cji(p12_hat.vec,p13_hat.vec,lambda_h,k_h)
  Kh <- outer(t_eval, tt, function(zz,z) dEpan((zz - z) / h))
  list_weighed_matcov_NB_t<- array(NA,c(2,2,J)) #array
  for (index_t in 1:J){
    denom<-(sum(Kh[index_t,]))^2
    list_Wji_t <- c()
    for (j in 1:dim(Kh)[2]){
      Khj<-Kh[index_t,j]
      for (i in 1:j){
        Khi <- Kh[index_t,i]
        KhjKhi <- (Khj * Khi)/denom
        list_Wji_t <- c(list_Wji_t, KhjKhi)
      }
    }
    list_weighed_matcov_NB_t[1,1,index_t] <- (1/I) * 1/(J^2) * sum(list_Wji_t * 2*uw_matcovNB_a )
    list_weighed_matcov_NB_t[1,2,index_t] <- (1/I) * 1/(J^2) * sum(list_Wji_t * 2*uw_matcovNB_c )  
    list_weighed_matcov_NB_t[2,1,index_t] <- (1/I) * 1/(J^2) * sum(list_Wji_t * 2*uw_matcovNB_c)   
    list_weighed_matcov_NB_t[2,2,index_t] <- (1/I) * 1/(J^2) * sum(list_Wji_t * 2*uw_matcovNB_b) 
    cat("time", index_t, "weighed","\n")
  }
  
  if (graphiques==TRUE){
    plot(tt,list_weighed_matcov_NB_t[1,1,])
    plot(tt,list_weighed_matcov_NB_t[1,2,]) 
    plot(tt,list_weighed_matcov_NB_t[2,2,]) 
  }
  return(list_weighed_matcov_NB_t)
}

  
# KEEP IN LIKE COMENT UNTIL NB IS DEVELOPED FOR ALL GRIDPOINTS
# #############  Matcov Na for all points of the map at the same time ###############
# ### ( Work in progress, not yet implemented as it is not finished for matcovNB)  ##
# ###################################################################################
# 
# matcovNA <- function(Kh.mat,p12.mat,p13.mat,lam.mat,k.mat,J,N){
#   
#   # les imputs .mat doivent être des matrices
#   p12.mat <- as.matrix(p12.mat)
#   p13.mat <- as.matrix(p13.mat)
#   lam.mat <- as.matrix(lam.mat)
#   k.mat <- as.matrix(k.mat)
#   
#   p14.mat <- p1rfarW_temps(lam.mat,k.mat,r.mat=matrix(4,J,N))$p1r
#   p15.mat <- p1rfarW_temps(lam.mat,k.mat,r.mat=matrix(5,J,N))$p1r
#   
#   Kh2.mat <- Kh.mat^2
#   
#   covNA_A.mat <- matcovNA_A(Kh2.mat,p12.mat,p13.mat,J,N)
#   covNA_B.mat <- matcovNA_B(Kh2.mat,p13.mat,p15.mat,J,N)
#   covNA_C.mat <- matcovNA_C(Kh2.mat,p12.mat,p13.mat,p14.mat,J,N)
#   
#   matcovNA<-list()
#   for (gridpoint in 1:N){
#     array_matcovNA_gridpoint <- array(NA,c(2,2,J))
#     array_matcovNA_gridpoint[1,1,] <- covNA_A.mat[,gridpoint]
#     array_matcovNA_gridpoint[2,2,] <- covNA_B.mat[,gridpoint]
#     array_matcovNA_gridpoint[1,2,] <- covNA_C.mat[,gridpoint]
#     array_matcovNA_gridpoint[2,1,] <- covNA_C.mat[,gridpoint]
#     
#     matcovNA <- c(matcovNA, list(array_matcovNA_gridpoint) )
#   }
#   return(matcovNA) # list of -> (#1 num [1:2, 1:2, 1:30]), #2...) 
# }
# 
# matcovNA_A <- function(Kh2.mat,p12.mat,p13.mat,J,N){
#   covNA_A_map <- matrix(NA,nrow = J, ncol = N)
#   for (gridpoint in 1:N){ # on doit créer un vecteur de variances pour chaque gridpoint
#     p1rmix.vec <- p13.mat[,gridpoint] - p12.mat[,gridpoint]^2
#     Khp1r.mat <- sweep(Kh2.mat, MARGIN=1, p1rmix.vec, `*`)# by row
#     Khp1r.vec <- rowSums(Khp1r.mat) 
#     # Khp1r.vec est un vecteur car à chaque pas de temps la variance sera different
#     covNA_A_map[,gridpoint]<-Khp1r.vec 
#   }
#   return(covNA_A_map) # chaque colonne est associé à une position differente
# }
# 
# matcovNA_B <- function(Kh2.mat,p13.mat,p15.mat,J,N){
#   covNA_B_map <- matrix(NA,nrow = J, ncol = N)
#   for (gridpoint in 1:N){ # on doit créer un vecteur de variances pour chaque gridpoint
#     p1rmix.vec <- p15.mat[,gridpoint] - p13.mat[,gridpoint]^2
#     Khp1r.mat <- sweep(Kh2.mat, MARGIN=1, p1rmix.vec, `*`) # by row
#     Khp1r.vec <- rowSums(Khp1r.mat) 
#     # Khp1r.vec est un vecteur car à chaque pas de temps la variance sera different
#     covNA_B_map[,gridpoint]<-Khp1r.vec 
#   }
#   return(covNA_B_map) # chaque colonne est associé à une position differente
# }
# 
# matcovNA_C <- function(Kh2.mat,p12.mat,p13.mat,p14.mat,J,N){
#   covNA_C_map <- matrix(NA,nrow = J, ncol = N)
#   for (gridpoint in 1:N){ # on doit créer un vecteur de variances pour chaque gridpoint
#     p1rmix.vec <- p14.mat[,gridpoint] - p12.mat[,gridpoint] * p13.mat[,gridpoint]
#     Khp1r.mat <- sweep(Kh2.mat, MARGIN=1, p1rmix.vec, `*`)
#     Khp1r.vec <- rowSums(Khp1r.mat) 
#     # Khp1r.vec est un vecteur car à chaque pas de temps la variance sera different
#     covNA_C_map[,gridpoint]<-Khp1r.vec 
#   }
#   return(covNA_C_map) # chaque colonne est associé à une position differente
# }
# 
# #########################  Work in progress for matcov NB #########################
# ###################################################################################
# ###################################################################################
# 
# matcovp12p13_t <- function(X.mat, Z.mat,tt,t_eval,r.mat,h){
#   
#   X.mat <- as.matrix(X.mat)
#   Z.mat <- as.matrix(Z.mat)
#   r.mat <- as.matrix(r.mat)
#   
#   J <- dim(Z.vec)[1]
#   N <- dim(Z.vec)[2]
#   I <- dim(X.vec)[1]
#   
#   p12p13 <- P12_P13_estimation(X.mat,Z.mat,tt,t_eval,h,kern=dEpan)
#   GmZ.mat <- p12p13$matGmZ
#   p12_hat.vec <- p12p13$matp12
#   p13_hat.vec <- p12p13$matp13
#   
#   
#   theta.mat <-weibullGMM_NonStationaire (GmZ.mat, tt, t_eval, h, kern=dEpan, truevalues=NULL)
#   lambda.mat <- theta.mat[[1]]
#   k.mat <- theta.mat[[2]]
#   
#   Kh <- outer(t_eval, tt, function(zz,z) dEpan((zz - z) / h))
#   
#   cov_NA <- matcovNA(Kh, p12_hat.vec, p13_hat.vec, lambda.mat, k.mat, J, N) # output, list of arrays. at each t, var is different
#   
#   cov_NB <- (1/I) * matcovNB(Kh,p12_hat.vec, p13_hat.vec,lambda_t,k_t) # add Kh . out out list of matrix. a chaque t, var est le même
#   # mais olus simple si on crée des arrays aussi 
#   
#   cov_N<-list()
#   for (gridpoint in 1:N){
#     array_matcovN_gridpoint <- cov_NA[[gridpoint]]+cov_NB[[gridpoint]]
#     cov_N <- c(cov_N, list(array_matcovN_gridpoint) )# list of arrays
#   }
#   return(cov_N)
# }
# 
# #######
# # matcovNB <- function(Kh.mat,p12.mat,p13.mat,lam.mat,k.mat,J,N)
# matcovNB <- function(Kh.mat,p12.mat, p13.mat,lam.mat,k.mat,J,N){
#   
#   # Calcul "poids"
#   Wji_t.vec <- list()
#   for (index_t in 1:dim(Kh)[1]){
#     denom<-(sum(Kh[index_t,]))^2
#     Wji_indext <- c()
#     for (j in 1:dim(Kh)[2]){
#       Khj<-Kh[index_t,j]
#       for (i in 1:j){
#         Khi <- Kh[index_t,i]
#         KhjKhi <- Khj * Khi
#         Wji_indext <- c(Wji_indext, KhjKhi)
#       }
#     }
#     Wji_t.vec [[index_t]] <- Wji_indext
#   }
#   WjtWit.mat = do.call(cbind, WjtWit)
#   
#   covNB_A.mat <- matcovNB_aji(WjtWit.mat, p12.mat, lam.mat, k.mat,J,N)
#   covNB_B.mat <- matcovNB_bji(WjtWit.mat, p13.mat, lam.mat, k.mat,J,N)
#   covNB_C.mat <- matcovNB_cji(WjtWit.mat, p12.mat, p13.mat, lam.mat, k.mat,J,N) #matrix where col:dif grid, row=t_evalx
#   
#   matcovNB<-list()
#   for (gridpoint in 1:N){
#     array_matcovNB_gridpoint <- array(NA,c(2,2,J))
#     array_matcovNB_gridpoint[1,1,] <- covNB_A.mat[,gridpoint]
#     array_matcovNB_gridpoint[2,2,] <- covNB_B.mat[,gridpoint]
#     array_matcovNB_gridpoint[1,2,] <- covNB_C.mat[,gridpoint]
#     array_matcovNB_gridpoint[2,1,] <- covNB_C.mat[,gridpoint]
#     
#     matcovNB <- c(matcovNB, list(array_matcovNB_gridpoint) )
#   }
#   return(matcovNB) # list of -> (#1 num [1:2, 1:2, 1:30]), #2...) 
# }
# 
# matcovNB_aji <- function(WjtWit.mat,p12.mat,lam.mat,k.mat,J,N){
#   cat("start: computation of A1ji","\n")
#   
#   #################
#   ################
#   covNA_A_map <- matrix(NA,nrow = J, ncol = N)
#   for (gridpoint in 1:N){ # on doit créer un vecteur de variances pour chaque gridpoint
#     p1rmix.vec <- p13.mat[,gridpoint] - p12.mat[,gridpoint]^2
#     Khp1r.mat <- sweep(Kh2.mat, MARGIN=1, p1rmix.vec, `*`)# by row
#     Khp1r.vec <- rowSums(Khp1r.mat) 
#     # Khp1r.vec est un vecteur car à chaque pas de temps la variance sera different
#     covNA_A_map[,gridpoint]<-Khp1r.vec 
#   }
#   return(covNA_A_map) # chaque colonne est associé à une position differente
#   ################
#   ##############
#   covNB_A_map <- matrix(NA,nrow = J, ncol = N)
#   for (gridpoint in 1:N){ 
#     for (j in 1:J){
#       lambda.j <- lambda_t[j]
#       k.j <- k_t[j]
#       p12_hat.j <- p12_hat_t[j]
#     }
#     ####
#     for (j in 1:length(p12_hat_t)){
#       cat("start: components for j= ",j,"\n")
#       lambda.j <- lambda_t[j]
#       k.j <- k_t[j]
#       p12_hat.j <- p12_hat_t[j]
#       
#       for (i in 1:j){
#         cat("start: components for i= ",i,"\n")
#         lambda.i <- lambda_t[i]
#         k.i <- k_t[i]
#         p12_hat.i <- p12_hat_t[i]
#         
#         matcov_NB_aji <-Mrfuncji(2,lambda.j,k.j,lambda.i,k.i) - p12_hat.j * p12_hat.i
#         
#         list_unweighed_matcov_NB_aji <- c(list_unweighed_matcov_NB_aji, matcov_NB_aji)
#         cat("done: unweighed list of components for ",j,", ",i,"\n")
#       }
#     }
#     return (list_unweighed_matcov_NB_aji)
#   }
#   
#   matcovNB_cji <- function(p12_hat_t,p13_hat_t,lambda_t,k_t,J,N){
#     
#     cat("start: computation of C1ji","\n")
#     
#     list_unweighed_matcov_NB_cji <- c()
#     
#     for (j in 1:length(p12_hat_t)){
#       cat("start: components for j= ",j,"\n")
#       lambda.j <- lambda_t[j]
#       k.j <- k_t[j]
#       p12_hat.j <- p12_hat_t[j]
#       p13_hat.j <- p13_hat_t[j]
#       
#       EGzjminGzjGzi_partieA1 <-calculEGzjminGziGzj_partieA1(lambda.j,k.j,lowerbnd=10^(-6),fac=0.5)
#       
#       for (i in 1:j){
#         cat("start: components for i= ",i,"\n")
#         lambda.i <- lambda_t[i]
#         k.i <- k_t[i]
#         p12_hat.i <- p12_hat_t[i]
#         p13_hat.i <- p13_hat_t[i]
#         
#         EGzjminGzjGzi_partieB <- calculEGzjminGzjGzi_partieB(lambda.j,k.j,lambda.i,k.i,lowerbnd=10^(-5),fac=0.5,tol=10^(-5))
#         
#         # Adding all parts
#         EGzjminGzjGzi <- EGzjminGzjGzi_partieA1 + p13_hat.j + EGzjminGzjGzi_partieB 
#         
#         matcov_NB_cji <- 2*(EGzjminGzjGzi + (p13_hat.j * p12_hat.i))
#         
#         list_unweighed_matcov_NB_cji <- c(list_unweighed_matcov_NB_cji, matcov_NB_cji)
#         cat("done: unweighed list of components for ",j,", ",i,"\n")
#       }
#     }
#     return (list_unweighed_matcov_NB_cji)
#   }
#   
#   matcovNB_bji <- function(p13_hat_t,lambda_t,k_t,J,N){
#     
#     
#     cat("start: computation of B1ji","\n")
#     
#     list_unweighed_matcov_NB_bji <- c()
#     
#     for (j in 1:length(p13_hat_t)){
#       cat("start: components for j= ",j,"\n")
#       lambda.j <- lambda_t[j]
#       k.j <- k_t[j]
#       p13_hat.j <- p13_hat_t[j]
#       
#       for (i in 1:j){
#         cat("start: components for i= ",i,"\n")
#         lambda.i <- lambda_t[i]
#         k.i <- k_t[i]
#         p13_hat.i <- p13_hat_t[i]
#         
#         NB_bji <- 4*(Mrfuncji(3,lambda.j,k.j,lambda.i,k.i) + p13_hat.j * p13_hat.i)
#         
#         list_unweighed_matcov_NB_bji <- c(list_unweighed_matcov_NB_bji, NB_bji)
#         cat("done: unweighed list for j= ", j," i= ",i,"\n")
#       }
#     }
#     return (list_unweighed_matcov_NB_bji)
#   }

##############   Functions that extract lam and k from CMIP simulations        ################
#############                                                                  ################

traj_from_data <- function(variable.df, grid_points, model.choice, run.choice,var="tmax"){ # tmax or pr
  
  # just extract trayectories: matx, matz and tt 
  
  if (var=="tmax"){
    hist_fin <- 2014
    rcp85_0 <- 2015
  }
  
  if (var=="pr"){
    hist_fin <- 2005
    rcp85_0 <- 2006
  }
  
  Z_historical <- variable.df %>% 
    select (institute,model, experiment, run, year, one_of(str_c(grid_points))) %>%
    filter(experiment == "historical" & model == model.choice & run == run.choice & between(year, 1850, hist_fin))  %>% 
    arrange(year) %>%
    select(!c(institute,model,experiment,run, year))
  Z_historical<-as.data.frame(Z_historical)
  Z_rcp85<- variable.df %>% 
    select (institute, model, experiment, run, year, one_of(str_c(grid_points))) %>%
    filter(experiment == "rcp85" & model == model.choice & run == run.choice & between(year, rcp85_0, 2100))  %>%
    arrange(year) %>% 
    select(!c(institute,model,experiment,run, year))
  Z_rcp85<-as.data.frame(Z_rcp85)
  matz <- bind_rows(Z_historical,Z_rcp85)
  matz <-as.matrix(matz)
  
  matx<- variable.df %>%
    select (institute, model, experiment, run, year, one_of(str_c(grid_points))) %>%
    filter(experiment == "historicalNat" & model == model.choice & run == run.choice)  %>%
    arrange(year) %>% select(!c(institute,model,experiment,run, year))
  matx<-as.matrix(matx)
  
  tt<-c(1:dim(matz)[1])
  
  return(list("matZ"=matz,"matX"=matx,"time.vec"=tt))
}

theta_map_from_data <-function(variable.df, model.choice, run.choice, varname, savemat=TRUE){
  
  ## for a given variable it creates the matrix matlam and matk,
  ## which contains the the values of the variables for each gridpoint
  ## and for each t from 1850 to 2100
  ## 
  ## Input
  ## variable.df: the data base containing Tmax trayectories for differents models, scenarios and runs, 
  ##          e.g. tmax_cmip6_yearmax.rds"
  ## model.choice: the name of the model in between "", e.g. IPSL-CM6A-LR"
  ## run.choice:  the name of the tun in between "", e.g. "r1i1p1f1"
  ## varname: name of the varible of study in between "", e.g. "tmax" or "pr"
  ##
  ## Output
  ## two files.rds, one with de values of lambda and the other with the values of k. the file's names are
  ## "matlam_",model.choice,"_",run.choice,".rds"
  ##
  ## Uses: matGZ_func, weibullGMM_NonStationaire
  ##
  
  stringYear<-toString(model.choice)
  stringYear<-toString(run.choice)
  
  Z_historical <- variable.df %>% 
    filter(experiment == "historical" & model == model.choice & run == run.choice & between(year, 1850, 2014))  %>% 
    arrange(year) %>%
    select(!c(institute,model,experiment,run, year))
  Z_historical<-as.data.frame(Z_historical)
  Z_rcp85<- variable.df %>% 
    filter(experiment == "rcp85" & model == model.choice & run == run.choice & between(year, 2015, 2100))  %>%
    arrange(year) %>% 
    select(!c(institute,model,experiment,run, year))
  Z_rcp85<-as.data.frame(Z_rcp85)
  matz <- bind_rows(Z_historical,Z_rcp85)
  matz <-as.matrix(matz)
  
  matx<- variable.df %>%
    filter(experiment == "historicalNat" & model == model.choice & run == run.choice)  %>%
    arrange(year) %>% select(!c(institute,model,experiment,run, year))
  matx<-as.matrix(matx)
  
  tt<-c(1:dim(matz)[1])
  
  matGmZ<-matGZ_func(matx,matz)
  mattetha <- weibullGMM_NonStationaire(matGmZ, tt, t_eval=tt, h=20, kern=dEpan, truevalues=NULL)
  matlam <- mattetha$lambdahat
  matk <- mattetha$khat
  
  if (savemat==TRUE){
    matlam <- as.data.frame(matlam)
    matk <- as.data.frame(matk)
    matlam_file_name <- paste("matlam_",varname,model.choice,"_",run.choice,".rds",sep ="")
    matk_file_name <- paste("matk_",varname,model.choice,"_",run.choice,".rds",sep ="")
    saveRDS(matlam, file=matlam_file_name)
    saveRDS(matk, file=matk_file_name)
  }

  return(list("matrix_lambda"=matlam,"matrix_k"=matk))

}

theta_point_traj_from_data <- function(variable.df, grid_points, model.choice, run.choice, varname,hopt=40,savemat=TRUE){
  
  ## for a given variable it creates the matrix matlam and matk,
  ## which contains the the values of the variables at the chosen grid points
  ## and for each t from 1850 to 2100
  ## 
  ## Input
  ## variable.df: the data base containing Tmax trayectories for differents models, scenarios and runs, 
  ##          e.g. tmax_cmip6_yearmax.rds"
  ## grid_points: vector with the numers of the column, e.g. c(8,500,245)
  ## model.choice: the name of the model in between "", e.g. IPSL-CM6A-LR"
  ## run.choice:  the name of the tun in between "", e.g. "r1i1p1f1"
  ## varname: name of the varible of study in between "", e.g. "tmax" or "pr"
  ##
  ## Output
  ##  - two files.rds, one with de values of lambda and the other with the values of k. the file's names are
  ## "matlam_",model.choice,"_",run.choice,".rds"
  ## - 4 matrix : "matrix_lambda","matrix_k","matZ","matX"
  ##
  ## Uses: matGZ_func, weibullGMM_NonStationaire
  ##
  
  Z_historical <- variable.df %>% 
    select (institute,model, experiment, run, year, one_of(str_c(grid_points))) %>%
    filter(experiment == "historical" & model == model.choice & run == run.choice & between(year, 1850, 2014))  %>% 
    arrange(year) %>%
    select(!c(institute,model,experiment,run, year))
  Z_historical<-as.data.frame(Z_historical)
  Z_rcp85<- variable.df %>% 
    select (institute, model, experiment, run, year, one_of(str_c(grid_points))) %>%
    filter(experiment == "rcp85" & model == model.choice & run == run.choice & between(year, 2015, 2100))  %>%
    arrange(year) %>% 
    select(!c(institute,model,experiment,run, year))
  Z_rcp85<-as.data.frame(Z_rcp85)
  matz <- bind_rows(Z_historical,Z_rcp85)
  matz <-as.matrix(matz)
  
  matx<- variable.df %>%
    select (institute, model, experiment, run, year, one_of(str_c(grid_points))) %>%
    filter(experiment == "historicalNat" & model == model.choice & run == run.choice)  %>%
    arrange(year) %>% select(!c(institute,model,experiment,run, year))
  matx<-as.matrix(matx)
  
  tt<-c(1:dim(matz)[1])
  
  matGmZ<-matGZ_func(matx,matz)
  mattetha <- weibullGMM_NonStationaire(matGmZ, tt, t_eval=tt, h=hopt, kern=dEpan, truevalues=NULL)
  matlam <- mattetha$lambdahat
  matk <- mattetha$khat
  
  if (savemat==TRUE){
    grid_points <- as.character(grid_points)
    grid_points_text<- paste("_",grid_points, sep="")
    grid_points_text <- paste(grid_points_text, collapse = "")
    grid_points_text <- toString(grid_points_text)
    matlam <- as.data.frame(matlam)
    matk <- as.data.frame(matk)
    matlam_file_name <- paste("matlam_",varname,model.choice,"_",run.choice,grid_points_text,".rds",sep ="")
    matk_file_name <- paste("matk_",varname,model.choice,"_",run.choice,grid_points_text,".rds",sep ="")
    saveRDS(matlam, file=matlam_file_name)
    saveRDS(matk, file=matk_file_name)
  }
  return(list("matrix_lambda"=as.matrix(matlam),"matrix_k"=as.matrix(matk),"matZ"=as.matrix(matz),"matX"=as.matrix(matx)))
}

##################     Functions allowing us to create maps         ############################
##################                                                  ############################

p1rfarW_yearT <- function(lam.mat,k.mat,r,indext,lowerbnd=10^(-5)){
  
  ## This function give containing de value of p1rt_hat at each grid point, for a given r and t
  ##
  ## Input :
  ## - lam.mat : matrix containing lambda parameters of all points at all times 
  ## - k.mat : matrix containing k parameters of all points at all times 
  ## - indext : numeric value, position in the vector of the year of interest
  ## - r :  numeric value, size record
  ##
  ## Outplut : a vector of length dim(lam.mat[2]) containing de values of p1rt at each grid point
  ## 
  ## Uses : laplaceWeibull
  
  lam.vec <- lam.mat[indext,]
  k.vec <- k.mat[indext,]
  
  Nbpoint = length(lam.vec)
  p1r.vec = rep(0, Nbpoint)
  
  for (i in 1:Nbpoint){# for each point
    p1r.vec[i] <- laplaceWeibull(r-1,lam.vec[i],k.vec[i],lowerbnd=lowerbnd)
  }  
  return(p1r.vec) # retourne un vector, mais on aimerait bien une matrice
}

library('gplots')
library('fields')
library('RColorBrewer')

Map_p1r <- function(p1r.vec,r,year){
  
  ## This function creates world maps of 36x72 where the color indicates de value of p1rt,
  ## when p1rt < 1/r the gridpoint is white
  ## Input :
  ## - p1r.vec : vector p1rt output of p1rfarW_yearT(), 
  ##             it must contain p1rt vaules for a given t at each gridpoint
  ## - r : record size (for deciding where the map is white)
  ## - year: number. for making the tittle of the plot
  ##
  ## Outplut : color map of dimension 36x72 

  p1r.vec[p1r.vec<1/r] <- NA
  gridmat_p1r <- matrix(data=p1r.vec,nrow=72,ncol=36)
  grid_lon<-seq(-177.5, 177.5, length.out = 72) 
  grid_lat<-seq(-87.5, 87.5, length.out = 36)
  rfcol <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  colsp <- rfcol(64)
  image.plot(grid_lon,grid_lat,gridmat_p1r,col = colsp,
             legend.args = list( text = expression(p["1rt"]),
                                 cex = 1,
                                 side = 3,
                                 line = .5),ylab="latitude",xlab="longitude")
  stringYear<-toString(year)
  title(stringYear,cex=0.5,line = .8)
  world(lwd = 1,add=TRUE) 
}

Map_far <- function(p1r.vec,r,year){
  
  ## This function creates world maps of 36x72 where the color indicates de value of far,
  ## when far < 0 the gridpoint is white
  ## Input :
  ## - p1r.vec : vector p1rt output of p1rfarW_yearT(), 
  ##             it must contain p1rt vaules for a given t at each gridpoint
  ## - r : for the computation of far from p1rt
  ## - year: number. for making the tittle of the plot
  ##
  ## Outplut : color map of dimension 36x72 
  
  far.vec<-1-(1/(r*p1r.vec))
  gridmat_far<-matrix(data=far.vec,nrow=72,ncol=36)
  gridmat_far[gridmat_far<0] <- NA
  grid_lon<-seq(-177.5, 177.5, length.out = 72) 
  grid_lat<-seq(-87.5, 87.5, length.out = 36)
  rfcol <- colorRampPalette(rev(brewer.pal(11,'Spectral')))
  colsp <- rfcol(64)
  image.plot(grid_lon,grid_lat,gridmat_far,col = colsp,
             legend.args = list( text = expression(far),
                                 cex = 1,
                                 side = 3,
                                 line = .5),ylab="latitude",xlab="longitude")
  stringYear<-toString(year)
  title(stringYear,cex=0.5,line = .8)
  world(lwd = 1,add=TRUE) 
}



