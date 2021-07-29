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

# Epanechnikov density distribution
dEpan <- function(x){
  k <- (3/4)* (1-x^2)
  k[-1>x] <- 0
  k[x>1] <- 0
  return (k)
} 

# Normal distribution
Kern_gauss <- function(x)exp(-x^2/2)

#### p12_t Non parametric computation ##############################################################
p12_NonPar <- function(X,Z,tt,t_eval,h,kern= dEpan){
  # G_Z computation
  G_emp <- ecdf(X)
  G_Z <- G_emp(Z)
  # W computation
  Kij <- outer(t_eval,tt,function(zz,z) kern((zz - z) / h))
  W <- Kij / rowSums(Kij)
  # p12_hat computation
  p12_hat <- W %*% G_Z
}

#### p13_t and p1r_t Non parametric computation ##############################################################

p13_NonPar <- function(X,Z,tt,t_eval,h,kern= dEpan){
  # G_Z computation
  G_emp_carre <- ecdf(X)
  G_Z_carre <- (G_emp_carre(Z))^2
  # W computation
  Kij <- outer(t_eval,tt,function(zz,z) kern((zz - z) / h))
  W <- Kij / rowSums(Kij)
  # p12_hat computation
  p13_hat <- W %*% G_Z_carre
}

p1r_NonPar <- function(X,Z,r,tt,t_eval,h,kern= dEpan){
  # G_Z computation
  G_emp_carre <- ecdf(X)
  G_Z_r <- (G_emp_carre(Z))^(r-1)
  # W computation
  Kij <- outer(t_eval,tt,function(zz,z) kern((zz - z) / h))
  W <- Kij / rowSums(Kij)
  # p12_hat computation
  p13_hat <- W %*% G_Z_r
}

################################### Confidence interval ###########################################

#### Auxiliary fonctions 
dEpan_2<- function(x){
  k <- ((3/4)* (1-x^2))^2
  k[-1>x] <- 0
  k[x>1] <- 0
  return (k)
}

kern_dens <- function(tt,t_eval,h,kern=dEpan) {
  fh <- numeric(length(t_eval))
  for (i in seq_along(tt))
    for (j in seq_along(t_eval))
      fh[i] <- fh[i] + kern((t_eval[i] - tt[j])/h) 
    fh <- fh / ( h * length(tt))
    return(fh)
}

#### p12_t's Confidence interval computation 
IC <- function(X,Z,tt,t_eval,h,kern=dEpan) {
  n <- length(Z)
  # computation of G_Z
  G_emp <- ecdf(X)
  G_Z <- G_emp(Z)
  #Computation of p12_hat
  p12_hat<-p12_NonPar(X,Z,tt,t_eval,h,kern)
  # integral of squared kernel
  K22<-integrate(dEpan_2,-1,1)$value
  # density distribution of t 
  f_t<-kern_dens(tt,tt,h)
  # variance of t 
  sigma_GZ <- as.numeric(var(G_Z-p12_hat))
  # Var(A):
  VarA <- (sigma_GZ * K22 * f_t )/(n*h)
  # Var(B):
  Khj <- outer(t_eval, tt, function(zz,z) dEpan((zz - z) / h))
  Khi <- outer(t_eval, tt, function(zz,z) dEpan((zz - z) / h))
  E <- mean(outer(G_Z,G_Z, pmin) - outer(G_Z,G_Z,"*"))
  VarB_num <- as.numeric(rowSums(Khj) %*% rowSums(Khi))
  VarB_denom <- length(tt)*(rowSums(Khj))^2
  VarB <- VarB_num*E/VarB_denom
  # sigma = Var(A)+ Var(B)
  sigma_m <- VarA + VarB
  # Computation of CI
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
      GZ_predict <-p12_NonPar(X,Z_tr,tt_tr,t_eval,h_using) 
      CV_err[i]<-(GZ_test - GZ_predict[i])^2
    }
    CV_err_h[j]<-(mean(CV_err))^0.5
  }
  return(CV_err_h) #vector of errors
}

#### General fonction #############################################################################

# Estimate p12_t and p13_t for each trajectory {(X)_t,(Z)_t}
P12_P13_estimation <- function(matX,matZ,tt,t_eval,h,kern=dEpan){
  
  matX <- as.matrix(matX)
  matZ <- as.matrix(matZ) 
  dimnsZ=dim(matZ)
  matp12 <- matrix(rep(0,prod(dimnsZ)),nrow=dimnsZ[1],ncol=dimnsZ[2])
  matp13 <- matrix(rep(0,prod(dimnsZ)),nrow=dimnsZ[1],ncol=dimnsZ[2])
  
  if (dimnsZ[1]!=length(tt)){
    stop("Z trajectory and tt must have the same length")
  }
  
  for (j in 1:dimnsZ[2]){
    X <- matX[,j]
    Z <- matZ[,j]
    matp12[,j] <- p12_NonPar(X,Z,tt,t_eval,h,kern)
    matp13[,j] <- p13_NonPar(X,Z,tt,t_eval,h,kern)
  }
  return(list("matp12"=matp12,"matp13"=matp13))
}

################################################################################################
############################ lambda_t and k_t computation ######################################
###############################  (multiple méthods) ############################################
################################################################################################

library(stats4); library(gmm); library(stats); library(np); library(EWGoF)

#### functions for E(G(Z)^j)

funcLaplace <- function(x,m,lam,k,a){ 
  (1/a) * exp( -(m*lam/a^(1/k))*(-log(x))^(1/k) ) * x^(1/a - 1) 
}

laplaceWeibull <- function(j,lambda,k,lowerbnd=10^(-6),upperbnd=1,fac=1,tol=10^(-5)){
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
  (-j/a^((1/k)+1)) * (-log(u))^(1/k) * exp( -(j*lam/a^(1/k))*(-log(u))^(1/k) ) * u^(1/a - 1) 
}

dgjoverdlambdafunc <- function(j,lambda,k,lowerbnd=10^(-6),fac=0.5){
  vala=fac*(j*lambda)^k
  I <- integrate(f=foncdgjoverdlambda,lower=lowerbnd,upper=1,subdivisions=1000L,
                 j=j,lam=lambda,k=k,a=vala,
                 stop.on.error = FALSE)
  return(I$value) 
}

foncdgjoverdk <- function(u,j,lam,k,a){
  (-lam/k^2) * log( (1/a)*(-log(u)) ) * foncdgjoverdlambda(u,j,lam,k,a)
}

dgjoverdkfunc <- function(j,lambda,k,lowerbnd=10^(-6),fac=1){
  vala=fac*(j*lambda)^k
  I <- integrate(f=foncdgjoverdk,lower=lowerbnd,upper=1,subdivisions=1000L,
                 j=j,lam=lambda,k=k,a=vala,stop.on.error = FALSE)
  return(I$value) 
}

################################################################################################
############################ fsolve method #####################################################

functionp12p13 <- function(theta,vecx){
  lambdaval=theta[1] ; kval = theta[2]
  p12 = vecx[1] ; p13 = vecx[2]
  p12val <- laplaceWeibull(j=1,lambda=lambdaval,k=kval,lowerbnd=10^(-8))
  p13val <- laplaceWeibull(j=2,lambda=lambdaval,k=kval,lowerbnd=10^(-8))
  M <- c( p12 - p12val ,  p13 - p13val )
  return(M)
}

fg <- function(theta){functionp12p13(theta,vecx=phat)}

J12 <- function(x){jacobianFunctiong12(x[1],x[2])}

weibullFsolve <- function(matp12,matp13,truevalues=NULL){
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

## fsolve_modif() is a modification of fsolev(), where we impose:
# [lambda,k]_xnew <- ([lambda,k]_x0)/2 when lamba<0
## We do not use the jacobian matrix

fsolve_modif <- function (f, x0, J = NULL, maxiter = 100, tol = .Machine$double.eps^(0.5), 
          ...) 
{
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

fgoptim<- function (theta,vecx) crossprod(functionp12p13(theta,vecx))

# fg2 <- function (theta){fgoptim(theta,vecx=phat)}
# Gradfg <- t(JX)%*%(functionp12p13(theta,phat)[[1]]) + crossprod(fg,JX) 

weibullOptim<- function(matp12,matp13,truevalues=NULL){
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

FastTestforp1r_gumbel <- function(tt,h,r,m,n,N,sigX.vec,sigZ.vec,muX.vec,muZ.vec,methode="optim"){
  theta_theo <- 1 / exp(muZ.vec)
  p1r_t_theo <- 1 / (1 + (r-1)*theta_theo)
  matX <- matrix(nrow=m,ncol=N)
  matZ <- matrix(nrow=n,ncol=N)
  matp <- matrix(nrow=n,ncol=N)
  for (j in 1:N){
    matX[,j] <- rgev(m, loc = muX.vec, scale = sigX.vec, shape = 0)
    matZ[,j] <- rgev(n, loc = muZ.vec, scale = sigZ.vec, shape = 0)
  }
  if (methode=="optim"){
    matp<-P12_P13_estimation(matX,matZ,tt,tt,h,kern=dEpan)
    thetahat<-weibullOptim(matp$matp12,matp$matp13,truevalues=NULL)
  }
  if (methode=="gmm"){
    GZ<-matGZ_func(matX,matZ)
    thetahat<-weibullGMM_NonStationaire(GZ, tt, tt, h, kern=dEpan, truevalues=NULL)
  }
  p1rfar<-p1rfarW_temps(thetahat[[1]],thetahat[[2]],matrix(r,ncol=N,nrow=n))
  p1r_mean <- rowMeans(p1rfar$p1r)
  plot(tt,p1r_t_theo,type="l",col="red",xlab="time",ylab="p1r_t", main=paste("evolution over time of p1r with r=",r))
  lines(tt,p1r_mean)
  for (i in 1:N){
    lines(tt,p1rfar$p1r[,i],col="gray")
  }
  return(list("matp1r"=p1rfar$p1r,"p1r_mean"=p1r_mean, "p1rmean"=p1r_mean))
}

################################################################################################
############################ Trajectory simulations ############################################
################################################################################################

#### Simulate W-class trajectories with non-sationnary Z 

simulWclass_nonStationnary_general <-function(m,N,ksiX,sigX=1,muX=0,ksiZ=NULL,muZ=NULL,sigZ=NULL,unknown="location", graph=TRUE){
  
  ## For Xt, ksi parameter must be given, sigma=1 and mu=0 by default, but can be modified
  ## For Zt, 2 from 3 parameters must be given, the third one must be unknown. We decide which one we want to compute in the option "unknown"
  ## Xt length m, Zt length will be taken from the length of its shape parameter ksi. If ksi is unknown, we will use the length of sigma.
  ## N is the number of Xt and Zt trajectories that we create.
  
  if (ksiX!=0){
    if (unknown == "location"){
      # We take our Z's sample size from the length of shape input
      n <- length(ksiZ)
      # We remind that if our shape input is just a value, the sample size will be 1
      if (n == 1){
        warning("Z's shape parameters length give us the size of this strajectory")
      }
      # Verification. Shape parameters must be of the same sign
      for (i in 1:n){
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
        sigZ <- rep(sigZ,n)
      }
      if (length(sigZ)!=n){
        stop("Z's parameters are not from the same length")
      }
      # support value
      support <- muX - sigX/ksiX
      # We impose the condition of equal support
      muZ <- support + sigZ/ksiZ
    }
    if (unknown == "scale"){
      # We take our sample size from the length of shape inputs
      n <- length(ksiZ)
      # We remind that if our shape inputs are just a value, the sample size will be 1
      if (n == 1){
        warning("Z's shape parameters length will give us the size of this strajectory")
      }
      # Verification. Shape parameters must be of the same sign
      for (i in 1:n){
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
        muZ <- rep(muZ,n)
      }
      if (length(muZ)!=n){
        stop("Z's parameters are not from the same length")
      }
      # support value
      support <- muX - sigX/ksiX
      # We impose the condition of equal support
      sigZ <- (support - muZ)*(-ksiZ)
    }
    if (unknown == "shape"){
      # We take our sample size from the length of scale inputs
      n <- length(sigZ)
      # We remind that if our shape inputs are just a value, the sample size will be 1
      if (n == 1){
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
        muZ <- rep(muZ,n)
      }
      if (length(muZ)!=n){
        stop("Z's parameters are not from the same length")
      }
      
      # support value
      support <- muX - sigX/ksiX
      # We impose the condition of equal support
      ksiZ <- sigZ/(muZ - support)
      # Verification. Shape parameters must be of the same sign
    }
    
      for (i in 1:n){
        ksiZ_i<- as.numeric(ksiZ[i])
        if ((sign(ksiZ_i)==sign(ksiX))==FALSE){
          stop("X and Z can not be W-class, shape parameters must have the same sign")
        }
      }
      
    matX<-matrix(nrow=m,ncol=N)
    matZ<-matrix(nrow=n,ncol=N)
    
    for(j in 1:N){
      matX[,j] <- rgev(m, loc=muX, scale= sigX, shape=ksiX)
      for(i in 1:n)
        matZ[i,j] <-rgev(1, loc=muZ[i], scale=sigZ[i], shape=ksiZ[i])}
    
    lam_t <- ((sigX/ksiX)/(sigZ/ksiZ))^(1/ksiX)
    k_t <- ksiX/ksiZ
  }
  
  if (ksiX==0){
    n <- length(sigZ)
    # We remind that if our shape inputs are just a value, the sample size will be 1
    # matX<-matrix(nrow=m,ncol=N)
    # matZ<-matrix(nrow=n,ncol=N)
    if (n == 1){
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
    ksiZ<-rep(0,n)
  
    # Verification that X's parameters are just values
    if (length(muX)!=1 | length(sigX)!=1 | length(ksiX)!=1){
      stop("X's parameters must be of just 1 value")
      }
    # Z's parameters must be of the same length
    if (length(muZ) == 1) {
      muZ <- rep(muZ,n)
      }
    if (length(sigZ) == 1) {
      sigZ <- rep(sigZ,n)
      }
    if (length(muZ)!=n){
      stop("Z's parameters are not from the same length")
      }
    if (length(sigZ)!=n){
      stop("Z's parameters are not from the same length")
    }
    
    k_t <- sigX/sigZ
    lam_t <- exp((muX-muZ)/sigX)
    
  }
  
  matX<-matrix(nrow=m,ncol=N)
  matZ<-matrix(nrow=n,ncol=N)
  
  for(j in 1:N){
    matX[,j] <- rgev(m, loc=muX, scale= sigX, shape=ksiX)
    for(i in 1:n){
      matZ[i,j] <-rgev(1, loc=muZ[i], scale=sigZ[i], shape=ksiZ[i])}
    
  }
  if (graph==TRUE){
    
    n <-length(ksiZ)

    muZ_t <- c(muZ[1],muZ[n/4],muZ[n/2], muZ[3*n/4], muZ[n])
    sigZ_t <- c(sigZ[1],sigZ[n/4],sigZ[n/2], sigZ[3*n/4],sigZ[n])
    ksiZ_t <- c(ksiZ[1],ksiZ[n/4],ksiZ[n/2], ksiZ[3*n/4], ksiZ[n])

    plot(density(matX[,1]),col="red",main= expression(paste(X[t]," and ",Z[t]," distribution evolution")))

    Z1 = rgev(n, loc = muZ_t[1], scale = sigZ_t[1], shape = ksiZ_t[1])
    lines(density(Z1),col="gray")
    Z2 = rgev(n, loc = muZ_t[2], scale = sigZ_t[2], shape = ksiZ_t[2])
    lines(density(Z2),col="green")
    Z3 = rgev(n, loc = muZ_t[3], scale = sigZ_t[3], shape = ksiZ_t[3])
    lines(density(Z3),col="blue")
    Z4 = rgev(n, loc = muZ_t[4], scale = sigZ_t[4], shape = ksiZ_t[4])
    lines(density(Z4),col="purple")
    legend(26.0,1.5,legend=c(expression(X[t]),expression(Z[t=1]),expression(Z[t=50]),expression(Z[t=100]),expression(Z[t=200])),col=c("red","gray","green","blue","purple"),lty=1:2,cex=0.7)
    
  }
  
  W=matrix(nrow=n, ncol=N)
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
  # size, mux, sigx and sigz must be scalar. # muz must be a vector
  X = rgev(size, loc = mux, scale = sigx, shape = 0)
  muz_t= c(muz[1],mu[length(muz)/4],mu[length(muz)/2], mu[3*length(muz)/4], mu[length(muz)])
  plot(density(X),col="black",)
  for (i in 1:length(muz_t)){
    Z = rgev(size, loc = muz_t[i], scale = sigz, shape = 0)
    lines(density(Z),col="gray")
  }
}

# Simulation of multiple stationnary {(X)t,(Z)t} W-class trajectories
simulWclass <- function(m,n,N,ksiX,ksiZ,sigX,supportauto=TRUE,muX=0,muZ=0,sigZ=0,graph=TRUE){
  
  if (sigZ == 0) { sigZ = sigX } 
  if (supportauto == TRUE){ muZ <- muX + sigZ/ksiZ - sigX/ksiX }
  matX <- matrix(rep(0,m*N),nrow=m,ncol=N) 
  matZ <- matrix(rep(0,n*N),nrow=n,ncol=N) 
  for (j in 1:N){
    matX[,j] = muX + (sigX/ksiX)*( (-log(runif(m)))^(-ksiX) - 1 )  
    matZ[,j] = muZ + (sigZ/ksiZ)*( (-log(runif(n)))^(-ksiZ) - 1 )  
  }
  if (graph==TRUE){
    plot(density(matX[,1]),
         main="Kernel density estimates \nfor the first X sample (black, counterfactual),\n and the first Z sample (red, factual)",
         cex.main=0.8,
         xlim=c(-4,12))
    lines(density(matZ[,1]),col="red")
  }
  
  # Calculation of lambda an k 
  return (list("matX"=matX,
               "matZ"=matZ,
               "lam"=((sigX/ksiX)/(sigZ/ksiZ))^(1/ksiX),
               "k"=ksiX/ksiZ,
               "mu_Z"=muZ,
               "sig_Z"=sigZ,
               "ksi_Z"=ksiZ))
  
}


################################################################################################
################################### CI for lambda_t and k_t ####################################
################################################################################################

##################### Variance - covariamce matrix of {p12_t, p13_t} ###########################
#####################                                                ###########################

#### VarCovar matrix ###########################################################################
matcovp12p13_t <- function(X.vec, Z.vec,tt,t_eval,h){
  
  ## computates varcovarp12p13_t = varcovarNA_t + varcovarNB_t
  
  ## Requires : matcovNA_t(), matcovNB_t()
  J <- length(Z.vec)
  I <- length(X.vec)
  
  list_varcovarN<- array(NA,c(2,2,J))
  
  G_emp <- ecdf(X.vec)
  GmZ <- G_emp(Z.vec)
  
  p13_hat.vec<- p13_NonPar(X.vec,Z.vec,tt,t_eval,h)
  p12_hat.vec<-p12_NonPar(X.vec,Z.vec,tt,t_eval,h)
  
  list_varcovarNA <- matcovNA_t(X.vec, Z.vec, GmZ,tt,t_eval,h) ##
  list_varcovarNB <- matcovNB_t(p12_hat.vec, p13_hat.vec, GmZ,tt,t_eval,h)
  
  for (i in 1:J){
    list_varcovarN[,,i]<- list_varcovarNB[,,i] + list_varcovarNA[,,i]
  }
  
  return(list_varcovarN)
}

##### varcovarNB_t #############################################################################
matcovNB_t <- function(I,J,p12_hat_t, p13_hat_t,GmZ,tt,t_eval,h){
  
  ## fluctuations de G_hat autour de G
  
  ## Requires : matcovNB_aij_t(), matcovNB_bij_t(), matcovNB_cij_t()
  
  theta_t <-weibullGMM_NonStationaire (GmZ, tt, t_eval, h, kern=dEpan, truevalues=NULL)
  lambda_t <- theta_t[[1]]
  k_t <- theta_t[[2]]
  
  Aij <- matcovNB_aij_t(p12_hat.vec,p13_hat.vec,lambda_t,k_t)
  Bij <- matcovNB_bij_t(p12_hat.vec,p13_hat.vec,lambda_t,k_t)
  Cij <- matcovNB_cij_t(p12_hat.vec,p13_hat.vec,lambda_t,k_t)
  
  Kh <- outer(t_eval, tt, function(zz,z) dEpan((zz - z) / h))
  denom<-rowSums(Kh)^2
  
  list_NB_t<- array(NA,c(2,2,J))
  
  for (i in 1:dim(Kh)[1]){
    Wij <- (h/(I*J)) * c(outer(Kh[i,], Kh[i,], "*")) 
    
    list_NB_t[1,1,i]<- sum(Wij*Aij) /denom[i]
    list_NB_t[1,2,i]<- sum(Wij*Cij) /denom[i]
    list_NB_t[2,1,i]<- sum(Wij*Cij) /denom[i]
    list_NB_t[2,2,i]<- sum(Wij*Bij) /denom[i]
  }
  return (list_NB_t) 
}

matcovNB_aij_t <- function(p12_hat_t,p13_hat_t,lambda_t,k_t){
  ## list of matcovNB_ij[1,1]
  
  liste_matcovNB_ij_11 <- c()
  
  for (j in 1:length(p12_hat_t)){
    lambda.j <- lambda_t[j]
    k.j <- k_t[j]
    p12_hat.j <- p12_hat_t[j]
    
    for (i in 1:length(p12_hat_t)){
      lambda.i <- lambda_t[i]
      k.i <- k_t[i]
      p12_hat.i <- p12_hat_t[i]
      
      NB_aij <-Mrfuncij(2,lambda.j,k.j,lambda.i,k.i) + p12_hat.j * p12_hat.i
      
      liste_matcovNB_ij_11 <- c(liste_matcovNB_ij_11, NB_aij)
    }
  }
  return (liste_matcovNB_ij_11)
}

matcovNB_cij_t <- function(p12_hat_t,p13_hat_t,lambda_t,k_t){
  ## list of matcovNB_ij[1,2] and matcovNB_ij[2,1]
  
  liste_matcovNB_ij_12 <- c()
  
  for (j in 1:length(p12_hat_t)){
    lambda.j <- lambda_t[j]
    k.j <- k_t[j]
    p12_hat.j <- p12_hat_t[j]
    p13_hat.j <- p13_hat_t[j]
    
    EGzjminGziGzj_partie1 <-partiedeEG1minG1G2func(lambda.j,k.j,lowerbnd=10^(-6),fac=0.5)
      
    for (i in 1:length(p12_hat_t)){
      lambda.i <- lambda_t[i]
      k.i <- k_t[i]
      p12_hat.i <- p12_hat_t[i]
      
      EGzjminGziGzj <- EGzjminGziGzj_partie1 + (1/2) * p12_hat.j * p12_hat.i + p13_hat.j
      NB_cij <- EGzjminGziGzj + (p13_hat.j * p12_hat.i)
    
      liste_matcovNB_ij_12 <- c(liste_matcovNB_ij_12, NB_cij)
    }
  }
  return (liste_matcovNB_ij_12)
}

matcovNB_bij_t <- function(p12_hat_t,p13_hat_t,lambda_t,k_t){
  ## list of matcovNB_ij[2,2]
  
  liste_matcovNB_ij_22 <- c()
  
  for (j in 1:length(p12_hat_t)){
    lambda.j <- lambda_t[j]
    k.j <- k_t[j]
    p13_hat.j <- p13_hat_t[j]
    
    for (i in 1:length(p12_hat_t)){
      lambda.i <- lambda_t[i]
      k.i <- k_t[i]
      p13_hat.i <- p13_hat_t[i]
      
      NB_bij <-Mrfuncij(3,lambda.j,k.j,lambda.i,k.i) + p13_hat.j * p13_hat.i
      
      liste_matcovNB_ij_22 <- c(liste_matcovNB_ij_22, NB_bij)
    }
  }
  return (liste_matcovNB_ij_22)
}


funcLaplaceterij <- function(v,m,lam.j,k.j,lam.i,k.i,a.j,a.i,lowerbnd=10^(-5),tol=10^(-5)){
  
  nv=length(v)
  mprime.j = (m-1)*lam.j / a.j^(1/k.j)
  mprime.i = (m-1)*lam.i / a.i^(1/k.i)
  facteur1.j = exp( -mprime.j*(-log(v))^(1/k.j) ) 
  facteur1.i = exp( -mprime.i*(-log(v))^(1/k.i) )
  facteur2.j=rep(0,nv)
  facteur2.i=rep(0,nv)
  for (i in 1:nv){
    facteur2.j[i] = laplaceWeibull( j=m, lambda=lam.j, k=k.j, lowerbnd=(10^-6)*v[i],
                                    upperbnd= v[i]^(1/a.j), tol=tol )  # borninf is guaranteed to be both close to 0, AND < v[i] 
  }
  for (i in 1:nv){
    facteur2.i[i] = laplaceWeibull( j=m, lambda=lam.i, k=k.i, lowerbnd=(10^-6)*v[i],
                                    upperbnd= v[i]^(1/a.i), tol=tol )  # borninf is guaranteed to be both close to 0, AND < v[i] 
  }
  facteur3.j=v^(1/a.j - 1)
  facteur3.i=v^(1/a.i - 1)
  return((1/a.j)*facteur1.j*facteur2.j*facteur3.j + (1/a.i)*facteur1.i*facteur2.i*facteur3.i)
}

Mrfuncij <- function(r,lambda.j,k.j,lambda.i,k.i,lowerbnd=10^(-5),fac=0.5,tol=10^(-5)){
  
  vala.j=fac*((r-1)*lambda.j)^k.j
  vala.i=fac*((r-1)*lambda.i)^k.i
  I <- integrate(f=funcLaplaceterij,
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

partiedeEG1minG1G2func <- function(lambda,k,lowerbnd=10^(-6),fac=0.5){
  # This function computes a part of E(G(Z1)*min(G(Z1),G(Z2))),
  # an expression appearing inside the limit covariance matrix of (\hat{p}_12,\hat{p}_13).
  # It is used inside function matcovNn().
  vala=fac*(2*lambda)^k
  I <- integrate(f=foncpartiedeEG1minG1G2,lower=lowerbnd,upper=1,subdivisions=1000L,
                 lam=lambda,k=k,a=vala,
                 stop.on.error = FALSE)
  return(I$value) 
}

foncpartiedeEG1minG1G2 <- function(x,lam,k,a){
  # Utilitary function used inside function partiedeRG1minG1G2func()
  (1/a) * x^( 2/a - 1) * exp( -(2*lam/a^(1/k)) * (-log(x))^(1/k) )
}

##### varcovarNA_t #############################################################################

matcovNA_t <- function(X.vec, Z.vec, GmZ,tt,t_eval,h){
  
  ## without G_hat 
  
  ## Requires : kern_dens(), K22(), p13_NonPar(),  p1r_NonPar()
  
  J <- length(Z.vec)
  
  f_t<-kern_dens(tt,t_eval,h) 
  
  list_varcovarNA <- array(NA,c(2,2,J))
  
  Aij <- matcovNA_aij_t(X.vec, Z.vec, GmZ,tt,t_eval,h) / f_t
  Bij <- matcovNA_bij_t(X.vec, Z.vec, GmZ,tt,t_eval,h) / f_t
  Cij <- matcovNA_cij_t(X.vec, Z.vec, GmZ,tt,t_eval,h) / f_t
  
  for (i in 1:J){
    list_varcovarNA[1,1,i]<-Aij[i]
    list_varcovarNA[1,2,i]<-Cij[i]
    list_varcovarNA[2,1,i]<-Cij[i]
    list_varcovarNA[2,2,i]<-Bij[i]
  }
  
  return (list_varcovarNA)
}

matcovNA_aij_t <- function(X.vec,Z.vec,tt,t_eval,h){
  ## matcovNA[1,1]_t
  J <- length(Z.vec)
  p13_hat<-p13_NonPar(X.vec,Z.vec,tt,t_eval,h)
  K22<-integrate(dEpan_2,-1,1)$value
  f_t<-kern_dens(tt,t_eval,h) #
  NA_aij_t <- (p13_hat * K22 * f_t)/(h*J)
  return (NA_aij_t)
}

matcovNA_bij_t <- function(X.vec,Z.vec,tt,t_eval,h){
  ## matcovNA[2,2]_t
  J <- length(Z.vec)
  p15_hat<-p1r_NonPar(X.vec,Z.vec,5,tt,t_eval,h)
  K22<-integrate(dEpan_2,-1,1)$value
  f_t<-kern_dens(tt,t_eval,h) #
  NA_bij_t <- (p15_hat * K22 * f_t)/(h*J)
  return (NA_bij_t)
}

matcovNA_cij_t <- function(X.vec,Z.vec, GmZ, tt,t_eval,h){
  ## matcovNA[1,2]_t & matcovNA[2,1]_t
  J <- length(Z.vec)
  p14_hat <- p1r_NonPar(X.vec,Z.vec,4,tt,t_eval,h)
  K22 <- integrate(dEpan_2,-1,1)$value
  f_t<-kern_dens(tt,t_eval,h) 
  Kij <- outer(t_eval,tt,function(zz,z) kern((zz - z) / h))
  r12 <- (Kij %*% GmZ)/J
  r13 <- (Kij %*% (GmZ^2))/J
  NA_cij_t <- ((p15_hat * K22 * f_t)/(h * J)) + (r12 * r13)
  return (NA_cij_t)
}

##### varcovarN_t #############################################################################

matcovN_t <- function(X.vec, Z.vec, GmZ,tt,t_eval.vec,h){
  J <- length(Z.vec)
  list_matcovN <- array(NA,c(2,2,J))
  list_matcovNA <- matcovNA_t(X.vec, Z.vec, GmZ,tt,t_eval,h)
  list_matcovNB <- matcovNB_t(X.vec, Z.vec,GmZ,tt,t_eval,h)
  for (i in 1:J){
    list_matcovN[,,i]<- list_matcovNA[,,i] + list_matcovNB[,,i]
  }
  return (list_matcovN)
}

##### varcovarp_t #############################################################################

matcovp_t <- function(matcovN.vec, X.vec, Z.vec, GmZ,tt,t_eval.vec,h){
  J <- legth(Z.vec)
  theta_t <-weibullGMM_NonStationaire (GmZ, tt, t_eval, h, kern=dEpan, truevalues=NULL)
  lam_t <- theta_t[[1]]
  k_t <- theta_t[[2]]
  Jacov <- jacobianFunctiong12(lam_t,k_t)
  list_matcovp <- array(NA,c(2,2,J))
  
  for (i in 1:J){
   # list_matcovp[,,i]<- inv( inv(jaconv[i]) %*% matcovN.vec[i] %*% jaconv[i] )
  }
  return (list_matcovp)
}


  
  