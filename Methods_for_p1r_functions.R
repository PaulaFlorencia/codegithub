##########################################################################################################################################
##########################################################################################################################################
library(gmm)
#############################################################   Useful functions  ########################################################
p1rWs_gamma <- function(r, lam, xi_G, xi_F, p_G){ # accepts a vector in r
  a <- abs(1/xi_F)
  b <- abs(xi_G)
  m <- r-1
  delta <- (-log(1-p_G))^(1/a)
  A <-  ( a*b*delta^(a-1) * exp(-(delta^a)) ) / ( (1-p_G)*(lam*m)^b ) 
  B <- ( a*b ) / ( (1-p_G)*(lam*m)^(2*b) )
  C <- ( (a-1) * (delta^(a-2)) * exp(-(delta^a)) ) - (  a * (delta^(2*a-2)) * exp(-(delta^a)) ) 
  prop2 <- A*gamma(b) + B*C*gamma(2*b)
  return(prop2)
}

p1r_NonP <- function(matGm,r.vec){
  r.vec <- as.numeric(r.vec)
  matGm <- as.matrix(matGm) 
  Ncol <- dim(matGm)[2] 
  nr <- length(r.vec)
  p1r.mat <- matrix(0,nrow=nr,ncol=Ncol)
  for (j in 1:Ncol){
    GnZ <- matGm[,j]
    p1r.vec <- rep(0,nr)
    for(i in 1:nr){
      p1r.vec[i] <- mean(GnZ^(r.vec[i] - 1))
    }
    p1r.mat[,j] <- p1r.vec
  }
  return(p1r.mat)
}

p1rfarW <- function(lam.vec,k.vec,r.vec,lowerbnd=10^(-5)){
  m=length(lam.vec)  
  p1r.vec=rep(0,m)
  for (j in 1:m){
    p1r.vec[j] <- laplaceWeibull(r.vec[j]-1,lam.vec[j],k.vec[j],lowerbnd=lowerbnd)
  }  
  far.vec <- 1 - 1 / (r.vec * p1r.vec)
  return( list("p1r"=p1r.vec,"far"=far.vec) )
} # VERIFY (parametric p1r under Wclass assumtion)

G <- function(u, xi_G, sig_G, mu_G){
  stopifnot(xi_G < 0)
  exp(-(1 + ksi_G * (u - mu_G)/sig_G)^(-1/ksi_G))
}

d_Ws <- function(u, lambda, p_G, xi_G, xi_F){
  #stopifnot(xi_G < 0 & xi_F < 0)
  delta <- (-log(1 - p_G))^abs(xi_F)
  abs(xi_G) / ((1 - p_G) * abs(xi_F) * lambda^abs(xi_G)) *
    exp( -((u / lambda)^abs(xi_G) + delta)^(1/abs(xi_F)) ) *
    ((u / lambda)^abs(xi_G) + delta)^(1/abs(xi_F) - 1) *
    u^(abs(xi_G) - 1)
}

p1r_star_integrand <- function(u, r, lambda, p_G, xi_G, xi_F){
  # stopifnot(xi_G < 0 & xi_F < 0)
  stopifnot(r >= 2)
  exp(-(r-1)*u) * d_Ws(u, lambda = lambda, p_G = p_G, xi_G = xi_G, xi_F = xi_F)
}

p1r_star <- function(r.vec, lam, p_G, xi_G, xi_F){
  m=length(r.vec)
  p1rstar.vec=rep(0,m)
  for (j in 1:m){
    # cat("r  = ", r.vec[j],", ")
    p1rstar.vec[j] <- Ws_integration(r.vec[j],lam, p_G, xi_G, xi_F)
  }  
  return(p1rstar.vec)
}

Ws_integration <- function(r, lambda, p_G, xi_G, xi_F,tol=1e-5){
  # cat("xi_F= ",xi_F,", xi_G= ",xi_G," lam= ",lambda,"\n")
  I <- integrate(f=p1r_star_integrand, lower=1e-8, upper=Inf,
                 r=r,
                 lambda=lambda,
                 p_G = p_G,
                 xi_G = xi_G,
                 xi_F = xi_F
  )
  resultat <- I$value
  return(resultat)
}

p1rWs_class <- function(r.vec, lambda, p_G, xi_G, xi_F){
  p1rstar.vec <- p1r_star(r.vec, lambda, p_G, xi_G, xi_F)
  p1rWs.vec <- p_G + (1-p_G)*p1rstar.vec
  return( list("p1rWs"=p1rWs.vec,"p1rstar"=p1rstar.vec) )
}

# A utilitary function required for the M-estimation in Ws
auxfunc <- function(theta, vecx, p_G){
  xiFval = theta[1]; xiGval = theta[2]; lambdaval = theta[3]
  p12val <- p1r_star(r=2, lam=lambdaval, p_G=p_G, xi_G=xiGval, xi_F=xiFval)
  p13val <- p1r_star(r=3, lam=lambdaval, p_G=p_G, xi_G=xiGval, xi_F=xiFval)
  p14val <- p1r_star(r=4, lam=lambdaval, p_G=p_G, xi_G=xiGval, xi_F=xiFval)
  M <- cbind( vecx - p12val ,  vecx^2 - p13val, vecx^3 - p14val )
  return(M)
}
# M-estimation
WsGMMestim <- function(matGm,p_G,truevalues=NULL){
  p_G <- as.numeric(p_G)
  matGm <- as.matrix(matGm) 
  Ncol <- dim(matGm)[2] 
  
  xiFhat.vec<- rep(0,Ncol) 
  xiGhat.vec<- rep(0,Ncol)
  lambdahat.vec <- rep(0,Ncol)
  
  if ( is.null(truevalues) ){
    startvalues <- matrix(c(-1,-1,1),nrow=3,ncol=Ncol)
    
  } else {
    startvalues <- truevalues
  }
  
  for (j in 1:Ncol){
    GnhatZ = matGm[,j] 
    p_G = p_G[j]
    fg <- function(theta,vecx){auxfunc(theta,vecx,p_G)}
    EGMMWs<- gmm(g=fg , 
                 x=GnhatZ, 
                 t0=startvalues[,j],
                 optfct = "nlminb",
                 lower=c(-Inf,-Inf,0),upper=c(-1e-10,-1e-10,Inf),
                 onlyCoefficients = TRUE
    )
    xiFhat.vec[j] <- EGMMWs$coefficients[[1]] 
    xiGhat.vec[j] <-EGMMWs$coefficients[[2]] 
    lambdahat.vec[j] <- EGMMWs$coefficients[[3]] 
  }   
  return( list("xiFhat"=xiFhat.vec,"xiGhat"=xiGhat.vec,"lambdahat"=lambdahat.vec,"p_Ghat"=p_G) )
}

##########################################################   Auxiliary functions  ########################################################
# This functions are the same that the ones used on record's paper
laplaceWeibull <- function(j,lambda,k,lowerbnd=10^(-6),upperbnd=1,fac=1,tol=10^(-5)){
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

funcLaplace <- function(x,m,lam,k,a){ 
  (1/a) * exp( -(m*lam/a^(1/k))*(-log(x))^(1/k) ) * x^(1/a - 1) 
}

functiong <- function(theta,vecx){
  lambdaval=theta[1] ; kval = theta[2]
  p12val <- laplaceWeibull(j=1,lambda=lambdaval,k=kval,lowerbnd=10^(-8))
  p13val <- laplaceWeibull(j=2,lambda=lambdaval,k=kval,lowerbnd=10^(-8))
  M <- cbind( vecx - p12val ,  vecx^2 - p13val )
  return(M)
}

GZestimation <- function(matX,matZ,methodGhatm="ecdfmodif",bandwth=0,b=0.05){
  matX <- as.matrix(matX)
  matZ <- as.matrix(matZ) 
  dimnsZ=dim(matZ) 
  matGm <- matrix(rep(0,prod(dimnsZ)),nrow=dimnsZ[1],ncol=dimnsZ[2])
  for (j in 1:dimnsZ[2]){
    X <- matX[,j]
    Z <- matZ[,j]
    if (methodGhatm=="npudist"){
      if (bandwth == 0){
        matGm[,j] <- npudist(tdat=X,edat=Z)$dist  
      } else {
        matGm[,j] <- npudist(tdat=X,edat=Z,bws=bandwth)$dist  
      }
    } 
    if (methodGhatm=="ecdf"){
      matGm[,j] <- ecdf(X)(Z) 
    }
    if (methodGhatm=="ecdfmodif"){
      matGm_temp <- ecdf(X)(Z) 
      m <- nrow(matX)
      matGm[,j] <- (m*matGm_temp+b)/(m+1) 
    }
  }
  return(matGm)  
}

weibullGMMestim <- function(matGm,truevalues=NULL){
  matGm <- as.matrix(matGm) 
  Ncol <- dim(matGm)[2] 
  lambdahat.vec <- rep(0,Ncol)
  khat.vec <- rep(0,Ncol)
  if ( is.null(truevalues) ){
    startvalues <- matrix(rep(1,2*Ncol),nrow=2,ncol=Ncol)
    # starting values, if not precised, are set to 1 and 1 for lambda and k
  } else {
    startvalues <- truevalues
  }
  for (j in 1:Ncol){
    #cat(j,":")
    GnhatZ = matGm[,j]  
    EGMMweibull <- gmm(g=functiong , 
                       x=GnhatZ , 
                       t0=startvalues[,j],
                       optfct = "nlminb",
                       lower=c(10^(-8),10^(-8)),upper=c(Inf,Inf),
                       onlyCoefficients = TRUE
    )
    
    lambdahat.vec[j] <- EGMMweibull$coefficients[[1]] 
    khat.vec[j] <- EGMMweibull$coefficients[[2]] 
    #cat(" ",round(lambdahat.vec[j],4),",")
    #cat(round(khat.vec[j],4),"(GMM)\n")
  }   
  return( list("lambdahat"=lambdahat.vec,"khat"=khat.vec) ) 
}

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
  return(list("matX"=matX,
              "matZ"=matZ,
              "lam"=((sigZ/ksiZ)/(sigX/ksiX))^(1/ksiX),
              "k"=ksiX/ksiZ,
              "muZ"=muZ,
              "muX"=muX,
              "ksiZ"=ksiZ,
              "ksiX"=ksiX,
              "sigZ"=sigZ,
              "sigX"=sigX
  ) 
  )  
} 
