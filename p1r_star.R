#########################  p1r estimation  ######################################
######################  ( p_G & p1r* method)  ###################################

# FUNCTIONS AT THE END OF THE SCRIPT 
# section 1 : data 
# section 2: parameters estimation
# section 3: record probability estimation
# section 4: summary of results

######################  Tests  ##################################################

### 1. samples examples ###
J <- 30; I <- 150; N <- 1
# we must chose one of the following siumlations
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.2,sigX=1,sigZ=1,
                    supportauto=FALSE,muX=0,muZ=3)

data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.2,sigX=1,sigZ=1,
                    supportauto=FALSE,muX=0,muZ=4)

data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.5,sigX=1,sigZ=1,
                    supportauto=FALSE,muX=0,muZ=3.5)

data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.5,sigX=1,sigZ=1,
                    supportauto=FALSE,muX=0,muZ=4)

data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.6,sigX=1,sigZ=1.5,
                    supportauto=FALSE,muX=0,muZ=4)

data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.6,sigX=1,sigZ=1.5,
                    supportauto=FALSE,muX=0,muZ=3.5)

data <- simulWclass(m=I,n=J,N=N,ksiX=-.2,ksiZ=-.25,sigX=1,sigZ=1,
                    supportauto=FALSE,muX=0,muZ=2)

data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.01,sigX=1.2,sigZ=1.5,
                    supportauto=FALSE,muX=4,muZ=5.7)

data <- simulWclass(m=I,n=J,N=N,ksiX=-0.8,ksiZ=-0.8,sigX=1,sigZ=1,
                    supportauto=FALSE,muX=0,muZ=1)

data <- simulWclass(m=I,n=J,N=N,ksiX=-1,ksiZ=-1,sigX=1,sigZ=1,##
                    supportauto=FALSE,muX=0,muZ=1)

Z <- data$matZ; X <- data$matX

### 2. theoretic and estimated parameters ###
# theoretic values of the parameters ksi,sigma, mu, lambda and k of G,F and W 
lam <- data$lam; k <- data$k; ksi_G <- data$ksiX; ksi_F <- data$ksiZ # (from here we can know the theoretic parameters of Ws)
mu_G <- data$muX; mu_F <- data$muZ; sig_G <- data$sigX; sig_F <- data$sigZ
# theoretic limit points of G and F (x_G and z_F)
x_G <- mu_G - sig_G/ksi_G; z_F <- mu_F - sig_F/ksi_F
# actual limit points of our sample (x_max and z_max)
x_max <- max(X); z_max <- max(Z)
# subsamble of z, values that are smaller than x_G
Zs <- Z[ Z < x_G]
# theoretic size of the subsample
NG <- length(Zs)
# p_G = P(Z>x_G)
p_G <- (J-NG)/J
# subsamble of z, values that are smaller than x_max
Zshat <- Z[ Z < x_max]
# actual size of the subsample
Nstar <- length(Zshat)
# G(Z) using theoretic CDF
matGZ_cdf <- CDF_G(Z, ksi_G, sig_G, mu_G)
# G(Zs) using theoretic CDF
matGZs_cdf <- CDF_G(Zs, ksi_G, sig_G, mu_G)
# G(Z) using empirical CDF
matGZ_ecdf <- GZestimation(X,Z)
# G(Zs) using empirical CDF
matGZs_ecdf <- GZestimation(X,Zs)
# G(Zshat) using empirical CDF
matGZshat_ecdf <- GZestimation(X,Zshat)
# p_Ghat, estimation of p_max = P(Z>x_max)
p_Ghat <- (J-Nstar)/J
# estimated parameters of Ws using CDF of G and the THEORETIC p_G
thetaWs_pG_cdf <- WsGMMestim(matGZs_cdf,p_G,truevalues=NULL)
ksi_F_pG_cdf <- thetaWs_pG_cdf$xiFhat; ksi_G_pG_cdf <- thetaWs_pG_cdf$xiGhat; lam_pG_cdf <- thetaWs_pG_cdf$lambdahat
# estimated parameters of Ws using ECDF of G and the THEORETIC p_G
thetaWs_pG_ecdf<-WsGMMestim(matGZs_ecdf,p_G,truevalues=NULL)
ksi_F_pG_ecdf <- thetaWs_pG_ecdf$xiFhat; ksi_G_pG_ecdf <- thetaWs_pG_ecdf$xiGhat; lam_pG_ecdf <- thetaWs_pG_ecdf$lambdahat
# estimated parameters of Ws using ECDF of G and the ESTIMATOR p_Ghat
thetaWs_pGhat_ecdf<-WsGMMestim(matGZshat_ecdf,p_Ghat,truevalues=NULL)
ksi_F_pGhat_ecdf <- thetaWs_pGhat_ecdf$xiFhat; ksi_G_pGhat_ecdf <- thetaWs_pGhat_ecdf$xiGhat; lam_pGhat_ecdf <- thetaWs_pGhat_ecdf$lambdahat

### 3. p1r estimation ###
r = 100 # record length choice
# p1r under Wclass assumption, using theoretic parametersn
p1rWclass <- p1rfarW(lam,k,r)$p1r

# p1r estimation under WCLASS assumption (using ECDF of G)
Wclassthetahat_ecdf <- weibullGMMestim(matGZ_ecdf)
lamW_hat_ecdf <- Wclassthetahat_ecdf$lambdahat; kW_hat_ecdf <- Wclassthetahat_ecdf$khat
p1rWclass_hat_ecdf <- p1rfarW(lamW_hat_ecdf,kW_hat_ecdf,r)$p1r

# p1r estimation under WCLASS assumption (using CDF of G)
Wclassthetahat_cdf <- weibullGMMestim(matGZ_cdf)
lamW_hat_cdf <- Wclassthetahat_cdf$lambdahat; kW_hat_cdf <- Wclassthetahat_cdf$khat
p1rWclass_hat_cdf <- p1rfarW(lamW_hat_cdf,kW_hat_cdf,r)$p1r

# p1r estimator using non-parametric method (using ecdf of G)
p1r_nonpar <- p1rfar_NPestim(matGZ_ecdf,r.vec=rep(r,N))$p1rhat

# p1r estimated as p1r_tilde ( p1r_tilde = p_G + (1-p_G)*p1rt_star ), using THEORETIC p_G and CDF of G
p1r_starhat_pG_cdf <- p1r_star(r, lam_pG_cdf, p_G, ksi_G_pG_cdf, ksi_F_pG_cdf )
p1r_tildehat_pG_cdf <- p_G + (1-p_G)*p1r_starhat_pG_cdf

# p1r estimated as p1r_tilde ( p1r_tilde = p_G + (1-p_G)*p1rt_star ), using THEORETIC p_G and ECDF of G
p1r_starhat_pG_ecdf <- p1r_star(r, lam_pG_ecdf, p_G, ksi_G_pG_ecdf, ksi_F_pG_ecdf )
p1r_tildehat_pG_ecdf <- p_G + (1-p_G)*p1r_starhat_pG_ecdf

# p1r estimated as p1r_tilde ( p1r_tilde = p_G + (1-p_G)*p1rt_star ), using ESTIMATED p_G and ECDF of G
p1r_starhat_pGhat_ecdf <- p1r_star(r, lam_pGhat_ecdf, p_Ghat, ksi_G_pGhat_ecdf, ksi_F_pGhat_ecdf )
p1r_tildehat_pGhat_ecdf <- p_Ghat + (1-p_Ghat)*p1r_starhat_pGhat_ecdf

# p1r computed as p1r_tilde ( p1r_tilde = p_G + (1-p_G)*p1rt_star ), using theoretic parameters
p1r_star_theo <- p1r_star(r, lam, p_G, ksi_G, ksi_F)
p1r_tilde_theo <- p_G + (1-p_G)*p1r_star_theo

### 4. Summary of results ###
cat("Theoretical parameters:","\n",
    "(theo) paramaters of G:","\n", "ksi_G = ",ksi_G, " sig_G = ",sig_G," mu_G = ", mu_G,"\n",
    "(theo) parameters of F","\n", "ksi_F = ",ksi_F, " sig_F = ",sig_F," mu_F = ", mu_F,"\n",
    "(theo) parameters of W:","\n","lam = ", lam, " k = ",k,"\n",
    "\n",
    "points limites:","\n",
    "(theo) x_G = ", x_G, " (estimated) x_max  = ", x_max,"\n" ,
    "(theo) z_F = ",z_F, " (estimated) z_max = ", z_max,"\n",
    "\n",
    "sample size for computing Ws:","\n",
    "(theo) NG = ",NG, ", part of Z = ", J/NG,"\n",
    "(estimated) Nstar = ",Nstar, ", part of Z = ", J/Nstar,"\n",
    "\n",
    "p_G thoretic and estimated:","\n",
    "(theo) p_G = P(Z>x_G) = ",p_G,"\n",
    "(estimated) p_Ghat = P(Z>x_max) = ",p_Ghat,"\n",
    "\n",
    "W (lam,k) theoretic and estimated parameters under Wclass assumption (we are not in Wclass) :","\n",
    "(theo) lam = ",lam," k = ",k,"\n",
    " (estimated under Wclass asumption using cdf) lamW_hat = ",lamW_hat_cdf," k = ", kW_hat_cdf, "\n",
    " (estimated under Wclass asumption using ecdf) lamW_hat = ",lamW_hat_ecdf," k = ", kW_hat_ecdf, "\n",
    "\n",
    "Ws (ksi_F,ksi_G,lam) theoretic and estimated parameters :","\n",
    "(theo) ksi_G = ",ksi_G, " ksi_F = ",ksi_F," lam = ", lam,"\n",
    "(estimated by MMO using p_G and cdf) ksi_F_pG = ",ksi_F_pG_cdf," ksi_G_pG = ",ksi_G_pG_cdf," lam_pG = ",
    lam_pG_cdf,"\n",
    "(estimated by MMO using p_G and ecdf) ksi_F_pG = ",ksi_F_pG_ecdf," ksi_G_pG = ",ksi_G_pG_ecdf," lam_pG = ",
    lam_pG_ecdf,"\n",
    "\n",
    "(estimated by MMO using p_Ghat and ecdf) ksi_F_pGhat = ",ksi_F_pGhat_ecdf," ksi_G_pGhat = ",ksi_G_pGhat_ecdf,
    " lam_pGhat = ",lam_pGhat_ecdf,"\n",
    "\n",
    "p1r computed as p1r_tilde ( p1r_tilde = p_G + (1-p_G)*p1rt_star ), using theoretic parameters","\n",
    "p1r_tilde = ",p1r_tilde_theo, " p1r* = ",p1r_star_theo ,"\n",
    "\n",
    "p1r estimated as p1r_tilde ( p1r_tilde = p_G + (1-p_G)*p1rt_star ), using estimated p_G :","\n",
    "(p_Ghat & ecdf) p1r = p1r_tilde_hat = ", p1r_tildehat_pGhat_ecdf," p1r*_hat = ", p1r_starhat_pGhat_ecdf, " p_G = ", p_Ghat,"\n",
    "\n",
    "p1r estimated as p1r_tilde ( p1r_tilde = p_G + (1-p_G)*p1rt_star ), using theoretic p_G : ","\n",
    "(p_G & ecdf) p1r = p1r_tilde_hat = ", p1r_tildehat_pG_ecdf," p1r*_hat = ", p1r_starhat_pG_ecdf, " p_G = ", p_G,"\n",
    "\n",
    "(p_G & cdf) p1r = p1r_tilde_hat = ", p1r_tildehat_pG_cdf," p1r*_hat = ", p1r_starhat_pG_cdf, " p_G = ", p_G,"\n",
    "\n",
    "p1r estim non-parametric method : ","\n",
    "p1r = p1r_nonpar = ",p1r_nonpar, "\n",
    "\n",
    "p1r estimation under Wclass (we are not in Wclass) : ","\n",
    "(ecdf) p1r_hat = p1rWclass_hat = ",p1rWclass_hat_ecdf,"lam = ",lamW_hat_ecdf,"k = " , kW_hat_ecdf,"\n",
    "(cdf) p1r_hat = p1rWclass_hat = ",p1rWclass_hat_cdf,"lam = ",lamW_hat_cdf,"k = " , kW_hat_cdf, "\n",
    "\n",
    "Theoretic p1r under Wclass assumption (we are not in Wclass), using theoretic lam & k : ","\n",
    "p1r = p1rW = ", p1rWclass ,"\n",
    "\n",
    "Summary of estimations p1r using cdf: ","\n",
    "(p_G) p1r = p1r_tilde_hat = ", p1r_tildehat_pG_cdf,"\n",
    "p1r_hat = p1rWclass_hat = ",p1rWclass_hat_cdf, "(we are not in Wclass) ", "\n",
    "\n",
    "Summary of estimations p1r using ecdf: ","\n",
    "(p_Ghat) p1r_hat = p1r_tilde_hat = ", p1r_tildehat_pGhat_ecdf,"\n",
    "(p_G) p1r_hat = p1r_tilde_hat = ", p1r_tildehat_pG_ecdf,"\n",
    "p1r_hat = p1rWclass_hat = ",p1rWclass_hat_ecdf, "(we are not in Wclass) ", "\n",
    "\n",
    "Non parametric estimation", "\n",
    "p1r_hat = p1r_nonpar = ",p1r_nonpar, "\n",
    "\n",
    "Theoretic p1r", "\n",
    "p1r = p1r_tilde = ",p1r_tilde_theo, "\n",
    "p1r = p1rW = ",p1rWclass, "(we are not in Wclass) ","\n"
    
    )

######################  Functions  ##############################################
# CDF G
CDF_G <- function(u, xi_G, sig_G, mu_G){
  stopifnot(xi_G < 0 )
  s <- (u - mu_G)/sig_G
  exp ( -( 1+xi_G*s )^(-1/(xi_G))
        )
}
CDF_G_plot <- function(x) CDF_G(x, ksi_G, sig_G, mu_G)
# density  
d_Ws <- function(u, lambda, p_G, xi_G, xi_F){
  #stopifnot(xi_G < 0 & xi_F < 0)
  delta <- (-log(p_G))^abs(xi_F)
  abs(xi_G) / (p_G * abs(xi_F) * lambda^abs(xi_G)) *
    exp( -((u / lambda)^abs(xi_G) + delta)^(1/abs(xi_F)) ) *
    ((u / lambda)^abs(xi_G) + delta)^(1/abs(xi_F) - 1) *
    u^(abs(xi_G) - 1)
}
# Integrand in p1r_star
p1r_star_integrand <- function(u, r, lambda, p_G, xi_G, xi_F){
  #stopifnot(xi_G < 0 & xi_F < 0)
  stopifnot(r >= 2)
  delta <- (-log(p_G))^abs(xi_F)
  exp(-(r-1)*u) * d_Ws(u, lambda = lambda, p_G = p_G, xi_G = xi_G, xi_F = xi_F)
}
# p1r_star
p1r_star <- function(r, lambda, p_G, xi_G, xi_F,tol=10^(-5)){
  cat("xi_F= ",xi_F,", xi_G= ",xi_G," lam= ",lambda,"\n")
  I <- integrate(f=p1r_star_integrand, lower=0, upper=Inf,
                 r=r,
                 lambda=lambda,
                 p_G = p_G,
                 xi_G = xi_G,
                 xi_F = xi_F
  )
  resultat <- I$value
  return(resultat)
}
# A utilitary function required for the M-estimation 
auxfunc <- function(theta, vecx, p_G){
  xiFval = theta[1]; xiGval = theta[2]; lambdaval = theta[3]
  p12val <- p1r_star(r=2, lambda=lambdaval, p_G=p_G, xi_G=xiGval, xi_F=xiFval)
  p13val <- p1r_star(r=3, lambda=lambdaval, p_G=p_G, xi_G=xiGval, xi_F=xiFval)
  p14val <- p1r_star(r=4, lambda=lambdaval, p_G=p_G, xi_G=xiGval, xi_F=xiFval)
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
    fg <- function(theta,vecx){auxfunc(theta,vecx,p_G=p_G)}
    EGMMWs<- gmm(g=fg , 
                 x=GnhatZ, 
                 t0=startvalues[,j],
                 optfct = "nlminb",
                 lower=c(-Inf,-Inf,0),upper=c(0,0,Inf),
                 onlyCoefficients = TRUE
    )
    xiFhat.vec[j] <- EGMMWs$coefficients[[1]] 
    xiGhat.vec[j] <-EGMMWs$coefficients[[2]] 
    lambdahat.vec[j] <- EGMMWs$coefficients[[3]] 
  }   
  return( list("xiFhat"=xiFhat.vec,"xiGhat"=xiGhat.vec,"lambdahat"=lambdahat.vec,"p_Ghat"=p_G) )
}
# p1r paramametric estimation under Wclass assumption (from records paper)
p1rfarW <- function(lam.vec,k.vec,r.vec,lowerbnd=10^(-5)){
  m=length(lam.vec)  
  p1r.vec=rep(0,m)
  for (j in 1:m){
    p1r.vec[j] <- laplaceWeibull(r.vec[j]-1,lam.vec[j],k.vec[j],lowerbnd=lowerbnd)
  }  
  far.vec <- 1 - 1 / (r.vec * p1r.vec)
  return( list("p1r"=p1r.vec,"far"=far.vec) )
}

### Function from records paper ###
p1rfar_NPestim <- function(matGm,r.vec,a=1){
  matGm <- as.matrix(matGm) 
  Ncol <- dim(matGm)[2] 
  nn <- dim(matGm)[1] 
  p1r <- rep(0,Ncol)
  p12rm1 <- p1r ; Mr <- p1r ; tau2r <- p1r ; sigma2r <- p1r
  
  mat.r <- matrix(rep(r.vec,nn),nrow=nn,ncol=Ncol,byrow=TRUE) 
  
  Grmoins2 = matGm^(mat.r - 2) 
  Grmoins1 = Grmoins2 * matGm 
  p1r = apply(Grmoins1,2,mean)  
  p12rm1 = apply(Grmoins1*Grmoins1,2,mean) 
  tau2r = p12rm1 - p1r^2 
  
  for (j in 1:Ncol){
    GnZ <- matGm[,j]
    Grm2 <- Grmoins2[,j]
    productsGrm2 <- outer(Grm2,Grm2)    # by default, the product is applied
    # and productsGrm2 is thus a nn x nn matrix, containing the Gn(Zi)^(r-2)*Gn(Zl)^(r-2)
    minimaG <- outer(GnZ,GnZ,pmin)  # contains the min(Gn(Zi),Gn(Zl)), and is also nn x nn
    # (rem : pmin is indeed required, not simply min !) 
    matMr <- productsGrm2 * minimaG
    Mr[j] <- mean(matMr)  
    # this is indeed the mean of the nn^2 values of the matrix matMr
  }
  
  sigma2r = tau2r + a*(r.vec-1)^2 * (Mr - p1r^2)  
  # here comes up the value a, in practice the ratio sqrt(n/m)
  
  return( list("p1rhat"=p1r,"p12rm1hat"=p12rm1,"tau2rhat"=tau2r,
               "Mrhat"=Mr,"sigma2rhat"=sigma2r) )
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

### Auxiliary functions from records paper ###
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

