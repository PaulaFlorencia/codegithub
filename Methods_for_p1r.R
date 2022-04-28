source("Methods_for_p1r_functions.R")
# library(gmm)

##########################################################################################################################################
##########################################################################################################################################
# Create data 

J <- 5000; I <- 5000; N <- 1

# Same support
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.2,sigX=1,sigZ=1, # TESTER SI CET EXAMPLE CREER DES TRAJ AVEC p_G PETIT (FAUT VERIFIER), TESTER DES AUTRES AUSSI
                    supportauto=TRUE,muX=0) # p_G = P(Z>x_G) =  0
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.3,sigX=1,sigZ=1.1,
                    supportauto=TRUE,muX=0)
# Small p_G
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.2,sigX=1,sigZ=1, #ca plante quand pG_hat est < 0.5 , meme si en realité pg=0.07
                    supportauto=FALSE,muX=0,muZ=3) #
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.2,sigX=1,sigZ=1, 
                    supportauto=FALSE,muX=0,muZ=1.5) 
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.2,sigX=1,sigZ=1.2,
                    supportauto=FALSE,muX=0,muZ=0) 
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.3,sigX=1,sigZ=1.2,
                    supportauto=FALSE,muX=0,muZ=2.5)
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.3,sigX=1,sigZ=1.2,
                    supportauto=FALSE,muX=0,muZ=3)
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.25,ksiZ=-0.1,sigX=1,sigZ=1,
                    supportauto=FALSE,muX=0,muZ=1)
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.9,ksiZ=-0.9,sigX=1,sigZ=1, # test. eliminer après
                    supportauto=FALSE,muX=0,muZ=3)
# Not that small p_G
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.5,sigX=1,sigZ=1,
                    supportauto=FALSE,muX=0,muZ=4) # p_G = P(Z>x_G) =  0.227
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.6,sigX=1,sigZ=1.5,
                    supportauto=FALSE,muX=0,muZ=4) # p_G = P(Z>x_G) =  0.332
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.2,sigX=1,sigZ=1,
                    supportauto=FALSE,muX=0,muZ=3.5)
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.1,ksiZ=-0.9,sigX=1,sigZ=1,# test. eliminer après
                    supportauto=FALSE,muX=0,muZ=3)
# Big p_G
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.2,sigX=1,sigZ=1, 
                    supportauto=FALSE,muX=0,muZ=5) # p_G = P(Z>x_G) =  0.622
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.8,ksiZ=-0.8,sigX=1,sigZ=1,
                    supportauto=FALSE,muX=0,muZ=2)
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.4,ksiZ=-0.4,sigX=1,sigZ=1.1,
                    supportauto=FALSE,muX=0,muZ=3.5)
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.2,ksiZ=-0.2,sigX=1,sigZ=1, 
                    supportauto=FALSE,muX=0,muZ=6) 
data <- simulWclass(m=I,n=J,N=N,ksiX=-0.1,ksiZ=-0.2,sigX=1,sigZ=1, # test. eliminer après
                    supportauto=FALSE,muX=0,muZ=6)
data <- simulWclass(m=I,n=J,N=N,ksiX=-1,ksiZ=-0.6,sigX=1,sigZ=1, # test. eliminer après
                    supportauto=FALSE,muX=0,muZ=2)

Z <- data$matZ; X <- data$matX # trajectories chosen 

##########################################################################################################################################
##########################################################################################################################################
################################## 1. Case where we know the parameter of both distributions  ############################################

# Here we know the real xi, sig, mu (then also k , lam & p_G)  of our trajectories
lam <- data$lam; k <- data$k; ksi_G <- data$ksiX; ksi_F <- data$ksiZ
mu_G <- data$muX; mu_F <- data$muZ; sig_G <- data$sigX; sig_F <- data$sigZ
# length of record 
r <- 100; rvec <- c(10:100) 

### Using Wclass method  
########################################################################################################################################## 
p1rWclass <- p1rfarW(lam,k,r)$p1r 
cat("for r=100 and under Wclass assumtion, record probability is : ", p1rWclass) 
p1rWclass.vec <- p1rfarW(rep(lam,length(rvec)),rep(k,length(rvec)),rvec)$p1r 

### Using Ws method ( p1r = p_G + (1-p_G)*p1r_star )
########################################################################################################################################## 
x_G <- mu_G - sig_G/ksi_G; z_F <- mu_F - sig_F/ksi_F
Zs <- Z[ Z < x_G] 
Ns <- length(Zs) 
cdfF <- function(x) exp(-(1 + ksi_F * (x - mu_F)/sig_F)^(-1/ksi_F))
p_Gbar <- cdfF(x_G)
p_G <- 1-p_Gbar
cat("p_G = P(Z>x_G) = ", p_G) 
p1rWsclass_andstar <- p1rWs_class(r, lam, p_G, ksi_G, ksi_F)
p1rWsclass <- p1rWsclass_andstar$p1rWs
p1rstar <- p1rWsclass_andstar$p1rstar
cat("for r=100 and using Ws_class method, record probability is: ", p1rWsclass) 
p1rWs_prop2<-p1rWs_gamma(r, lam,ksi_G, ksi_F, p_G)
cat("pour r =  100 : p1rstar = ",p1rstar," - p1rstar_prop2 = ", p1rWs_prop2 ,"\n")
p1rWsclass_andstar.vec <- p1rWs_class(as.numeric(rvec), lam, p_G, ksi_G, ksi_F)  # Same but for different r
p1rWsclass.vec <- p1rWsclass_andstar.vec$p1rWs
p1rstar.vec <- p1rWsclass_andstar.vec$p1rstar
p1rWs_prop2.vec<-p1rWs_gamma(as.numeric(rvec), lam,ksi_G, ksi_F, p_G)

### Non-parametric estimation method
########################################################################################################################################## 
matGZ_ecdf <- GZestimation(X,Z)
p1r_nonpar <- p1r_NonP(matGZ_ecdf,r)
cat("for r=100 the estimated record probability using non-parametric method is : ", p1r_nonpar) 
p1r_nonpar.vec <- p1r_NonP(matGZ_ecdf,rvec)

### Comparison of three methods 
########################################################################################################################################## 
cat("Using r = 100","\n","p1rW: ",p1rWclass," - p1rWs: ",p1rWsclass," - p1r_nonParam: ",p1r_nonpar,"\n",", where  p_G = ", p_G)
plot(p1r_nonpar.vec,type="l",ylim=c(0,1),xlab = "r",ylab="p1r",main = " theo and nonpar p1r") # plot of p1r, for r= 10, 20, 30,..., 100
lines(p1rWclass.vec,col="blue")
lines(p1rWsclass.vec,col="red")
legend("topleft",legend=c("p1rW",  "p1rWs", "p1rNonPar"), col=c("blue", "red", "black"),lty=c(1,1,1),cex=0.8)


# Comparison between p1rstar (theo) and prop2 pour p1rstar (theo)
########################################################################################################################################## 
cat("pour r =  100 : p1rstar = ",p1rstar," - p1rstar_prop2 = ", p1rWs_prop2 ,"\n")
plot(p1rstar.vec,type="l",ylim=c(0,max(p1rstar.vec)),xlab="r",ylab="p1r_star",main = " p1r_star & prop2") # plot of p1r, for r= 10, 20, 30,..., 100
lines(p1rWs_prop2.vec,col="green")
legend("topleft",legend=c("prop2",  "p1r_star"),col=c("blue","black"),lty=c(1,1),cex=0.8)

##########################################################################################################################################
##########################################################################################################################################
############### 2. Case where we have just the trajectories and we know nothing about it's parameters (real context) #####################

### Estimation under Wclass assumtion
########################################################################################################################################## 
matGZ_ecdf <- GZestimation(X,Z)
Wclass_theta <- weibullGMMestim(matGZ_ecdf) # we use the ecdf computated before
lamW_hat <- Wclass_theta$lambdahat; 
kW_hat <- Wclass_theta$khat
p1rWclass_hat <- p1rfarW(lamW_hat,kW_hat,r)$p1r
cat("for r=100 the record probability using Wclass parametric method: ", p1rWclass_hat ," (unknown parameters)","\n") 
p1rWclass_hat.vec <- p1rfarW(rep(lamW_hat,length(rvec)),rep(kW_hat,length(rvec)),rvec)$p1r # for each r between 10 and 100. We will plot it later (in a comparision plot)

### Estimation using Ws method
########################################################################################################################################## 
x_max <- max(X) # Here we no longer know the actual support su we estimate it as max(x)
Zs_hat <- Z[ Z < x_max]
Ns_hat <- length(Zs_hat)
N_G_hat <-J-Ns_hat
p_G_hat <- (N_G_hat)/J
cat("p_G_hat = P(Z>x_max) = ", p_G_hat ) 
matGZs_ecdf <- GZestimation(X,Zs_hat) 
Wsclass_theta <- WsGMMestim(matGZs_ecdf,p_G_hat,truevalues=NULL) 
#######################
# p1rstarintegarnd <- p1r_star_integrand_visu(r.vec, lam, p_G, xi_G, xi_F)
#######################
ksiFWs_hat <- Wsclass_theta$xiFhat;
ksiGWs_hat <- Wsclass_theta$xiGhat; 
lamWs_hat <- Wsclass_theta$lambdahat
p1rWsclass_andstar_hat <- p1rWs_class(r, lamWs_hat, p_G_hat, ksiGWs_hat, ksiFWs_hat)
p1rWsclass_hat <- p1rWsclass_andstar_hat$p1rWs
p1rstar_hat <- p1rWsclass_andstar_hat$p1rstar
cat("for r=100 the record probability using Ws parametric method: ", p1rWsclass_hat," (unknown parameters)","\n") 
p1rWsclass_andstar_hat.vec <- p1rWs_class(rvec, lamWs_hat, p_G_hat, ksiGWs_hat, ksiFWs_hat)
p1rWsclass_hat.vec <- p1rWsclass_andstar_hat.vec$p1rWs  #
p1rstar_hat.vec <- p1rWsclass_andstar_hat.vec$p1rstar
p1rWs_prop2_hat <- p1rWs_gamma(r, lamWs_hat,ksiGWs_hat, ksiFWs_hat, p_G_hat)
p1rWs_prop2_hat.vec <- p1rWs_gamma(as.numeric(rvec), lamWs_hat,ksiGWs_hat, ksiFWs_hat, p_G_hat)

### Comparison of estimations
########################################################################################################################################## 
cat("p1rW_hat: ",p1rWclass_hat," - p1rWs_hat: ",p1rWsclass_hat," - p1r_nonParam: ",p1r_nonpar,"\n",
    "Theoretical p1rWs is : p1rWs = ",p1rWsclass, "\n"," p_G = ",p_G," - p_Ghat = ",p_G_hat,"\n") # using r = 100
plot(p1r_nonpar.vec,type="l",ylim=c(0,1),xlab="r",ylab="p1r",main = " p1r estim ") # plot of p1r, for r= 10, 20, 30,..., 100
lines(p1rWclass_hat.vec,col="blue")
lines(p1rWsclass_hat.vec,col="green") 
lines(p1rWsclass.vec,col="red") 
legend("topleft",legend=c("p1rNonPar",  "p1rW_hat", "p1rWs_hat","p1rWs"), 
       col=c("black", "blue", "green","red"),lty=c(1,1,1,1),cex=0.8)


# Comparison between p1rstar and prop2_p1rstar,with estimated parameters
########################################################################################################################################## 
cat("pour r =  100 : p1rstar_hat = ",p1rstar_hat," - p1rstar_prop2_hat = ", p1rWs_prop2_hat ,"\n")
# trajectory from r= 10 to r=100
plot(p1rstar_hat.vec,type="l",ylim=c(0,max(p1rstar_hat.vec)),xlab="r",ylab="p1r_star",main = " p1r_starhat & prop2hat") # plot of p1r, for r= 10, 20, 30,..., 100
lines(p1rWs_prop2_hat.vec,col="green")
legend("topleft",legend=c("prop2",  "p1r_star"),
       col=c("green","black"),lty=c(1,1),cex=0.8)

##########################################################################################################################################
##########################################################################################################################################
###################################### 3. Differences between "hibrid" estimations of p1rWs ##############################################

# 3.1 p1rWs using theoretical p_G (not P_G_hat) (but still Zs_hat)
########################################################################################################################################## 
matGZs_ecdf <- GZestimation(X,Zs_hat) 
Wsclass_theta_pg <- WsGMMestim(matGZs_ecdf,p_G,truevalues=NULL) 
ksiFWs_hat_pg <- Wsclass_theta_pg$xiFhat;
ksiGWs_hat_pg <- Wsclass_theta_pg$xiGhat; 
lamWs_hat_pg <- Wsclass_theta_pg$lambdahat
p1rWsclass_andstar_hat_pg <- p1rWs_class(r, lamWs_hat_pg, p_G, ksiGWs_hat_pg, ksiFWs_hat_pg)
p1rWsclass_hat_pg <- p1rWsclass_andstar_hat_pg$p1rWs
p1rstar_hat_pg <- p1rWsclass_andstar_hat_pg$p1rstar
cat("for r=100 the record probability using Ws parametric method: ", p1rWsclass_hat_pg," (unknown parameters but known cdf)","\n") 
p1rWsclass_andstar_hat_pg.vec <- p1rWs_class(rvec, lamWs_hat_pg, p_G, ksiGWs_hat_pg, ksiFWs_hat_pg)
p1rWsclass_hat_pg.vec <- p1rWsclass_andstar_hat_pg.vec$p1rWs  
p1rstar_hat_pg.vec <- p1rWsclass_andstar_hat_pg.vec$p1rstar

# 3.2 p1rWs estimation using  real CDF of G (instead of the ECDF)
########################################################################################################################################## 
CDFmatGZ <- G(Zs_hat, ksi_G, sig_G, mu_G)
Wsclass_theta_cdf <- WsGMMestim(CDFmatGZ,p_G_hat,truevalues=NULL) 
ksiFWs_hat_cdf <- Wsclass_theta_cdf$xiFhat;
ksiGWs_hat_cdf <- Wsclass_theta_cdf$xiGhat; 
lamWs_hat_cdf <- Wsclass_theta_cdf$lambdahat
p1rWsclass_andstar_hat_cdf <- p1rWs_class(r, lamWs_hat_cdf, p_G_hat, ksiGWs_hat_cdf, ksiFWs_hat_cdf)
p1rWsclass_hat_cdf <- p1rWsclass_andstar_hat_cdf$p1rWs
p1rstar_hat_cdf <- p1rWsclass_andstar_hat_cdf$p1rstar
cat("for r=100 the record probability using Ws parametric method: ", p1rWsclass_hat_cdf," (unknown parameters but known cdf)","\n") 
p1rWsclass_andstar_hat_cdf.vec <- p1rWs_class(rvec, lamWs_hat_cdf, p_G_hat, ksiGWs_hat_cdf, ksiFWs_hat_cdf)
p1rWsclass_hat_cdf.vec <- p1rWsclass_andstar_hat_cdf.vec$p1rWs  
p1rstar_hat_cdf.vec <- p1rWsclass_andstar_hat_cdf.vec$p1rstar

# 3.1 p1rWs using real CDF
########################################################################################################################################## 

# Comparison
########################################################################################################################################## 
cat("p1rWs_hat: ",p1rWsclass_hat," - p1r_nonParam: ",p1r_nonpar," - p1rWs_hat_cdf: ",p1rWsclass_hat_cdf," - p1rWs_hat_cdf_pg:",p1rWsclass_hat_pg,"\n",
    "Theoretical p1rWs is : p1rWs = ",p1rWsclass, "\n") # using r = 100
plot(p1rWsclass.vec,type="l",ylim=c(0,1),xlab="r",ylab="p1rWs",main = " p1rWs",col="red") # plot of p1r, for r= 10, 20, 30,..., 100
lines(p1rWclass_hat.vec,col="blue")
lines(p1rWsclass_hat_cdf.vec,col="green") 
lines(p1rWsclass_hat_pg.vec,col="black")
legend("topleft",legend=c("p1rWstheo",  "p1rWs_hat", "p1rWs_hat_cdf","p1rWs_hat_pg"),
       col=c("red", "blue", "green","black"),lty=c(1,1,1,1),cex=0.8)

##########################################################################################################################################
##########################################################################################################################################
##########################################      4. plot p1rWs vs p_G   ###################################################################

# Start parameters 
r=100
pG.vec <- seq(0, 0.9, by=0.01) 
rvec <- seq(10, 100, by=1) 

ksi_G = -0.5; ksi_F = -0.5; sig_G=1; sig_F=1 # On peut changer cette ligne comme on veux
lam = ((sig_F/ksi_F)/(sig_G/ksi_G))^(-1/ksi_G)

J <- 1000; I <- 1000; N <- 1
# mu_F <- 3 # for NonPar

#  4.1 Theoretic pWs  for differents p_G et r (map & plot)
##########################################################################################################################################
pWs_r_PG.mat <- x<-matrix(0, ncol=length(pG.vec), nrow=length(rvec))
for (i in 1:length(pG.vec)){
  for(j in 1:length(rvec)){
    cat("pGF= ",pG.vec[i],", r= ",rvec[j],"\n")
    pWs_r_PG.mat[j,i] <- p1rWs_class(rvec[j], lam, pG.vec[i], ksi_G, ksi_F)$p1rWs
  }
}
image.plot(rvec,pG.vec,pWs_r_PG.mat, main="p1rWs",xlab="r",ylab="pG")
plot(rvec, pWs_r_PG.mat[,50],type="l")
plot(pG.vec, pWs_r_PG.mat[50,],type="l",ylim=c(0,1))

pWs_vs_PG <- rep(0,length(pG.vec))
for (i in 1:length(pG.vec)){
  pWs_vs_PG[i] <- p1rWs_class(r, lam, pG.vec[i], ksi_G, ksi_F)$p1rWs
}
plot(pG.vec,pWs_vs_PG,ylim=c(0,1), col="blue",type="l",main="p1rWs vs pG",ylab="p1rWs",xlab="pG")

# pW method, with p_G=0 assumtion (plot)
##########################################################################################################################################
k = ksi_G/ksi_F
p1rWclass <- p1rfarW(lam,k,r)$p1r 
p1rWclass_plot <- rep(p1rWclass,length(pG.vec))
lines(pG.vec,p1rWclass_plot) # Expected result because we are just changing location parameters and they do not avec p1rWclass

# 4.2 Non-parametric 
##########################################################################################################################################
Nn =100 # We will plot de maen of 100 estimations
mu_F = 10
pGbar.vec <- 1-pG.vec
QuantF <- function(p) mu_F + (sig_F/ksi_F) * ( ( (-log(p))^(-ksi_F) ) -1)
x_G.vec <- QuantF(pGbar.vec)
# p1r_nonpar <- rep(0,length(pG.vec))
mu_G.vec <- x_G.vec + sig_G/ksi_G

p1r_nonpar <- rep(0,length(pG.vec))
for(i in 1:length(x_G.vec)){
  mu_G <- mu_G.vec[i]
  data <- simulWclass(m=I,n=J,N=Nn,ksiX=ksi_G,ksiZ=ksi_F,sigX=sig_G,sigZ=sig_F, 
                       supportauto=FALSE,muX=mu_G, muZ = mu_F, graph=FALSE) 
  Zmat <- data$matZ; Xmat <- data$matX 
  matGZ_ecdf <- GZestimation(Xmat,Zmat)
  p1r_nonpar[i] <- rowMeans(p1r_NonP(matGZ_ecdf,r))
} # Il manque une boucle pour fair ca plusieures fois et prendre la moyenne

# Comparison beetwen methods (plot)
##########################################################################################################################################
plot(pG.vec,p1rWclass_plot,ylim=c(0,1), col="black", main='p1r vs pG', xlab="pG",ylab="p1r",type="l") # NonPar
lines(pG.vec,p1r_nonpar,col="blue") # W
lines(pG.vec,pWs_vs_PG,col="red",) # W
legend("topleft", legend=c("p1rWs", "p1rNonpar","p1rW"),col=c("red", "blue","black"),lty=c(1,1,1),cex=0.8)


# difference betwwen p1rW et p1rNonpar  (maps)
##########################################################################################################################################

Nn = 10
mu_F = 10
pGbar.vec <- 1-pG.vec
QuantF <- function(p) mu_F + (sig_F/ksi_F) * ( ( (-log(p))^(-ksi_F) ) -1)
x_G.vec <- QuantF(pGbar.vec)
# p1r_nonpar <- rep(0,length(pG.vec))
mu_G.vec <- x_G.vec + sig_G/ksi_G

p1r_Nonpar.mat <- matrix(0, ncol=length(pG.vec), nrow=length(rvec))
for(j in 1:length(rvec)){
  r = rvec[j]
  for(i in 1:length(x_G.vec)){
    mu_G <- mu_G.vec[i]
    data <- simulWclass(m=I,n=J,N=Nn,ksiX=ksi_G,ksiZ=ksi_F,sigX=sig_G,sigZ=sig_F, 
                         supportauto=FALSE,muX=mu_G, muZ = mu_F, graph=FALSE) 
    Zmat <- data$matZ; Xmat <- data$matX 
    matGZ_ecdf <- GZestimation(Xmat,Zmat)
    p1r_Nonpar.mat [j,i] <- rowMeans(p1r_NonP(matGZ_ecdf,r))
  }
}
image.plot(rvec,pG.vec,p1r_Nonpar.mat)

difp1r.mat <-pWs_r_PG.mat-p1r_Nonpar.mat

rfcol <- colorRampPalette(brewer.pal(9,'RdBu'))
colsp <- rfcol(64)
r_breaks <- seq(-0.2,0.2,length.out=65)
image.plot(rvec,pG.vec,difp1r.mat,col = colsp,breaks = r_breaks) # p1RNonPar nous donne des valeurs biaises plus petits

# 4.2 Comparison between  pWs & pWs_hat vs p_G
########################################################################################################################################## 

#  p1r methods vs pG (plot)
##########################################################################################################################################
r=100
pG.vec <- seq(0, 0.9, by=0.01) 
rvec <- seq(10, 100, by=1) 

pGbar.vec <- 1-pG.vec
ksi_G = -0.8; ksi_F = -0.8; sig_G=1; sig_F=1 # avec -0.5 l'integration nous donne un erreur
mu_F = 10
QuantF <- function(p) mu_F + (sig_F/ksi_F) * ( ( (-log(p))^(-ksi_F) ) -1)
x_G.vec <- QuantF(pGbar.vec)
mu_G.vec <- x_G.vec + sig_G/ksi_G
lam = ((sig_F/ksi_F)/(sig_G/ksi_G))^(-1/ksi_G)

J <- 1000; I <- 1000; N <- 1

Nn=3
p1rWshat.mat <- matrix(0, ncol=Nn, nrow=length(pG.vec))
p1rstarhat.mat <- matrix(0, ncol=Nn, nrow=length(pG.vec))
for(i in 1:length(x_G.vec)){
  mu_G <- mu_G.vec[i]
  data <- simulWclass(m=I,n=J,N=Nn,ksiX=ksi_G,ksiZ=ksi_F,sigX=sig_G,sigZ=sig_F, 
                       supportauto=FALSE,muX=mu_G, muZ = mu_F, graph=FALSE) 
  Zmat <- data$matZ; Xmat <- data$matX
  
  for (j in 1:Nn){
    Z <- Zmat[,j]
    X <- Xmat[,j]
    x_max <- max(X) 
    Zs_hat <- Z[ Z < x_max]
    Ns_hat <- length(Zs_hat)
    N_G_hat <-J-Ns_hat
    p_G_hat <- (N_G_hat)/J
    matGZs_ecdf <- GZestimation(X,Zs_hat) 
    Wsclass_theta <- WsGMMestim(matGZs_ecdf,p_G_hat,truevalues=NULL) 
    ksiFWs_hat <- Wsclass_theta$xiFhat;
    ksiGWs_hat <- Wsclass_theta$xiGhat; 
    lamWs_hat <- Wsclass_theta$lambdahat
    p1rWsclass_andstar_hat <- p1rWs_class(r, lamWs_hat, p_G_hat, ksiGWs_hat, ksiFWs_hat)
    p1rWshat.mat[i,j] <- p1rWsclass_andstar_hat$p1rWs
    p1rstarhat.mat[i,j] <- p1rWsclass_andstar_hat$p1rstar
  }
  
  matGZ_ecdf <- GZestimation(Xmat,Zmat)
  p1r_nonpar[i] <- rowMeans(p1r_NonP(matGZ_ecdf,r))
} 
p1rWs_hat <- rowMeans(p1rWshat.mat)
p1rstar_hat <- rowMeans(p1rstarhat.mat)

pWs_vs_PG <- rep(0,length(pG.vec))
for (i in 1:length(pG.vec)){
  pWs_vs_PG[i] <- p1rWs_class(r, lam, pG.vec[i], -0.8, -0.8)$p1rWs
}

plot(pG.vec,pWs_vs_PG,ylim=c(0,1),col="red",type="l",main="p1r vs pG",ylab="p1rWs",xlab="pG")
lines(pG.vec,p1r_nonpar,col="blue") 
lines(pG.vec,p1rWs_hat, col="green")
legend("topleft",legend=c("p1rWstheo",  "p1rWs_hat", "p1rNonPar"),
       col=c("red", "green", "blue"),lty=c(1,1,1),cex=0.8)

#  p1r methods vs pG (map)
##########################################################################################################################################

# p1rWs hat
p1rWshat.mat <- matrix(0, ncol=length(pG.vec), nrow=length(rvec))
p1rstarhat.mat <- matrix(0, ncol=length(pG.vec), nrow=length(rvec))
for(i in 1:length(x_G.vec)){
  mu_G <- mu_G.vec[i]
  data <- simulWclass(m=I,n=J,N=Nn,ksiX=ksi_G,ksiZ=ksi_F,sigX=sig_G,sigZ=sig_F, 
                      supportauto=FALSE,muX=mu_G, muZ = mu_F, graph=FALSE) 
  Zmat <- data$matZ; Xmat <- data$matX
  
  p1rWshat_pG.mat <- matrix(0, ncol=Nn, nrow=length(pG.vec))
  p1rstarhat_pG.mat <- matrix(0, ncol=Nn, nrow=length(pG.vec))
  for (j in 1:Nn){
    Z <- Zmat[,j]
    X <- Xmat[,j]
    x_max <- max(X) 
    Zs_hat <- Z[ Z < x_max]
    Ns_hat <- length(Zs_hat)
    N_G_hat <-J-Ns_hat
    p_G_hat <- (N_G_hat)/J
    matGZs_ecdf <- GZestimation(X,Zs_hat) 
    Wsclass_theta <- WsGMMestim(matGZs_ecdf,p_G_hat,truevalues=NULL) 
    ksiFWs_hat <- Wsclass_theta$xiFhat;
    ksiGWs_hat <- Wsclass_theta$xiGhat; 
    lamWs_hat <- Wsclass_theta$lambdahat
    p1rWsclass_andstar_hat <- p1rWs_class(rvec, lamWs_hat, p_G_hat, ksiGWs_hat, ksiFWs_hat)
    p1rWshat_pG.mat[,j] <- p1rWsclass_andstar_hat$p1rWs
    p1rstarhat_pG.mat[,j] <- p1rWsclass_andstar_hat$p1rstar
  }
  p1rWshat.mat[,i] <- rowMeans(p1rWshat_pG.mat)
  p1rstarhat.mat[,i]<- rowMeans(p1rstarhat_pG.mat)
  } 
image.plot(rvec,pG.vec,p1rWshat.mat, main="p1rWs_hat",xlab="r",ylab="pG")

#p1rWs theo
pWs_r_PG.mat <- x<-matrix(0, ncol=length(pG.vec), nrow=length(rvec))
for (i in 1:length(pG.vec)){
  for(j in 1:length(rvec)){
    cat("pGF= ",pG.vec[i],", r= ",rvec[j],"\n")
    pWs_r_PG.mat[j,i] <- p1rWs_class(rvec[j], lam, pG.vec[i], ksi_G, ksi_F)$p1rWs
  }
}
image.plot(rvec,pG.vec,pWs_r_PG.mat, main="p1rWs",xlab="r",ylab="pG")

#diff
difp1r.mat <-pWs_r_PG.mat-p1rWshat.mat
rfcol <- colorRampPalette(brewer.pal(9,'RdBu'))
colsp <- rfcol(64)
r_breaks <- seq(-0.2,0.2,length.out=65)
image.plot(rvec,pG.vec,difp1r.mat,col = colsp,breaks = r_breaks,
           main="pWs - pWshat.",xlab="r",ylab="pG") # le diff  est des deux cotés








