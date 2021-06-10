###################### Exemples ###############################
library(evd)
# RNG
#set.seed(1)
set.seed(Sys.time())

# Easy visualisation distributions
size <- 20
tt <- seq.int(size)
mu = seq(0, 5, length.out = size) 
X_gev = rgev(size, loc = 0, scale = 1, shape = 0)
Z_gev = rgev(size, loc = 2, scale = 1, shape = 0)
Z_gev_trend = rgev(size, loc = mu, scale = 1, shape = 0)

X_norm = rnorm(size,sd=2)
Z_norm = rnorm(size,0.4,2)

############# Example 1
# Gumbel - (set.seed(292))

# Simulations
size <-  250*4
rp <- 50
tt <- seq.int(size)/size
mu = 2 + seq(0, 5, length.out = size) 
X = rgev(size * 5, loc = 0, scale = 1, shape = 0)
Z = rgev(size, loc = mu, scale = 1, shape = 0)
theta_theo <- 1 / exp(mu)
p12_theo <- 1 / (1 + theta_theo)

############ Example 2
# Frechet - set.seed(304)

# Simulations
size <-  250 * 4
rp <- 50
tt <- seq.int(size)/size
sigma <- seq(1, 2, length.out = size)
xi = 0.5
x = rgev(size * 1/4, loc = 1, scale = xi, shape = xi)
z = rgev(size, loc = sigma, scale = xi * sigma, shape = xi)
theta_theo <- sigma^(-1 / xi)
p12_theo <- 1 / (1 + theta_theo)
G_theo <- function(z) pgev(z, loc = 1, scale = xi, shape = xi)
far_theo_rp <- (1 - theta_theo) * (1 - 1/rp)

#########################  Kernels exemples ##########################################
# Epanechnikov density distribution
dEpan <- function(x){
  k <- (3/4)* (1-x^2)
  k[-1>x] <- 0
  k[x>1] <- 0
  return (k)
} 
# Normal distribution
Kern_gauss <- function(x)exp(-x^2/2)

###################### p12 Non parametric computation ################
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

################################### Computation of IC ################################
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

IC <- function(X,Z,tt,t_eval,h,kern=dEpan) {
  # computation of G_Z
  G_emp <- ecdf(X)
  G_Z <- G_emp(Z)
  # integral of squared kernel
  K22<-integrate(dEpan_2,-1,1)$value
  # density distribution of t 
  f_t<-kern_dens(tt,tt,h)
  # variance of t 
  sigma_t <- var(tt)
  # Var(A):
  VarA <- (sigma_t * K22 * f_t )^0.5
  # Var(B):
  Khj <- outer(t_eval, tt, function(zz,z) dEpan((zz - z) / h))
  Khi <- outer(t_eval, tt, function(zz,z) dEpan((zz - z) / h))
  E <- mean(outer(G_Z,G_Z, pmin) - outer(G_Z,G_Z,"*"))
  VarB_num <- as.numeric(rowSums(Khj) %*% rowSums(Khi))
  VarB_denom <- length(tt)*(rowSums(Khj))^2
  VarB <- VarB_num*E/VarB_denom
  # sigma = Var(A)+ Var(B)
  sigma_m <- VarA + VarB
  #Computation of p12_hat
  p12_hat<-p12_NonPar(X,Z,tt,t_eval,h,kern)
  # Computation of CI
  s <- (sigma_m)^1/2
  error <- qnorm(0.975)*s/sqrt(length(tt))
  left <- p12_hat-error
  right <- p12_hat+error
  return(list(low=left,high=right))
}

##### Computation and plot of p12_hat and CI for 1 couple (X,Z)
library(ggplot2)
p12_hat<-p12_NonPar(X,Z,tt,tt,0.11)
ic<-IC(X,Z,tt,tt,0.11) 
df <-data.frame(x=tt,y=p12_hat,z1=ic$low, z2=ic$high,heo=p12_theo)
p12plot<-ggplot(df,aes(x=tt,y=p12_hat)) + geom_line(colour="red") 
p12plot<- p12plot + geom_line(aes(y=p12_theo), colour="black")
p12plot <- p12plot + geom_ribbon(aes(ymin=ic$low, ymax=ic$high), linetype=2, alpha=0.1) 
p12plot <- p12plot + ggtitle("p12 evolution over time") + ylab("p12") + xlab("time")
p12plot

################################# Quality of convergence #################################

# Computation of multiple samples (X,Z)
N <- 10
sigma <- seq(1, 2, length.out = size)
xi = 0.5
samplex <- NULL; for (i in 1:N){samplex <- cbind(samplex,rgev(size * 1/4, loc = 1, scale = xi, shape = xi))}
samplez <- NULL; for (i in 1:N){samplez <- cbind(samplez,rgev(size, loc = sigma, scale = xi * sigma, shape = xi))}

# Computation of p12_hat_sample for each (Xi,Zi)
p12_hat_samples <- matrix(0,size,N); ic_samples_high <- matrix(0,size,N); ic_samples_low <- matrix(0,size,N)
for (i in 1:N){
  p12_hat_samples[,i] <- p12_hat_samples[,i] + p12_NonPar(samplex[,i],samplez[,i],tt,tt,0.11)
  ic_samples_high[,i] <- ic_samples_high[,i] + as.vector(IC(samplex[,i],samplez[,i],tt,tt,0.11)$high)
  ic_samples_low[,i] <- ic_samples_low[,i] + as.vector(IC(samplex[,i],samplez[,i],tt,tt,0.11)$low)
}

# Computation on average for big N
p12_hat_moyen <- rowMeans(as.matrix(p12_hat_samples)); ic_samples_high <- rowMeans(as.matrix(ic_samples_high)); ic_samples_low <-rowMeans(as.matrix(ic_samples_low))

# Plot
library(ggplot2)
df <-data.frame(x=tt,y=p12_hat_moyen,z1=ic_samples_low , z2=ic_samples_high ,theo=p12_theo)
p12plot<-ggplot(df,aes(x=tt,y=p12_hat_moyen)) + geom_line(colour="red") 
p12plot<- p12plot + geom_line(aes(y=p12_theo), colour="black")
p12plot <- p12plot + geom_ribbon(aes(ymin=ic_samples_low, ymax=ic_samples_high), linetype=2, alpha=0.1) 
p12plot <- p12plot + ggtitle("p12 evolution over time") + ylab("p12") + xlab("time")
p12plot

# Error between p12_theo and p12_hat
p12_error <- abs(p12_theo - p12_hat_moyen)
plot(tt, p12_error, main="Error of estimation", ylab="error", xlab="temps")

# Incertitude of estimation
ic_error <- ic_samples_high - ic_samples_low
plot(tt, ic_error, main="Incertitude of estimation", ylab=" difference", xlab="temps")

# max error evolution within N (Global error)
Nn=50
max_error_evolution <- rep(0,Nn)
for (j in 1:Nn) {
  samplex <- NULL; for (i in 1:j){samplex <- cbind(samplex,rgev(size * 1/4, loc = 1, scale = xi, shape = xi))}
  samplez <- NULL; for (i in 1:j){samplez <- cbind(samplez,rgev(size, loc = sigma, scale = xi * sigma, shape = xi))}
  p12_hat_samples <- matrix(0,size,j)
  for (i in 1:j){
    p12_hat_samples[,i] <- p12_hat_samples[,i] + p12_NonPar(samplex[,i],samplez[,i],tt,tt,0.11)
  }
  p12_hat_moyen <- rowMeans(p12_hat_samples)
  p12_error <- abs(p12_theo - p12_hat_moyen[j])
  max_error_evolution[j] <- max(p12_error)
}
plot(1:Nn,max_error_evolution)

# Conditional error: 

Nn=50
matrix_p12_moyen <- NULL
mean_p12 <- matrix(0,size,Nn)
conditional_error <- matrix(0,size,Nn)
for (j in 1:Nn) {
  samplex <- NULL; for (i in 1:j){samplex <- cbind(samplex,rgev(size * 1/4, loc = 1, scale = xi, shape = xi))}
  samplez <- NULL; for (i in 1:j){samplez <- cbind(samplez,rgev(size, loc = sigma, scale = xi * sigma, shape = xi))}
  p12_hat_samples <- NULL; for (i in 1:j){p12_hat_samples<-cbind(p12_hat_samples,p12_NonPar(samplex[,i],samplez[,i],tt,tt,0.11))}
  mean_p12[,j] <- rowMeans(p12_hat_samples)
  conditional_error[,j] <- abs(mean_p12[,j] - p12_theo)
}

matplot(x= tt , y= as.matrix(conditional_error), type='l', pch=1, 
        col= 2:5, xlab='tt', ylab = 'error')
###########
Nn=c(5,50,500)
mean_p12 <- matrix(0,size,length(Nn))
conditional_error <- matrix(0,size,length(Nn))
for (j in 1:length(Nn)) {
  samplex <- NULL; for (i in 1:Nn[j]){samplex <- cbind(samplex,rgev(size * 1/4, loc = 1, scale = xi, shape = xi))}
  samplez <- NULL; for (i in 1:Nn[j]){samplez <- cbind(samplez,rgev(size, loc = sigma, scale = xi * sigma, shape = xi))}
  p12_hat_samples <- NULL; for (i in 1:Nn[j]){p12_hat_samples<-cbind(p12_hat_samples,p12_NonPar(samplex[,i],samplez[,i],tt,tt,0.11))}
  mean_p12[,j] <- rowMeans(p12_hat_samples)
  conditional_error[,j] <- abs(mean_p12[,j] - p12_theo)
}

matplot(x= tt , y= conditional_error, type='l', pch=1, 
        col= 2:5, xlab='tt', ylab = 'error')

for (j in Nn) {print(j)}
# rouge,vert, bleu 

################## Distributions with random location parameter ###############################
# distributions
ui<-rnorm(2);bi<-rnorm(2,mean=1,sd=0.1)
X1<-rgev(size, loc=ui[1], scale=1, shape=0); X2<-rgev(size, loc=ui[2], scale=1, shape=0)
Z1<-rgev(size, loc = ui[1]+bi[1]*tt, scale = 1, shape = 0); Z2<-rgev(size, loc = ui[2]+bi[2]*t, scale = 1, shape = 0)
# p12 computation
p12_1<-p12_NonPar(X1,Z1,tt,tt,0.2); p12_2<-p12_NonPar(X2,Z2,tt,tt,0.2);ic_1<-IC(X1,Z1,tt,tt,0.2) ;ic_2<-IC(X2,Z2,tt,tt,0.2) 
#plot
df_Ex <-data.frame(temp=tt,x=p12_1,y=p12_2,z1=ic_1$low,z2=ic_1$high,z3=ic_2$low,z4=ic_2$high)
p12plot<-ggplot(df_Ex,aes(x=tt,y=p12_1))+ geom_line(colour="red")+geom_ribbon(aes(ymin=ic_1$low, ymax=ic_1$high), linetype=2, alpha=0.1)+
  geom_line(aes(y=p12_2), colour="black")+geom_ribbon(aes(ymin=ic_2$low, ymax=ic_2$high), linetype=2, alpha=0.1)
p12plot

############################### Bandwith selection #################################

#####  Some plots of p12_hat and G(Z) computation with different bandwiths

h1=0.01; h2=0.1; h3=1

# Gumbel
p12_g1<-p12_NonPar(X,Z,tt,tt,h1); p12_g2<-p12_NonPar(X,Z,tt,tt,h2); p12_g3<-p12_NonPar(X,Z,tt,tt,h3)
G_emp_X <- ecdf(X); G_Z <- G_emp(Z)
df_g<-data.frame(x=tt, y1=p12_g1, y2=p12_g2, y3=p12_g3, theo=p12_theo, g=G_Z)
hchoix_g<-ggplot(df_g,aes(x=tt,y=G_Z)) + geom_point(colour="#999999",size=0.5, shape=23)+
  geom_line(aes(y=p12_theo), colour="black")+
  geom_line(aes(x=tt,y=p12_g1,colour="red"))+geom_line(aes(x=tt,y=p12_g2,colour="blue"))+
  geom_line(aes(x=tt,y=p12_g3,colour="green"))
hchoix_g

# Frechet
p12_f1<-p12_NonPar(x,z,tt,tt,h1); p12_f2<-p12_NonPar(x,z,tt,tt,h2); p12_f3<-p12_NonPar(x,z,tt,tt,h3)
G_emp_x<-ecdf(x); G_z<-G_emp(z)
df_f<-data.frame(x=tt,y1=p12_f1,y2=p12_f2,y3=p12_f3,theo=p12_theo,g=G_z)
hchoix_f<-ggplot(df_h_2,aes(x=tt,y=G_z)) + geom_point(colour="#999999",size=0.5, shape=23)+
  geom_line(aes(y=p12_theo), colour="black")+
  geom_line(aes(x=tt,y=p12_f1,colour="red"))+geom_line(aes(x=t,y=p12_f2,colour="blue"))+
  geom_line(aes(x=tt,y=p12_f3,colour="green"))
hchoix_f

##### Optimal bandwidth computation

# Computation of vector of errors

h_seq<-seq(from=0.01, to=1,by=0.02)

size <-  250 * 4
tt <- seq.int(size)/size
t_eval <- seq.int(size)/size

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
error_vector<-CV_error(X,Z,tt,t_eval,h_seq,kern= dEpan)

# Error plot 
plot(x=h_seq, y=error_vector, type="b", lwd=3, col="blue",
     xlab="Smoothing bandwidth", ylab="LOOCV prediction error")

h_opt <- h_seq[which(CV_err_h == min(CV_err_h))]; h_opt # optimal bandwith = 0.11
opt_err <-min(CV_err_h); opt_err # optimal erreur = 0.03589835



############################# Add "real" variable Y #################################


###################### X stationary and Z with trend
size <-  250*4
tt <- seq.int(size)/size
# numer of models 
m <- 2
# paraters of the models
ui<-rnorm(m) 
bi<-rnorm(m,mean=1,sd=0.1)

scale_Z_real=1
shape_Z_real=0

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

# "real" factual(Z) and couterfactual(X) run
Z_real<-Z_real_comp(size,scale=1,shape=0)
X_real<-rgev(size*5, loc=mean, scale=1, shape=0)
# run from model 1
X1<-rgev(size*5, loc=ui[1], scale=1, shape=0)
Z1<-rgev(size, loc = ui[1]+bi[1]*tt, scale = 1, shape = 0)
# run from model 2
X2<-rgev(size*5, loc=ui[2], scale=1, shape=0)
Z2<-rgev(size, loc = ui[2]+bi[2]*tt, scale = 1, shape = 0)
# "real" p12 , p12_hat and CI of each model 
p12_real<-p12_NonPar(X_real,Z_real,tt,tt,0.11); p12_1<-p12_NonPar(X1,Z1,tt,tt,0.11); p12_2<-p12_NonPar(X2,Z2,tt,tt,0.11)
ic_1<-IC(X1,Z1,tt,tt,0.11) ;ic_2<-IC(X2,Z2,tt,tt,0.11) 
# plot
df_Ex <-data.frame(temp=tt, x=p12_real, y1=p12_1, y2=p12_2,z1=ic_1$low,z2=ic_1$high,z3=ic_2$low,z4=ic_2$high)
p12plot<-ggplot(df_Ex,aes(x=tt,y=p12_1))+ geom_line(colour="red")+geom_ribbon(aes(ymin=ic_1$low, ymax=ic_1$high), linetype=2, alpha=0.1)+
  geom_line(aes(y=p12_2), colour="black")+geom_ribbon(aes(ymin=ic_2$low, ymax=ic_2$high), linetype=2, alpha=0.1) +  geom_line(y=p12_real, colour="blue")
p12plot

dev.off()

################# Monde avant "effet" du forçage anthropique 

# run from model 1
X1<-rgev(size*5, loc=ui[1], scale=1, shape=0); Z1<-rgev(size, loc = ui[1], scale = 1, shape = 0)
# run from model 2
X2<-rgev(size*5, loc=ui[2], scale=1, shape=0); Z2<-rgev(size, loc = ui[2], scale = 1, shape = 0)

p12_1<-p12_NonPar(X1,Z1,tt,tt,0.11); p12_2<-p12_NonPar(X2,Z2,tt,tt,0.11)
ic_1<-IC(X1,Z1,tt,tt,0.11) ;ic_2<-IC(X2,Z2,tt,tt,0.11) 

p12hat_moyen<- (p12_1 + p12_2)/2

df_Ex <-data.frame(temp=tt, y1=p12_1, y2=p12_2, y3= p12hat_moyen, z1=ic_1$low,z2=ic_1$high,z3=ic_2$low,z4=ic_2$high)
p12plot<-ggplot(df_Ex,aes(x=tt,y=p12_1))+ geom_line(colour="red")+geom_ribbon(aes(ymin=ic_1$low, ymax=ic_1$high), linetype=2, alpha=0.1)+
  geom_line(aes(y=p12_2), colour="blue")+geom_ribbon(aes(ymin=ic_2$low, ymax=ic_2$high), linetype=2, alpha=0.1) +  geom_line(y=p12hat_moyen, colour="green")+
  geom_hline(yintercept = 0.5) + ggtitle("p02 estimation over time") + ylab("p02") + xlab("time")
p12plot

################################## p13_t computation #######################################

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

p13_hat<-p13_NonPar(X,Z,tt,tt,0.11)
plot(tt, p13_hat)
p12_hat<-p12_NonPar(X,Z,tt,tt,0.11)
plot(tt, p12_hat)

dev.off()

############# Using non-parametric kernel computed p12, p13 for p1r parametric computation ##################

library(stats4); library(gmm); library(stats); library(np); library(EWGoF)

# simulate multiple stationnary {(X)t,(Z)t} W-class trajectories

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
              "lam"=((sigX/ksiX)/(sigZ/ksiZ))^(1/ksiX),
              "k"=ksiX/ksiZ) 
  )  
}

# Estimate p12t and p13t for each trajectory {(X)t,(Z)t}

GZestimation <- function(matX,matZ,tt,t_eval,h,kern=dEpan){
  
  matX <- as.matrix(matX)
  matZ <- as.matrix(matZ) 
  dimnsZ=dim(matZ)
  matp12 <- matrix(rep(0,prod(dimnsZ)),nrow=dimnsZ[1],ncol=dimnsZ[2])
  matp13 <- matrix(rep(0,prod(dimnsZ)),nrow=dimnsZ[1],ncol=dimnsZ[2])
  for (j in 1:dimnsZ[2]){
    X <- matX[,j]
    Z <- matZ[,j]
    matp12[,j] <- p12_NonPar(X,Z,tt,t_eval,h,kern= dEpan)
    matp13[,j] <- p13_NonPar(X,Z,tt,t_eval,h,kern= dEpan)
  }
  return(list("matp12"=matp12,"matp13"=matp13))
}
  

# funcLaplace and laplace Weibull are necessary to compute E(G(Z)^j)

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

functiong <- function(theta,vecx){
  lambdaval=theta[1] ; kval = theta[2]
  p12 = vecx[1] ; p13 = vecx[2]
  p12val <- laplaceWeibull(j=1,lambda=lambdaval,k=kval,lowerbnd=10^(-8))
  p13val <- laplaceWeibull(j=2,lambda=lambdaval,k=kval,lowerbnd=10^(-8))
  M <- c( p12 - p12val ,  p13 - p13val )
  return(M)
}

################### Estimation of lambda_t and k_t using different méthods

# 1.
# Using fsolve() with jacobian matrix. we will call this method "weibullGMMestim"

J12 <- function(x){jacobianFunctiong12(x[1],x[2])}

weibullGMMestim <- function(matp12,matp13,truevalues=NULL){
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
      
      fg <- function(theta){functiong(theta,vecx=phat)}
      
      EGMMweibull <- fsolve(fg,x0=startval_t,J=J12
      )
      
      lambdahat.vec[i] <- EGMMweibull$x[1]
      khat.vec[i] <- EGMMweibull$x[2] 
    }
    lambdahat.mat[,j] <- lambdahat.vec
    khat.mat[,j] <- khat.vec
  }   
  return( list("lambdahat"=lambdahat.mat,"khat"=khat.mat) ) 
}

# 2.
# Using fsolve3() without jacobian matrix. we will call this method "weibullGMMestim3"
# fsolve3() is a modification of fsolev(), where we impose [lambda,k]_xnew <- ([lambda,k]_x0)/2 when lamba<0

fsolve3<-function (f, x0, J = NULL, maxiter = 100, tol = .Machine$double.eps^(0.5), 
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
    sol = broyden3(f, x0, J0 = J(x0), maxiter = maxiter, 
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

broyden3<-function (Ffun, x0,J0 = NULL, ..., maxiter = 100, tol = .Machine$double.eps^(1/2)) 
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
  if (xnew[1]<0){xnew <- x0/2}
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

weibullGMMestim3 <- function(matp12,matp13,truevalues=NULL){
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
      
      EGMMweibull <- fsolve3(functiong,x0=startval_t,vecx=phat
                            )
      
      lambdahat.vec[i] <- EGMMweibull$x[1]
      khat.vec[i] <- EGMMweibull$x[2] 
    }
    lambdahat.mat[,j] <- lambdahat.vec
    khat.mat[,j] <- khat.vec
  }   
  return( list("lambdahat"=lambdahat.mat,"khat"=khat.mat) ) 
}

library(pracma)

# errors: we observe sometimes two type of errors. 
# error 1: when n is very big (1000) , k gets negative -> non-finite function value
# error 2: when n is little (20 or 30) -> Error in norm(s, "F") : 'A' must be a numeric matrix 


####### Optisation approcah, using quasi-Newton methods (BFGS and L-BFGS-B)

fg1<- function (theta,vecx) crossprod(functiong(theta,vecx)) #elle marche
# fg2 <- function (theta){fg1(theta,vecx=phat)} # necessary when adding 

# Gradfg <- t(JX)%*%(functiong(theta,phat)[[1]]) + crossprod(fg,JX) not ready

# choice of method
EOptimweibull <- optim(c(0.7, 1.3),fg1,method="BFGS",phat=c(0.86, 0.79))

EOptimweibull <- optim(c(0.7, 1.3),fg1,method="L-BFGS-B",phat=c(0.86, 0.79),lower = 0.001)

# weibullGMMestim4 estimates lambda_t, k_t using the optimisation method of choice

weibullGMMestim4 <- function(matp12,matp13,truevalues=NULL){
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
      EGMMweibull <- optim(par=startval_t,fn=fg1,method="L-BFGS-B",vecx=phat,lower = 0.001
      )
      #EGMMweibull <- optim(par=startval_t,fn=fg1,method="BFGS",vecx=phat
      #)
      
      lambdahat.vec[i] <- EGMMweibull$par[1]
      khat.vec[i] <- EGMMweibull$par[2] 
    }
    lambdahat.mat[,j] <- lambdahat.vec
    khat.mat[,j] <- khat.vec
  }   
  return( list("lambdahat"=lambdahat.mat,"khat"=khat.mat) ) 
}

####### Computation of p1r_t, far_t using {lambda_t,k_t}

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

#################### Observation of distributions in time #########################

# test ( just for Gumbel, X stationary, Z non-sationary with linear growth in location parameter)
plotd_time <- function(size,locat,mu,scal,shpe){
  X = rgev(size, loc = locat, scale = scal, shape = shpe)
  mu_t= c(mu[1],mu[length(mu)/4],mu[length(mu)/2], mu[3*length(mu)/4], mu[length(mu)])
  plot(density(X),col="black",)
  for (i in 1:length(mu_t)){
    Z = rgev(size, loc = mu_t[i], scale = 1, shape = 0)
    lines(density(Z),col="gray")
  }
}

n=1000
mu = seq(0, 5, length.out = n) 
plotd_time(n,0,mu,1,0)

mu = seq(0, 3, length.out = n) 
plotd_time(n,0,mu,1,0.2)


########################### Simulate W-class trajectories with non-sationnary Z #################

simulWclass_trend <- function(m,n,N,ksiX,ksiZ,sigX,tt.vec,b.vec,supportauto=TRUE,muX=0,muZ=0,sigZ=0,graph=TRUE){

  if (sigZ == 0) { sigZ = sigX } 
  muZ <- matrix(nrow=n,ncol=N)
  if (supportauto == TRUE){ 
    for (i in 1:n){
      muZ[i] <- muX + sigZ/ksiZ - sigX/ksiX - b.vec[i]*tt.vec[i]}
    }
  
  matX <- matrix(nrow=m,ncol=N) 
  matZ <- matrix(nrow=n,ncol=N) 
  
  for (j in 1:N){
    matX[,j] = muX + (sigX/ksiX)*( (-log(runif(m)))^(-ksiX) - 1 )  
    matZ[,j] = (muZ+b.vec*tt.vec) + (sigZ/ksiZ)*( (-log(runif(n)))^(-ksiZ) - 1 )  
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
              "lam"=((sigX/ksiX)/(sigZ/ksiZ))^(1/ksiX),
              "k"=ksiX/ksiZ) 
  )  
}

b=rep(1,1000)
tt <- seq.int(1000)/1000
data<-simulWclass_trend(m=1000,n=1000,N=1,ksiX=-0.20,ksiZ=-0.25,sigX=1,tt.vec=tt,b.vec=b,sigZ=1)


# W-class distributions in time

plotd_time_Wclass <- function(m,n,N,ksiX,ksiZ,sigX,tt.vec,b.vec,muX=0,muZ=0,sigZ=0){
  
  if (sigZ == 0) { sigZ = sigX } 
  muZ <- matrix(nrow=n,ncol=N)
   
  for (i in 1:n){
      muZ[i] <- muX + sigZ/ksiZ - sigX/ksiX - b.vec[i]*tt.vec[i]}
  
  matX <- matrix(nrow=m,ncol=N) 
  matZ <- matrix(nrow=n,ncol=N) 
  
  #for (j in 1:N){
   # matX[,j] = muX + (sigX/ksiX)*( (-log(runif(m)))^(-ksiX) - 1 )  
    #matZ[,j] = (muZ+b.vec*tt.vec) + (sigZ/ksiZ)*( (-log(runif(n)))^(-ksiZ) - 1 )  
    #}
  X = rgev(m, loc = muX, scale = sigX, shape = ksiX)
  muZ_t= c(muZ[1],muZ[length(muZ)/4],muZ[length(muZ)/2], muZ[3*length(muZ)/4], muZ[length(muZ)])
  tt.vec_t= c(tt.vec[1],tt.vec[length(muZ)/4],tt.vec[length(muZ)/2], tt.vec[3*length(muZ)/4], tt.vec[length(muZ)])
  b.vec_t= c(b.vec[1],b.vec[length(muZ)/4],b.vec[length(muZ)/2], b.vec[3*length(muZ)/4], b.vec[length(muZ)])
  
  plot(density(X),col="black",)
  for (i in 1:length(muZ_t)){
    Z = rgev(n, loc = muZ_t[i]+b.vec_t[i]*tt.vec_t[i], scale = sigZ, shape = ksiZ)
    lines(density(Z),col="gray")
  }
  
}

plotd_time_Wclass(m=1000,n=1000,N=1,ksiX=-0.05,ksiZ=-0.05,sigX=1,tt.vec=tt,b.vec=b,sigZ=1)
plotd_time_Wclass(m=1000,n=1000,N=1,ksiX=-0.20,ksiZ=-0.25,sigX=1,tt.vec=tt,b.vec=b,sigZ=1)
  
##############  Statistiques tests for tre transfomation fonction f(W(lambda,k))->W(1,1) ######

# N trajectories of Weibull sample of size n and parameters lambda,k
k=0.8
lambda= 0.2
rweibull(20,k,lambda)



