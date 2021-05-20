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

# we will need the functions from the file Fonctions-for-far-WeibullEstimation.R

#
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

functiong <- function(theta,vecx){
  lambdaval=theta[1] ; kval = theta[2]
  p12val <- laplaceWeibull(j=1,lambda=lambdaval,k=kval,lowerbnd=10^(-8))
  p13val <- laplaceWeibull(j=2,lambda=lambdaval,k=kval,lowerbnd=10^(-8))
  M <- cbind( vecx - p12val ,  vecx^2 - p13val )
  return(M)
}

#
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
  return(matGm)  # dans notre cas elle nous retoune G(Z)m empirique pour le couple (X,Z)
}

#
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