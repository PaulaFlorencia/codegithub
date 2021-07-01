################################################################################################
################################################################################################
### This file contains test executions of the functions in  p.R                              ###
################################################################################################
################################################################################################


################################################################################################
#### 1. p12_t observations                                                                  ####
################################################################################################


##############################   Test trajectories   ###########################################
#                          ( with known theoretic p1r)                                         #

library(evd)
library(ggplot2)

set.seed(Sys.time())


#### Gumbel trajectories  - (set.seed(292))

size <-  250*4
rp <- 50
tt <- seq.int(size)/size
mu = 2 + seq(0, 5, length.out = size) 
X = rgev(size * 5, loc = 0, scale = 1, shape = 0)
Z = rgev(size, loc = mu, scale = 1, shape = 0)
theta_theo <- 1 / exp(mu)
p12_theo <- 1 / (1 + theta_theo)

#### Frechet trajectories - set.seed(304)

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


##############################   p12_t and CI   ################################################
#                       (for just 1 trajectory , N = 1)                                        #

p12_hat<-p12_NonPar(X,Z,tt,tt,0.11)
ic<-IC(X,Z,tt,tt,0.11) 
df <-data.frame(x=tt,y=p12_hat,z1=ic$low, z2=ic$high,heo=p12_theo)
p12plot<-ggplot(df,aes(x=tt,y=p12_hat)) + geom_line(colour="red")  +
  geom_line(aes(y=p12_theo), colour="black")
p12plot <- p12plot+geom_ribbon(aes(ymin=ic$low, ymax=ic$high), linetype=2, alpha=0.1) 
p12plot <- p12plot +ggtitle("p12 evolution over time") + ylab("p12") + xlab("time")
p12plot


##############################  Quality of convergence   ########################################

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
p12_hat_moyen <- rowMeans(as.matrix(p12_hat_samples))
ic_samples_high <- rowMeans(as.matrix(ic_samples_high)); ic_samples_low <-rowMeans(as.matrix(ic_samples_low))

# Plot
library(ggplot2)
df <-data.frame(x=tt,y=p12_hat_moyen,z1=ic_samples_low , z2=ic_samples_high ,theo=p12_theo)
p12plot<-ggplot(df,aes(x=tt,y=p12_hat_moyen)) + geom_line(colour="red") + 
  geom_line(aes(y=p12_theo), colour="black")+ 
  geom_ribbon(aes(ymin=ic_samples_low, ymax=ic_samples_high), linetype=2, alpha=0.1) + 
  ggtitle("p12 evolution over time") + ylab("p12") + xlab("time")
p12plot

# Error between p12_theo and p12_hat
p12_error <- abs(p12_theo - p12_hat_moyen)
plot(tt, p12_error, main="Error of estimation", type="l",ylab="error", xlab="temps")

# Incertitude of estimation
ic_error <- ic_samples_high - ic_samples_low
plot(tt, ic_error, main="Incertitude of estimation",type="l", ylab=" difference", xlab="temps")

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

# Conditional error
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
matplot(x= tt , y= as.matrix(conditional_error), type='l', pch=1, col= 2:5, xlab='tt', ylab = 'error')

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

matplot(x= tt , y= conditional_error, type='l', pch=1, col= 2:5, xlab='tt', ylab = 'error')

for (j in Nn) {print(j)}
# rouge,vert, bleu 

################################# Bandwith selection ############################################

#####  Some plots of p12_hat and G(Z) computation with different bandwiths
h1=0.01; h2=0.1; h3=1

# Gumbel
p12_g1<-p12_NonPar(X,Z,tt,tt,h1); p12_g2<-p12_NonPar(X,Z,tt,tt,h2); p12_g3<-p12_NonPar(X,Z,tt,tt,h3)
G_emp_X <- ecdf(X); G_Z <- G_emp(Z)
df_g<-data.frame(x=tt, y1=p12_g1, y2=p12_g2, y3=p12_g3, theo=p12_theo, g=G_Z)
hchoix_g<-ggplot(df_g,aes(x=tt,y=G_Z)) + geom_point(colour="#999999",size=0.5, shape=23)+
  geom_line(aes(y=p12_theo), colour="black")+ geom_line(aes(x=tt,y=p12_g1,colour="red"))+
  geom_line(aes(x=tt,y=p12_g2,colour="blue"))+ geom_line(aes(x=tt,y=p12_g3,colour="green"))
hchoix_g

# Frechet
p12_f1<-p12_NonPar(x,z,tt,tt,h1); p12_f2<-p12_NonPar(x,z,tt,tt,h2); p12_f3<-p12_NonPar(x,z,tt,tt,h3)
G_emp_x<-ecdf(x); G_z<-G_emp(z)
df_f<-data.frame(x=tt,y1=p12_f1,y2=p12_f2,y3=p12_f3,theo=p12_theo,g=G_z)
hchoix_f<-ggplot(df_h_2,aes(x=tt,y=G_z)) + geom_point(colour="#999999",size=0.5, shape=23)+
  geom_line(aes(y=p12_theo), colour="black")+geom_line(aes(x=tt,y=p12_f1,colour="red"))+
  geom_line(aes(x=t,y=p12_f2,colour="blue"))+geom_line(aes(x=tt,y=p12_f3,colour="green"))
hchoix_f

##### Optimal bandwidth computation
h_seq<-seq(from=0.01, to=1,by=0.02)
# h_seq <-seq(from=0.1,to=1, by=0.1)
size <-  250 * 4
tt <- seq.int(size)/size
t_eval <- seq.int(size)/size

error_vector<-CV_error(x,z,tt,tt,h_seq,kern= dEpan)

# Error plot 
plot(x=h_seq, y=error_vector, type="b", lwd=3, col="blue",
     xlab="Smoothing bandwidth", ylab="LOOCV prediction error")
h_opt <- h_seq[which(CV_err_h == min(CV_err_h))]; h_opt # optimal bandwith = 0.11
opt_err <-min(CV_err_h); opt_err # optimal erreur = 0.03589835


################## Trajectories with random location parameter ###############################

#### X stationary and Z non stationary with linear trend un location parameter

ui<-rnorm(2); bi<-rnorm(2,mean=1,sd=0.1)
X1<-rgev(size, loc=ui[1], scale=1, shape=0); X2<-rgev(size, loc=ui[2], scale=1, shape=0)
Z1<-rgev(size, loc = ui[1]+bi[1]*tt, scale = 1, shape = 0); Z2<-rgev(size, loc = ui[2]+bi[2]*t, scale = 1, shape = 0)
# p12 computation
p12_1<-p12_NonPar(X1,Z1,tt,tt,0.2); p12_2<-p12_NonPar(X2,Z2,tt,tt,0.2);ic_1<-IC(X1,Z1,tt,tt,0.2) ;ic_2<-IC(X2,Z2,tt,tt,0.2) 
#plot
df_Ex <-data.frame(temp=tt,x=p12_1,y=p12_2,z1=ic_1$low,z2=ic_1$high,z3=ic_2$low,z4=ic_2$high)
p12plot<-ggplot(df_Ex,aes(x=tt,y=p12_1))+ geom_line(colour="red")+geom_ribbon(aes(ymin=ic_1$low, ymax=ic_1$high), linetype=2, alpha=0.1)+
  geom_line(aes(y=p12_2), colour="black")+geom_ribbon(aes(ymin=ic_2$low, ymax=ic_2$high), linetype=2, alpha=0.1)
p12plot

#### Multiple trajectories
size <-  250*4
tt <- seq.int(size)/size
N <- 2 #number of trajectories
ui<-rnorm(N) 
bi<-rnorm(N,mean=1,sd=0.1)
scale_Z_real=1
shape_Z_real=0
Z_real<-Z_real_comp(size,scale=1,shape=0) # "real" factual(Z) and couterfactual(X) run
X_real<-rgev(size*5, loc=mean, scale=1, shape=0)
X1<-rgev(size*5, loc=ui[1], scale=1, shape=0) # run from model 1
Z1<-rgev(size, loc = ui[1]+bi[1]*tt, scale = 1, shape = 0)
X2<-rgev(size*5, loc=ui[2], scale=1, shape=0) # run from model 2
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

################## World before anthropogenic forcing effets #################################

X1<-rgev(size*5, loc=ui[1], scale=1, shape=0); Z1<-rgev(size, loc = ui[1], scale = 1, shape = 0)# Trajectory world 1
X2<-rgev(size*5, loc=ui[2], scale=1, shape=0); Z2<-rgev(size, loc = ui[2], scale = 1, shape = 0)# Trajectory world 2

p12_1<-p12_NonPar(X1,Z1,tt,tt,0.11); p12_2<-p12_NonPar(X2,Z2,tt,tt,0.11)
ic_1<-IC(X1,Z1,tt,tt,0.11) ;ic_2<-IC(X2,Z2,tt,tt,0.11) 

p12hat_moyen<- (p12_1 + p12_2)/2

df_Ex <-data.frame(temp=tt, y1=p12_1, y2=p12_2, y3= p12hat_moyen, z1=ic_1$low,z2=ic_1$high,z3=ic_2$low,z4=ic_2$high)
p12plot<-ggplot(df_Ex,aes(x=tt,y=p12_1))+ geom_line(colour="red")+geom_ribbon(aes(ymin=ic_1$low, ymax=ic_1$high), linetype=2, alpha=0.1)+
  geom_line(aes(y=p12_2), colour="blue")+geom_ribbon(aes(ymin=ic_2$low, ymax=ic_2$high), linetype=2, alpha=0.1) +  geom_line(y=p12hat_moyen, colour="green")+
  geom_hline(yintercept = 0.5) + ggtitle("p02 estimation over time") + ylab("p02") + xlab("time")
p12plot

################################################################################################
#### 2. p13_t observations                                                                  ####
################################################################################################

p13_hat<-p13_NonPar(X,Z,tt,tt,0.11)
plot(tt, p13_hat)
p12_hat<-p12_NonPar(X,Z,tt,tt,0.11)
plot(tt, p12_hat)

dev.off()

################################################################################################
################################################################################################
### Extension to p1r_t                                                                       ###
################################################################################################
################################################################################################
library(fields); library(pracma)

################################################################################################
#### 3. lambda_t & k_t values, using different methods                                      ####
################################################################################################


##################  Tests for fsolve(), optim()  ###############################################
################################################################################################

####  Visualisation of teh solution  

startval_t <- c(0.7, 1.3) 
phat <-c(1/2, 1/3) 
phat <-c(0.86, 0.79) 
fg <- function(theta){functionp12p13(theta,vecx=phat)}

lambda <- seq(0.01, 1.5, 0.01)
k <- seq(0.01, 1.5, 0.01)
p12_mat <- matrix(nrow = length(lambda), ncol = length(k))
p13_mat <- matrix(nrow = length(lambda), ncol = length(k))
for(i in seq_along(lambda)){
  for(j in seq_along(k)){
    p12_mat[i, j] <- fg(c(lambda[i], k[j]))[1]
    p13_mat[i, j] <- fg(c(lambda[i], k[j]))[2]
  }
}
idx12 <- which(abs(p12_mat) == min(abs(p12_mat)), arr.ind =  TRUE);idx12
idx13 <- which(abs(p13_mat) == min(abs(p13_mat)), arr.ind =  TRUE);idx13
image.plot(lambda, k, p12_mat)
points(lambda[idx12[1]], k[idx12[2]], col="white", pch=20)
image.plot(lambda, k, p13_mat)
points(lambda[idx13[1]], k[idx13[2]], col="white", pch=20)

####  3. fsolve() and fsolve_modif() methods ###################################################

# Estimation of lambda_t k_t for only 1 time t
J12 <- function(x){jacobianFunctiong12(x[1],x[2])[[1]]
}
startval_t <- c(0.7, 1.3) 
phat <-c(1/2, 1/3) 
phat <-c(0.86, 0.79) 
fg <- function(theta){functionp12p13(theta,vecx=phat)}

# Using the jacobian
Fsolveweibull <- fsolve(fg,x0=startval_t,J=J12) # Error: C stack usage  7970896 is too close to the limit
# Using the modified function fsolve_modif()
Fsolveweibull_modif <- fsolve_modif(fg,x0=startval_t)


# Estimation of lambda_t, k_t for for two stationary trajectories X,Z
data <- simulWclass(m=100,n=100,N=2,ksiX=-0.20,ksiZ=-0.25,sigX=1,sigZ=1)
tt <- seq.int(100)/100
matp<-P12_P13_estimation(data$matX,data$matZ,tt,tt,11,kern=dEpan) 
theta<-weibullFsolve(matp$matp12,matp$matp13,truevalues=NULL)
theta3<-weibullFsolve_modif(matp$matp12,matp$matp13,truevalues=NULL) 

####  4. Optim() method ########################################################################

# Estimation of lambda_t k_t for only 1 time t
startval_t <- c(0.7, 1.3) 
#startval2_t <- c(0.8, 1.8) 
phat <-c(0.86, 0.79) 
#phat <-c(1/2, 1/3) 

EOptimweibull <- optim(startval_t,fgoptim,method="BFGS",vecx=phat) # Using method BFGS without bound
EOptimweibull <- optim(startval_t,fgoptim,method="L-BFGS-B",vecx=phat,lower = 0.001) # Using methos L-BFGS-B with bound (no negative parameters)

# Estimation of lambda_t, k_t for for stationary trajectories X,Z
data <- simulWclass(m=20,n=20,N=1,ksiX=-0.20,ksiZ=-0.25,sigX=1,sigZ=1)
tt <- seq.int(20)/20
matp<-P12_P13_estimation(data$matX,data$matZ,tt,tt,0.1,kern=dEpan) 
plot(tt,matp$matp13[,1],xlab="t",ylab="p12_hat")

thetahat<-weibullOptim (matp$matp12,matp$matp13,truevalues=NULL)#lam= 0.001 ,k= 1.387237
plot(tt,thetahat$lambdahat[,1], xlab="t",ylab="lambda_hat") # parameter change a llitle bit
plot(tt,thetahat$khat[,1],xlab="t",ylab="k_hat")

##########  Estimation on p1r
# p1.10.t for each trayectory
p1rfarestim <- p1rfarW(thetahat$lambdahat,thetahat$khat,matrix(10,ncol=2,nrow=20)) 
rvalues <- seq(5,30,by=5)
p1rfarchoix <- p1rfarW(rep(thetahat$lambdahat[1],6),rep(thetahat$khat[1],6),rvalues)

## computation of lambda_t, k_t for each time. Z_t non stationary and same support
size=20
tt <- seq.int(size)/size
mu = 2 + seq(0, 5, length.out = size) 
X = rgev(size * 5, loc = 0, scale = 1, shape = 0)
Z = rgev(size, loc = mu, scale = 1, shape = 0)
matp<-P12_P13_estimation(X,Z,tt,tt,0.11,kern=dEpan)
#plot(tt,matp$matp13,type="l", col="blue")
#lines(tt,matp$matp12,type="l")
thetahat<-weibullOptim (matp$matp12,matp$matp13,truevalues=NULL)
rvalues <- seq(5,30,by=5)
##########  Estimation on p1r
p1rfar1<-p1rfarW(rep(thetahat$lambdahat[1],6),rep(thetahat$khat[1],6),rvalues) # t1 and sequence of r
p1rfarestim <- p1rfarW(thetahat$lambdahat,thetahat$khat,matrix(10,ncol=2,nrow=20))

# p1r_t with r=4 for all times 
p1rfar<-p1rfarW_temps(thetahat[[1]],thetahat[[2]],matrix(4,ncol=1,nrow=size))
plot(tt,p1rfar$p1r,type="l",col="blue",xlab="t",ylab="p12_4")
lines(tt,matp$matp12,col="gray")
lines(tt,matp$matp13,col="gray")

################################################################################################
#### 4. testing the method with "known" W-class X and Z                                     ####
################################################################################################

####  Using only 1 trajectory {(X)t,(Z)t} (N=1) ################################################

size <-  200
tt <- seq.int(size)/size
muz = 2 + seq(0, 5, length.out = size)
mux = 0
sigmax = 1
sigmaz = 1
r=10 ; cat("r: ", r )
X = rgev(size * 5, loc = mux, scale = sigmax, shape = 0)
Z = rgev(size, loc = muz, scale = sigmaz, shape = 0)
# lambda_t and k_t theorique 
sigmax_t <- rep(sigmax,size)
sigmaz_t <- rep(sigmaz,size)
k_t_theo <- sigmax_t/sigmaz_t; cat("Thoretic values of k_t: ", k_t_theo )
muz_t <- muz
mux_t <- rep(mux, size)
lambda_t_theo <- exp((mux_t-muz_t)/sigmax_t);cat("Thoretic values of lambda_t: ", lambda_t_theo )
# p12_t and p13_t theorique
theta_theo <- 1 / exp(muz)
p12_t_theo <- 1 / (1 + theta_theo);cat("Thoretic values of p12_t: ", p12_t_theo )
p13_t_theo <- 1 / (1 + 2*theta_theo); cat("Thoretic values of p13_t: ", p13_t_theo  )
p1r_t_theo <- 1 / (1 + (r-1)*theta_theo); cat("Thoretic values of p1r_t: ", p1r_t_theo  )
# p12hat and p13hat
matp<-P12_P13_estimation(X,Z,tt,tt,0.11,kern=dEpan)
p12hat_t<-matp$matp12; cat("Estimated values of p12_t: ", p12hat_t  )
p13hat_t<-matp$matp13; cat("Estimated values of p3r_t: ", p13hat_t  )
err_p12 <- abs(p12_t_theo - p12hat_t); cat("dif p12: ", err_p12  )
err_p13 <- abs(p13_t_theo - p13hat_t); cat("dif p13: ", err_p13  )
# lambdahat_t and khat_t (with optim)
thetahat<-weibullOptim (matp$matp12,matp$matp13,truevalues=NULL)
lambdahat_t <- thetahat$lambdahat; cat("Estimated values of lambda_t: ", lambdahat_t )
khat_t <- thetahat$khat; cat("Estimated values of k_t: ", khat_t  )
optimval_theo<-fgoptim(c(lambda_t_theo,k_t_theo ),c(p12_t_theo,p13_t_theo)) # with theoretic p12t p13t
cat("theoretic min of function: ", optimval_theo )
optimval_hat<-fgoptim(c(lambdahat_t,khat_t ),c(p12_t_theo,p13_t_theo)) # with estiated p12t p13t
cat("estimated min of function: ", optimval_hat )
#[1] optimval_theo: 0.07771465, optimval_hat: 0.07810431
# Comparation between lambdahat_t and lambda_t_theo (same for k)
err_lambda <- abs(lambdahat_t-lambda_t_theo); cat("error in lambda_t estimation: ", err_lambda )
err_k <- abs(khat_t-k_t_theo); cat("error in k_t estimation ", err_k)
# il y a toujours beaucoup d'erreur sur k (qui en theorie est constant)
# computation of p1rhat_t
p1rfar<-p1rfarW_temps(thetahat[[1]],thetahat[[2]],matrix(r,ncol=1,nrow=size)) # r=4
cat("estimation of p1r_t: ", p1rfar$p1r )
# difference between p1r_T theoretic estimated 
err_p1r <- abs(p1rfar$p1r - p1r_t_theo); cat("error in p1r_t estimation: ", err_p1r)
#plot(tt,err_p1r,ylab="p1.5")
plot(tt,p1r_t_theo,type="l",col="red")
lines(tt,p1rfar$p1r)

p1r_hat_noyaux<-p1r_NonPar(X,Z,r,tt,tt,0.11,kern=dEpan)
lines (tt, p1r_hat_noyaux, col="blue")

####  Multiple trajectories {(X)t,(Z)t} (N=1) ##################################################

# Both Gumbel and Z with linear trend un muz

size <-  50
tt <- seq.int(size)/size
muz = 2 + seq(0, 5, length.out = size)
mux = 0
sigmax = 1
sigmaz = 1
testmoyen<-FastTestforp1r_gumbel(tt,0.11,5,20,20,10,rep(sigmax,size),rep(sigmaz,size),rep(mux,size),muz)
testmoyen<-FastTestforp1r_gumbel(tt,0.11,5,50*2,50,16,rep(sigmax,size*2),rep(sigmaz,size),rep(mux,size*2),muz)

plot(tt,testmoyen$matp1r[,1], type="l") # plot of 1 trayectory

# Fonction similar two (2.) but for a type of Frechet, in this case the trajectories are not necessay W-class
############# not ready :
size <-  20
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


size <-  20
tt <- seq.int(size)/size
sigma <- seq(1, 2, length.out = size)
xi = 0.5
testmoyen<-FastTestforp1r_frechet(tt,0.5,5,size*2,size,10,xi,sigma,rep(1,size*2),sigma)

FastTestforp1r_frechet <- function(tt,h,r,m,n,N,xi,sigma.vec,muX.vec,muZ.vec){
  theta_theo <- sigma.vec^(-1 / rep(xi,n))
  p12_theo <- 1 / (1 + theta_theo)
  
  matX <- matrix(nrow=m,ncol=N)
  matZ <- matrix(nrow=n,ncol=N)
  matp <- matrix(nrow=n,ncol=N)
  
  for (j in 1:N){
    matX[,j] <- rgev(m, loc = muX.vec, scale = rep(xi,m), shape = xi)
    matZ[,j] <- rgev(n, loc = sigma.vec, scale = rep(xi,n)*sigma.vec, shape =xi)
  }
  matp<-P12_P13_estimation(matX,matZ,tt,tt,h,kern=dEpan)
  
  thetahat<-weibullOptim (matp$matp12,matp$matp13,truevalues=NULL)
  p1rfar<-p1rfarW_temps(thetahat[[1]],thetahat[[2]],matrix(r,ncol=N,nrow=n))
  p1r_mean <- rowMeans(p1rfar$p1r)
  
  plot(tt,p1r_t_theo,type="l",col="red")
  lines(tt,p1r_mean)
  for (i in 1:N){
    lines(tt,p1rfar$p1r[,i],col="gray")
  }
  return(list("matp1r"=p1rfar$p1r,"matfar"=p1rfar$far,"p1rmean"=p1r_mean))
}


################################################################################################
#### 5. Simulation of W-class trajectories with F = non-stationary                          ####
################################################################################################

size=200
tt <- seq.int(size)/size
sigma <- seq(1, 2, length.out = size)
simul<-simulWclass_nonStationnary(m=size*2,n=size,N=1,ksiX=0.2,ksiZ=0.25,sigX=0.7,muX=0,muZ=sigma,sigZ=0,unknownMu=FALSE) 
G<-ecdf(simul$matX)
plot(tt,G(simul$matZ))


####  Visualisation of de evolution on the trajectories {(X)t,(Z)t} ############################
a<-length(simul$mu_Z)
muZ_t= c(simul$mu_Z[1],simul$mu_Z[a/4],simul$m_Z[a/2], simul$mu_Z[3*a/4], simul$mu_Z[a])
sigZ_t = c(sigma[1],sigma[a/4],sigma[a/2], sigma[3*a/4], sigma[a])
plot(density(simul$matX),col="red")
Z1 = rgev(size, loc = muZ_t[1], scale = sigZ_t[1], shape = 0.25)
lines(density(Z1),col="gray")
Z2 = rgev(size, loc = muZ_t[2], scale = sigZ_t[2], shape = 0.25)
lines(density(Z2),col="green")
Z3 = rgev(size, loc = muZ_t[3], scale = sigZ_t[3], shape = 0.25)
lines(density(Z3),col="blue")
Z4 = rgev(size, loc = muZ_t[4], scale = sigZ_t[4], shape = 0.25)
lines(density(Z4),col="purple")


size=20
tt <- seq.int(size)/size
# ksi!=0 , ksiZ constant, sigZ variable, mu unknown
sigma <- seq(1, 3, length.out = size)
simul1<-simulWclass_nonStationnary(m=size,N=1,ksiX=0.20,sigX=1,muX=0,ksiZ=rep(0.25,size),muZ=NULL,sigZ=sigma,unknown="location")
# ksi!=0 , ksiZ variable, sigZ variable, mu unknown
sigma <- seq(1, 3, length.out = size)
ksi <- seq(0.1,0.2,length.out = size)
simul2<-simulWclass_nonStationnary(m=size,N=1,ksiX=0.1,sigX=1,muX=0,ksiZ=ksi,muZ=NULL,sigZ=sigma,unknown="location")
# ksi!=0 , ksiZ constant, sigZ variable, sig unknown
mu <- seq(1, 2, length.out = size)
simul3<-simulWclass_nonStationnary(m=size,N=1,ksiX=0.20,sigX=1,muX=0,ksiZ=rep(0.25,size),muZ=mu,sigZ=NULL,unknown="scale")
# ksi!=0 , ksiZ variable, sigZ variable, sig unknown
sigma <- seq(1, 3, length.out = size)
ksi <- seq(0.1,0.2,length.out = size)
simul4<-simulWclass_nonStationnary(m=size,N=1,ksiX=0.1,sigX=1,muX=0.2,ksiZ=ksi,muZ=mu,sigZ=NULL,unknown="scale")

# ksi!=0 , muZ constant, sigZ variable, ksi unknown
sigma <- seq(1, 2, length.out = size)
simul5<-simulWclass_nonStationnary(m=size,N=1,ksiX=0.1,sigX=1,muX=0.2,ksiZ=NULL,muZ=1,sigZ=sigma,unknown="shape") # pourquoi lam const?
sigma <- seq(1, 3, length.out = size)
simul6<-simulWclass_nonStationnary(m=size,N=1,ksiX=0.1,sigX=1,muX=0.2,ksiZ=NULL,muZ=1,sigZ=sigma,unknown="shape")
# -> lambda est toujours constante
# ksi!=0 , muZ variable, sigZ variable, ksi unknown
mu <- seq(1, 2, length.out = size)
simul7<-simulWclass_nonStationnary(m=size,N=1,ksiX=0.25,sigX=1,muX=0.2,ksiZ=NULL,muZ=mu,sigZ=sigma,unknown="shape") 

# ksi==0 , muZ variable, sigZ variable, ksi unknown
sigma <- seq(1, 3, length.out = size)
mu <- seq(1, 2, length.out = size)
simul8<-simulWclass_nonStationnary(m=size,N=1,ksiX=0,sigX=1,muX=0.2,ksiZ=0,muZ=mu,sigZ=sigma) 



################################################################################################
################################################################################################
### lambda_t & k_tUsing GMM                                                                  ###
################################################################################################
################################################################################################
library(gmm)

# functions needed for gmm
# matGZ_func, weibullGMM_NonStationaire_startval_1, weibullGMM_NonStationaire, 
# function_gmm_noyaux,laplaceWeibull, funclaplace, dEpan.


################################################################################################
#### 6. p1r_t estimation of both GMN and Optim                                              ####
################################################################################################


####  Using only 1 trajectory {(X)t, (Z)t} (N=1)       #########################################

size <- 200
tt <- seq.int(size)/size
mu = 2 + seq(0, 5, length.out = size) 
X = rgev(size * 5, loc = 0, scale = 1, shape = 0)
Z = rgev(size, loc = mu, scale = 1, shape = 0)
r=10
h=0.11

split.screen(c(2,2))
# GMM
GZ<-matGZ_func(X,Z)
param<-weibullGMM_NonStationaire(GZ, tt, tt, h, kern=dEpan, truevalues=NULL)
p1r_gmm<-p1rfarW_temps(param[[1]],param[[2]],matrix(r,ncol=1,nrow=size))
# Optim
matp<-P12_P13_estimation(X,Z,tt,tt,h,kern=dEpan)
thetahat<-weibullOptim (matp$matp12,matp$matp13,truevalues=NULL) # lam= 1e-05 ,k= 0.2834847 
p1r_opt<-p1rfarW_temps(thetahat[[1]],thetahat[[2]],matrix(r,ncol=1,nrow=size))
# theo
theta_theo <- 1 / exp(mu)
p1r_t_theo <- 1 / (1 + (r-1)*theta_theo)
screen(4) # change to 2,3,4
plot(tt,p1r_t_theo,col="red",type="l")
lines(tt,p1r_gmm$p1r,col="green")
lines(tt,p1r_opt$p1r,col="blue")
# dev.off()

####  Multiple trajectories {(X)t, (Z)t} (N=1)         #########################################

size <- 250
tt <- seq.int(size)/size
mu = 2 + seq(0, 5, length.out = size)
r=4
plotoptim<-FastTestforp1r_gumbel(tt,h=0.11,r=r,m=size*4,n=size,N=16,sigX.vec=1,sigZ.vec=1,muX.vec=0,muZ.vec=mu,methode="optim")
plotgmm<-FastTestforp1r_gumbel(tt,h=0.11,r=r,m=size*4,n=size,N=16,sigX.vec=1,sigZ.vec=1,muX.vec=0,muZ.vec=mu,methode="gmm")

theta_theo <- 1 / exp(mu)
p1r_t_theo <- 1 / (1 + (r-1)*theta_theo)

plot(tt,p1r_t_theo,col="red",type="l")
lines(tt,plotoptim$p1r_mean, col="blue")
lines(tt,plotgmm$p1r_mean,col="green")

################################################################################################
#### 7. Testing GMM approach                                                                ####
################################################################################################

######### lambda_t & k_t  using GMM with different starting value methodology ##################

size <- 200
tt <- seq.int(size)/size
mu = 2 + seq(0, 5, length.out = size) 
X = rgev(size * 5, loc = 0, scale = 1, shape = 0)
Z = rgev(size, loc = mu, scale = 1, shape = 0)

h=0.11
GZ<-matGZ_func(X,Z)

## Starting from lambda_(t-1), k_(t-1) at each t>1: 
param<-weibullGMM_NonStationaire(GZ, tt, tt, h, kern=dEpan, truevalues=NULL)
# system.time(param<-weibullGMM_NonStationaire(GZ, tt, tt, h, kern=dEpan, truevalues=NULL))#lam= 0.002590114 ,k= 1.304888  
# system tye: 3.055 s
# Starting always from the same value
system.time(param2<-weibullGMM_NonStationaire_startval_1(GZ, tt, tt, h, kern=dEpan, truevalues=NULL))
# system time: 14.580 s
# output: matrix lambdahat and khat :lam= 1.649007e-06 ,k= 0.72061, lam= 0.001947633 ,k= 0.7214428 


####  p1r_t  estimation   ######################################################################
r=4
p1r_gmm<-p1rfarW_temps(param[[1]],param[[2]],matrix(r,ncol=1,nrow=size))

######### Deeper comparison of initialisation methods ############################################
##################################################################################################

set.seed(401)
# set.seed(Sys.time())

size <- 200
tt <- seq.int(size)/size
mu = 2 + seq(0, 5, length.out = size) 
X = rgev(size * 5, loc = 0, scale = 1, shape = 0)
Z = rgev(size, loc = mu, scale = 1, shape = 0)
h=0.11
# r=4
k_t_theo <- rep(1,size) # sigx/sigz 
lambda_t_theo <- exp((0-mu)/1) # exp((mux-muz)/sigx)
GZ<-matGZ_func(X,Z)

param<-weibullGMM_NonStationaire(GZ, tt, tt, h, kern=dEpan, truevalues=NULL)
param2<-weibullGMM_NonStationaire_startval_1(GZ, tt, tt, h, kern=dEpan, truevalues=NULL)

#### Observations : We don't have the sames results  
plot(tt,k_t_theo,type="l", main=" k estimation")
lines(tt,param$khat,col="blue")
lines(tt,param2$khat,col="red")

plot(tt,lambda_t_theo,type="l", main="lambda estimation")
lines(tt,param$lambdahat,col="blue")
lines(tt,param2$lambdahat,col="red")

######################### Visualisation of the system's solution ##################

library(fields)

fg_test <- function(theta,vecx){function_gmm_noyaux(theta,vecx=GZ,index=200,tt.vec=tt,t_eval.vec=tt,bandwidth=h)}

lambda <- seq(0.001, 0.15, 0.005)
k <- seq(0.01, 1.1, 0.01)

dif_p12_mat <- matrix(nrow = length(lambda), ncol = length(k))
dif_p13_mat <- matrix(nrow = length(lambda), ncol = length(k))

for(i in seq_along(lambda)){
  for(j in seq_along(k)){
    dif_p12_mat[i, j] <- mean(fg_test(c(lambda[i], k[j]))[,1])
    dif_p13_mat[i, j] <- mean(fg_test(c(lambda[i], k[j]))[,2])
  }
}
idx12 <- which(abs(dif_p12_mat) == min(abs(dif_p12_mat)), arr.ind =  TRUE);idx12
idx13 <- which(abs(dif_p13_mat) == min(abs(dif_p13_mat)), arr.ind =  TRUE);idx13

image.plot(lambda, k,dif_p12_mat)
points(lambda[idx12[1]], k[idx12[2]], col="white", pch=20)
image.plot(lambda, k, dif_p13_mat)
points(lambda[idx13[1]], k[idx13[2]], col="white", pch=20)

###########################  convergence of the method  ##########################################

fg_noyaux <- function(theta,vecx){function_gmm_noyaux(theta,vecx=GZ,index=150,tt.vec=tt,t_eval.vec=tt,bandwidth=h)}

EGMMweibull_NonStationary <- gmm(g=fg_noyaux,
                                 x=GZ,
                                 t0= c(1,1),
                                 optfct="nlminb",
                                 lower=c(10^(-8),10^(-8)),upper=c(Inf,Inf),
                                 onlyCoefficients = FALSE, #avant c'était TRUE
                                 control = list(trace=1 )
                                 
)

EGMMweibull_NonStationary

summary(EGMMweibull_NonStationary)


################################################################################################
#### 8. Visialusation of trajectory distributions evolution in time                         ####
################################################################################################

n=1000
mu = seq(0, 5, length.out = n) 
plotd_time(n,0,mu,1,0)

mu = seq(0, 3, length.out = n) 
plotd_time(n,0,mu,1,0.2)


################################################################################################
#### 9. Simulation of W-class trajectories with non-sationnary Z                            ####
################################################################################################

#### Examples of trayectories###################################################################

n=100
m=200
tt <- seq.int(n)/n
N=1
ksi_x<-0.2
sig_x<-0.7
mu_x<-0

# Example 1
simul<-simulWclass_nonStationnary_general(m,N,ksiX=ksi_x,sigX=1,muX=mu_x,ksiZ=rep(0.25,n),sigZ=seq(1.0, 1.1, length.out = n),unknown="location",graph=FALSE)
# Example 2
simul<-simulWclass_nonStationnary_general(m,N,ksiX=ksi_x,sigX=1,muX=mu_x,ksiZ=seq(0.25, 0.30, length.out = n),sigZ=seq(1.0, 1.1, length.out = n),unknown="location")
# Example 3
simul<-simulWclass_nonStationnary_general(m,N,ksiX=ksi_x,sigX=1,muX=mu_x,ksiZ=seq(0.25, 0.30, length.out = n),sigZ=seq(1.0, 1.5, length.out = n),unknown="location")
# Exemple 4
simul<-simulWclass_nonStationnary_general(m,N,ksiX=ksi_x,sigX=1,muX=mu_x,ksiZ=seq(0.25, 0.30, length.out = n),muZ=seq(1.0, 1.1, length.out = n),unknown="scale")
# Exemple 5
simul<-simulWclass_nonStationnary_general(m,N,ksiX=ksi_x,sigX=1,muX=mu_x,ksiZ=seq(0.20, 0.24, length.out = n),muZ=seq(1.0, 1.1, length.out = n),unknown="scale") 
# Exemple 6
simul<-simulWclass_nonStationnary_general(m,N,ksiX=ksi_x,sigX=1,muX=mu_x,ksiZ=seq(0.20, 0.23, length.out = n),muZ=rep(1.1,n),unknown="scale")
# Exemple 7
simul<-simulWclass_nonStationnary_general(m,N,ksiX=ksi_x,sigX=1,muX=mu_x,sigZ=seq(0.7, 1, length.out = n),muZ=rep(1.1,n),unknown="shape", graph="FALSE")
# Example 8
simul<-simulWclass_nonStationnary_general(m,N,ksiX=ksi_x,sigX=1,muX=mu_x,sigZ=seq(0.7, 1, length.out = n),muZ=seq(0, 3, length.out = n),unknown="shape", graph="FALSE")
# Example 9
simul<-simulWclass_nonStationnary_general(m,N,ksiX=0,sigX=1,muX=0,ksiZ = rep(0,n),sigZ=seq(1, 1.5, length.out = n), muZ=seq(0, 4, length.out = n),unknown="scale")
# Example 10
simul<-simulWclass_nonStationnary_general(m=200,N=1,ksiX=0,sigX=1,muX=0,ksiZ = rep(0,n),sigZ=rep(1,n),muZ=seq(0, 1, length.out = n), graph="FALSE")



################################################################################################
#### Trajectories compuatated from CMIP trajectory parameters                               ####
################################################################################################

##  TMAX #######################################################################################
n <- 100
m <- 100*2

# point 1
X <- rgev(m,loc=25.5,scale=0.12,shape=0.1)
plot(density(X))
Z <- rgev(n,loc=26.23,scale=0.11,shape=0.03)
lines(density(Z))
# point 70 
X <- rgev(m,loc=25.72,scale=0.13,shape=0.13)
plot(density(X))
Z <- rgev(m,loc=26.04,scale=0.12,shape=0.06)
lines(density(Z))
# point 200
X <- rgev(m,loc=25.27,scale=0.25,shape=0.16)
plot(density(X))
Z <- rgev(m,loc=25.67,scale=0.22,shape=0.25)
lines(density(Z))
# point 1000
X <- rgev(m,loc=31.78,scale=0.11,shape=0.29)
plot(density(X))
Z <- rgev(m,loc=32.31,scale=0.11,shape=0.1)
lines(density(Z))
# point 2000
X <- rgev(m,loc=30.20,scale=0.16,shape=0.16)
plot(density(X))
Z <- rgev(m,loc=30.92,scale=0.19,shape=0.46)
lines(density(Z))

##  prMAX #######################################################################################

# point 1 
X <- rgev(m,loc=4.45,scale= 1.42,shape= -0.04)
plot(density(X))
Z <- rgev(n,loc=5.76,scale=1.57,shape= - 0.01)
lines(density(Z))
# point 70
X <- rgev(m,loc=2.24,scale=0.06,shape= -0.08)
plot(density(X))
Z <- rgev(m,loc=3.09,scale=0.88,shape=-0.07)
lines(density(Z))
# point 200 
X <- rgev(m,loc=0.36,scale=0.08,shape=0.04)
plot(density(X))
Z <- rgev(m,loc=0.61,scale=0.15,shape=0.15)
lines(density(Z))
# point 1000 
X <- rgev(m,loc= 10.11,scale= 5.09,shape= -0.03)
plot(density(X))
Z <- rgev(m,loc=6.13,scale=3.31,shape= -0.47)
lines(density(Z))
# point 2000
X <- rgev(m,loc=9.05,scale=2.13,shape=0.07)
plot(density(X))
Z <- rgev(m,loc=9.07,scale=1.89,shape= -0.21) ## shape parameter of different sign
lines(density(Z))

#####  lambda_t  & k_t estimation ##############################################################

n <- 100
m <- 100*2

#### Test 1:
muX <- 30.20
sigX <- 0.16
ksiX <- 0.16

muZ <- seq(30.20, 30.92, length.out = n)
sigZ <- seq(0.16, 0.19, length.out = n)
ksiZ <- seq(0.16, 0.46, length.out = n)

#### Test 2:
muX <- 25.5
sigX <- 0.12
ksiX <- 0.1

muZ <- seq(25.5, 26.23, length.out = n)
sigZ <- seq(0.12, 0.11, length.out = n)
ksiZ <- seq(0.1, 0.03, length.out = n)

##### Estimation 
lambda_t_theo <- ((sigX/ksiX)/(sigZ/ksiZ))^(1/ksiX)
k_t_theo <- ksiX/ksiZ

matZ<-rep(0,n)
matX <- rgev(m, loc=muX, scale= sigX, shape=ksiX)
for(i in 1:n){
  matZ[i] <-rgev(1, loc=muZ[i], scale=sigZ[i], shape=ksiZ[i])}

plot(matZ,type="l")
lines(matX,col="blue")

GZ<-matGZ_func(matX,matZ)
plot(GZ)

tt <- seq.int(length(matZ))

h1_seq<-seq(from=2, to=50,by=1)
error_vector1<-CV_error(matX,matZ,tt,tt,h1_seq,kern= dEpan)
h1_opt <- h1_seq[which(error_vector1 == min(error_vector1))]; h1_opt 

param<-weibullGMM_NonStationaire(GZ, tt, tt, h1_opt, kern=dEpan, truevalues=NULL)
param2<-weibullGMM_NonStationaire_startval_1(GZ, tt, tt, h, kern=dEpan, truevalues=NULL)

plot(tt,k_t_theo,type="l", main=" k estimation")
lines(tt,param$khat,col="blue")
lines(tt,param2$khat,col="red")

plot(tt,lambda_t_theo,type="l", main="lambda estimation")
lines(tt,param$lambdahat,col="blue")
lines(tt,param2$lambdahat,col="red")

