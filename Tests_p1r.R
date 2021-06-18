##################  Tests for fsolve(), optim() and estimation of p1rt ############
###################################################################################

########################### Visualisation of teh solution ######################## 
library(fields)

startval_t <- c(0.7, 1.3) 
phat <-c(1/2, 1/3) 
phat <-c(0.86, 0.79) 
fg <- function(theta){functiong(theta,vecx=phat)}

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
#39 186
# 53  29
############### (A) fsolve() ###############################3
library(pracma)

# A.0 We test le fonction that estimates lambda_t k_t for only 1 time
J12 <- function(x){jacobianFunctiong12(x[1],x[2])[[1]]
}
startval_t <- c(0.7, 1.3) 
phat <-c(1/2, 1/3) 
phat <-c(0.86, 0.79) 
fg <- function(theta){functiong(theta,vecx=phat)}

# Using the jacobian
EGMMweibull <- fsolve(fg,x0=startval_t,J=J12) # Error: C stack usage  7970896 is too close to the limit
# Using the modified function fsolve3()
EGMMweibul3 <- fsolve3(fg,x0=startval_t)


# A.1 Computation of lambda_t, k_t for each time. X,Z stationary
data <- simulWclass(m=100,n=100,N=2,ksiX=-0.20,ksiZ=-0.25,sigX=1,sigZ=1)
tt <- seq.int(100)/100
matp<-P12_P13_estimation(data$matX,data$matZ,tt,tt,11,kern=dEpan) 
theta<-weibullGMMestim(matp$matp12,matp$matp13,truevalues=NULL)
theta3<-weibullGMMestim3(matp$matp12,matp$matp13,truevalues=NULL) # ... lam= 336.1192 ,k= -2826.279 error...
# also lam= 0.2871943 ,k= 1.36206 

###############(B) Optim() ####################################

# B.0
# We test le fonction that estimates lambda_t k_t for only 1 time

## lambda0 and k0 to start the algorithm -> startval_t
startval_t <- c(0.7, 1.3) 
#startval2_t <- c(0.8, 1.8) 

## Estimated  p12hat_t1 p13hat_t1 ->vecx
#phat <-c(1/2, 1/3) 
phat <-c(0.86, 0.79) 

#lmb:0.3,k:0.8 p13, p12->lmbd:03,k:1.3

### Using method BFGS without bound
EOptimweibull <- optim(c(0.7, 1.3),fg1,method="BFGS",vecx=c(0.86, 0.79))
#lam= 0.1301583 ,k= 0.7718406 . no convergence

### Using methos L-BFGS-B with bound (no negative parameters)
EOptimweibull <- optim(c(0.7, 1.3),fg1,method="L-BFGS-B",vecx=c(0.86, 0.79),lower = 0.001)
# lam= 0.08297448 ,k= 0.4102397 .  convergence
EOptimweibull <- optim(c(0.5, 0.8),fg1,method="L-BFGS-B",vecx=c(0.86, 0.79),lower = 0.001)

# B.1
## computation of lambda_t, k_t for each time. X,Z stationary 

data <- simulWclass(m=20,n=20,N=1,ksiX=-0.20,ksiZ=-0.25,sigX=1,sigZ=1)
tt <- seq.int(20)/20
matp<-P12_P13_estimation(data$matX,data$matZ,tt,tt,0.1,kern=dEpan) 
plot(tt,matp$matp13[,1],xlab="t",ylab="p12_hat")

thetahat<-weibullGMMestim4 (matp$matp12,matp$matp13,truevalues=NULL)#lam= 0.001 ,k= 1.387237
plot(tt,thetahat$lambdahat[,1], xlab="t",ylab="lambda_hat") # parameter change a llitle bit
plot(tt,thetahat$khat[,1],xlab="t",ylab="k_hat")

# p1.10.t for each trayectory
p1rfarestim <- p1rfarW(thetahat$lambdahat,thetahat$khat,matrix(10,ncol=2,nrow=20)) 
rvalues <- seq(5,30,by=5)
p1rfarchoix <- p1rfarW(rep(thetahat$lambdahat[1],6),rep(thetahat$khat[1],6),rvalues)


# B.2
## computation of lambda_t, k_t for each time. Z_t non stationary and same support

size=20
tt <- seq.int(size)/size
mu = 2 + seq(0, 5, length.out = size) 
X = rgev(size * 5, loc = 0, scale = 1, shape = 0)
Z = rgev(size, loc = mu, scale = 1, shape = 0)
matp<-P12_P13_estimation(X,Z,tt,tt,0.11,kern=dEpan)
#plot(tt,matp$matp13,type="l", col="blue")
#lines(tt,matp$matp12,type="l")
thetahat<-weibullGMMestim4 (matp$matp12,matp$matp13,truevalues=NULL)
rvalues <- seq(5,30,by=5)
p1rfar1<-p1rfarW(rep(thetahat$lambdahat[1],6),rep(thetahat$khat[1],6),rvalues) # t1 and sequence of r
p1rfarestim <- p1rfarW(thetahat$lambdahat,thetahat$khat,matrix(10,ncol=2,nrow=20))
# r=10 for all t.vec. Here we can see that here we are taking in consideration time

p1rfar<-p1rfarW_temps(thetahat[[1]],thetahat[[2]],matrix(4,ncol=1,nrow=size))
plot(tt,p1rfar$p1r,type="l",col="blue",xlab="t",ylab="p12_4")
lines(tt,matp$matp12,col="gray")
lines(tt,matp$matp13,col="gray")

#########

################ testing the method with "known" W-class F and CF #############

# 1. for 1 trayactory {(X)t,(Z)t} (N=1)
size <-  200
tt <- seq.int(size)/size
muz = 2 + seq(0, 5, length.out = size)
mux = 0
sigmax = 1
sigmaz = 1
r=10 ; cat("r: ", r )
X = rgev(size * 5, loc = mux, scale = sigmax, shape = 0)
Z = rgev(size, loc = muz, scale = sigmaz, shape = 0)
# 2. lambda_t and k_t theorique 
sigmax_t <- rep(sigmax,size)
sigmaz_t <- rep(sigmaz,size)
k_t_theo <- sigmax_t/sigmaz_t; cat("Thoretic values of k_t: ", k_t_theo )
muz_t <- muz
mux_t <- rep(mux, size)
lambda_t_theo <- exp((mux_t-muz_t)/sigmax_t);cat("Thoretic values of lambda_t: ", lambda_t_theo )
# 3. p12_t and p13_t theorique
theta_theo <- 1 / exp(muz)
p12_t_theo <- 1 / (1 + theta_theo);cat("Thoretic values of p12_t: ", p12_t_theo )
p13_t_theo <- 1 / (1 + 2*theta_theo); cat("Thoretic values of p13_t: ", p13_t_theo  )
p1r_t_theo <- 1 / (1 + (r-1)*theta_theo); cat("Thoretic values of p1r_t: ", p1r_t_theo  )
# 4. p12hat and p13hat
matp<-P12_P13_estimation(X,Z,tt,tt,0.11,kern=dEpan)
p12hat_t<-matp$matp12; cat("Estimated values of p12_t: ", p12hat_t  )
p13hat_t<-matp$matp13; cat("Estimated values of p3r_t: ", p13hat_t  )
err_p12 <- abs(p12_t_theo - p12hat_t); cat("dif p12: ", err_p12  )
err_p13 <- abs(p13_t_theo - p13hat_t); cat("dif p13: ", err_p13  )
# 5. lambdahat_t and khat_t (with optim)
thetahat<-weibullGMMestim4 (matp$matp12,matp$matp13,truevalues=NULL)
lambdahat_t <- thetahat$lambdahat; cat("Estimated values of lambda_t: ", lambdahat_t )
khat_t <- thetahat$khat; cat("Estimated values of k_t: ", khat_t  )
optimval_theo<-fg1(c(lambda_t_theo,k_t_theo ),c(p12_t_theo,p13_t_theo)) # with theoretic p12t p13t
cat("theoretic min of function: ", optimval_theo )
optimval_hat<-fg1(c(lambdahat_t,khat_t ),c(p12_t_theo,p13_t_theo)) # with estiated p12t p13t
cat("estimated min of function: ", optimval_hat )
#[1] optimval_theo: 0.07771465, optimval_hat: 0.07810431
# 6. Comparation between lambdahat_t and lambda_t_theo (same for k)
err_lambda <- abs(lambdahat_t-lambda_t_theo); cat("error in lambda_t estimation: ", err_lambda )
err_k <- abs(khat_t-k_t_theo); cat("error in k_t estimation ", err_k)
# il y a toujours beaucoup d'erreur sur k (qui en theorie est constant)
# 7. computation of p1rhat_t
p1rfar<-p1rfarW_temps(thetahat[[1]],thetahat[[2]],matrix(r,ncol=1,nrow=size)) # r=4
cat("estimation of p1r_t: ", p1rfar$p1r )
# 8. difference between p1r_T theoretic estimated 
err_p1r <- abs(p1rfar$p1r - p1r_t_theo); cat("error in p1r_t estimation: ", err_p1r)
#plot(tt,err_p1r,ylab="p1.5")
plot(tt,p1r_t_theo,type="l",col="red")
lines(tt,p1rfar$p1r)


# 2. for multiple trajectories {(X)t,(Z)t}

# Both Gumbel and linear trend un muz. We know two Gumbel have the same support

size <-  50
tt <- seq.int(size)/size
muz = 2 + seq(0, 5, length.out = size)
mux = 0
sigmax = 1
sigmaz = 1
testmoyen<-FastTestforp1r_gumbel(tt,0.11,5,20,20,10,rep(sigmax,size),rep(sigmaz,size),rep(mux,size),muz)
testmoyen<-FastTestforp1r_gumbel(tt,0.11,5,50*2,50,16,rep(sigmax,size*2),rep(sigmaz,size),rep(mux,size*2),muz)

plot(tt,testmoyen$matp1r[,1], type="l") # 1 trayectory

# 3. Fonction similar two (2.) but for a type of Frechet, in this case the trajectories are not necessay W-class

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
  
  thetahat<-weibullGMMestim4 (matp$matp12,matp$matp13,truevalues=NULL)
  p1rfar<-p1rfarW_temps(thetahat[[1]],thetahat[[2]],matrix(r,ncol=N,nrow=n))
  p1r_mean <- rowMeans(p1rfar$p1r)
  
  plot(tt,p1r_t_theo,type="l",col="red")
  lines(tt,p1r_mean)
  for (i in 1:N){
    lines(tt,p1rfar$p1r[,i],col="gray")
  }
  return(list("matp1r"=p1rfar$p1r,"matfar"=p1rfar$far,"p1rmean"=p1r_mean))
}

################ Simulation of W-class trajectories with F = non-stationary ###########3

size=200
tt <- seq.int(size)/size
sigma <- seq(1, 2, length.out = size)

simul<-simulWclass_nonStationnary(m=size*2,n=size,N=1,ksiX=0.2,ksiZ=0.25,sigX=0.7,muX=0,muZ=sigma,sigZ=0,unknownMu=FALSE) 
# obs: we still need to add a the fonctionality for testing if all parameters given creater W-class trajectory or not
# And also define warnings or errors in the focntion because here we are assuming the its goning to be "well used"

G<-ecdf(simul$matX)
plot(tt,G(simul$matZ))


# plot in time
# Fast visualisation of de evolution on the trajectories {(X)t,(Z)t}
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


###################### lambdahat_t and khat_t using gmm ###################

library(gmm)

size <- 200
tt <- seq.int(size)/size
mu = 2 + seq(0, 5, length.out = size) 
X = rgev(size * 5, loc = 0, scale = 1, shape = 0)
Z = rgev(size, loc = mu, scale = 1, shape = 0)

h=0.11

GZ<-matGZ_func(X,Z)

param<-weibullGMM_NonStationaire(GZ, tt, tt, h, kern=dEpan, truevalues=NULL)
# output: matrix lambdahat and khat :lam= 1.649007e-06 ,k= 0.72061
p1r_gmm<-p1rfarW_temps(param[[1]],param[[2]],matrix(r,ncol=1,nrow=size))


######### Comparison between de p1r_t estimation of both gmm an optim   #################

# (This quick comparisons are done with this Gumbels because we know they are W-class and we know p1r theorique)

####### Using only 1 trajectory {(X)t, (Z)t} (N=1)
split.screen(c(2,2))
size <- 200
tt <- seq.int(size)/size
mu = 2 + seq(0, 5, length.out = size) 
X = rgev(size * 5, loc = 0, scale = 1, shape = 0)
Z = rgev(size, loc = mu, scale = 1, shape = 0)
r=10
h=0.11

# GMM
GZ<-matGZ_func(X,Z)
param<-weibullGMM_NonStationaire(GZ, tt, tt, h, kern=dEpan, truevalues=NULL)
# output: matrix lambdahat and khat :lam= 1.649007e-06 ,k= 0.72061
p1r_gmm<-p1rfarW_temps(param[[1]],param[[2]],matrix(r,ncol=1,nrow=size))
# Optim
matp<-P12_P13_estimation(X,Z,tt,tt,h,kern=dEpan)
thetahat<-weibullGMMestim4 (matp$matp12,matp$matp13,truevalues=NULL) # lam= 1e-05 ,k= 0.2834847 
p1r_opt<-p1rfarW_temps(thetahat[[1]],thetahat[[2]],matrix(r,ncol=1,nrow=size))
# theo
theta_theo <- 1 / exp(mu)
p1r_t_theo <- 1 / (1 + (r-1)*theta_theo)
screen(1) # change for 2,3,4
plot(tt,p1r_t_theo,col="red",type="l")
lines(tt,p1r_gmm$p1r,col="green")
lines(tt,p1r_opt$p1r,col="blue")
# dev.off()

######## Multiple trajectories {(X)t, (Z)t} (N=1) , comparison

# attention, les vales de X,Z seront different à chaque fois quòn fait toure la fonction 
# les aleurs que prennent les deux méthodes seront differents
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


