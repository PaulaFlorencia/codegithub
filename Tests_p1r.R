##################  Tests for fsolve(), optim() and estimation of p1rt ############
###################################################################################

########################### Visualisation of teh solution ######################## 
library(fields)

startval_t <- c(0.7, 1.3) 
phat <-c(1/2, 1/3) 
phat <-c(0.86, 0.79) 
fg <- function(theta){functiong(theta,vecx=phat)}

lambda <- seq(0.01, 2, 0.01)
k <- seq(0.01, 2, 0.01)
p12_mat <- matrix(NA, nrow = length(lambda), ncol = length(k))
p13_mat <- matrix(NA, nrow = length(lambda), ncol = length(k))
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
data <- simulWclass(m=20,n=20,N=2,ksiX=-0.20,ksiZ=-0.25,sigX=1,sigZ=1)
tt <- seq.int(20)/20
matp<-GZestimation(data$matX,data$matZ,tt,tt,11,kern=dEpan) 
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


# B.1
## computation of lambda_t, k_t for each time. X,Z stationary 

data <- simulWclass(m=20,n=20,N=2,ksiX=-0.20,ksiZ=-0.25,sigX=1,sigZ=1)
tt <- seq.int(20)/20
matp<-GZestimation(data$matX,data$matZ,tt,tt,11,kern=dEpan) 
thetahat<-weibullGMMestim4 (matp$matp12,matp$matp13,truevalues=NULL)
# p1.10.t for each trayectory
p1rfarestim <- p1rfarW(thetahat$lambdahat,thetahat$khat,matrix(10,ncol=2,nrow=20)) 
rvalues <- seq(5,30,by=5)
p1rfarchoix <- p1rfarW(rep(thetahat$lambdahat[1],6),rep(thetahat$khat[1],6),rvalues)

# B.2
## computation of lambda_t, k_t for each time. Z_t non stationary and same support

size=1000
tt <- seq.int(size)/size
mu = 2 + seq(0, 5, length.out = size) 
X = rgev(size * 5, loc = 0, scale = 1, shape = 0)
Z = rgev(size, loc = mu, scale = 1, shape = 0)
matp<-GZestimation(X,Z,tt,tt,0.11,kern=dEpan)
thetahat<-weibullGMMestim4 (matp$matp12,matp$matp13,truevalues=NULL)
rvalues <- seq(5,30,by=5)

p1rfar1<-p1rfarW(rep(thetahat$lambdahat[1],6),rep(thetahat$khat[1],6),rvalues) # t1 and sequence of r
p1rfarestim <- p1rfarW(thetahat$lambdahat,thetahat$khat,matrix(10,ncol=2,nrow=1000))
# r=10 for all t.vec. Here we can see that here we are taking in consideration time

p1rfar<-p1rfarW_temps(thetahat[[1]],thetahat[[2]],matrix(4,ncol=1,nrow=size))
plot(p1rfar$p1r)
lines(matp$matp12)
lines(matp$matp13)

#########

################ testing the method with "known" W-class F and CF #############

# 1. for 1 trayactory {(X)t,(Z)t} (N=1)
size <-  100
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
matp<-GZestimation(X,Z,tt,tt,0.11,kern=dEpan)
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
FastTestforp1r_gev <- function(tt,h,r,m,n,N,sigX.vec,sigZ.vec,muX.vec,muZ.vec){
  
  theta_theo <- 1 / exp(muZ.vec)
  p1r_t_theo <- 1 / (1 + (r-1)*theta_theo)
  
  matX <- matrix(nrow=m,ncol=N)
  matZ <- matrix(nrow=n,ncol=N)
  matp <- matrix(nrow=n,ncol=N)
  
  for (j in 1:N){
    matX[,j] <- rgev(m, loc = muX.vec, scale = sigX.vec, shape = 0)
    matZ[,j] <- rgev(n, loc = muZ.vec, scale = sigZ.vec, shape = 0)
  }
  matp<-GZestimation(matX,matZ,tt,tt,h,kern=dEpan)
  
  thetahat<-weibullGMMestim4 (matp$matp12,matp$matp13,truevalues=NULL)
  p1rfar<-p1rfarW_temps(thetahat[[1]],thetahat[[2]],matrix(r,ncol=N,nrow=n))
  p1r_mean <- rowMeans(p1rfar$p1r)
  
  plot(tt,p1r_t_theo,type="l",col="red",xlab="time",ylab="p1r_t", main=paste("evolution over time of p1r with r=",r))
  lines(tt,p1r_mean)
  for (i in 1:N){
    lines(tt,p1rfar$p1r[,i],col="gray")
  }
  return(list("matp1r"=p1rfar$p1r,"p1r_mean"=p1r_mean, "p1rmean"=p1r_mean))
}

size <-  20
tt <- seq.int(size)/size
muz = 2 + seq(0, 5, length.out = size)
mux = 0
sigmax = 1
sigmaz = 1
testmoyen<-FastTestforp1r_gumbel(tt,0.11,5,20,20,10,rep(sigmax,size),rep(sigmaz,size),rep(mux,size),muz)
testmoyen<-FastTestforp1r_gumbel(tt,0.11,5,20*2,20,10,rep(sigmax,size*2),rep(sigmaz,size),rep(mux,size*2),muz)

plot(tt,testmoyen$matp1r[,1], type="l")

# 3. Fonction similar two (2.) but for a type of Frechet, in this case the trajectories are not necessay W-class


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

testmoyen<-FastTestforp1r_gumbel(tt,0.5,5,size*2,size,10,rep(xi,size*2),sigma,rep(1,size*2),sigma)

FastTestforp1r_frechet <- function(tt,h,r,m,n,N,xi.vec,sigma.vec,muX.vec,muZ.vec){
  
  theta_theo <- sigma.vec^(-1 / xi.vec)
  p12_theo <- 1 / (1 + theta_theo)
  
  matX <- matrix(nrow=m,ncol=N)
  matZ <- matrix(nrow=n,ncol=N)
  matp <- matrix(nrow=n,ncol=N)
  
  for (j in 1:N){
    matX[,j] <- rgev(m, loc = muX.vec, scale = xi.vec, shape = xi.vec)
    matZ[,j] <- rgev(n, loc = muZ.vec, scale = xi.vec*muZ.vec, shape = xi.vec)
  }
  matp<-GZestimation(matX,matZ,tt,tt,h,kern=dEpan)
  
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



