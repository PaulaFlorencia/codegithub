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

p1rfar1<-p1rfarW(rep(thetahat$lambdahat[1],6),rep(thetahat$khat[1],6),rvalues)
p1rfarestim <- p1rfarW(thetahat$lambdahat,thetahat$khat,matrix(10,ncol=2,nrow=1000))

p1rfar<-p1rfarW_temps(thetahat[[1]],thetahat[[2]],matrix(4,ncol=1,nrow=size))
plot(p1rfar$p1r)
lines(matp$matp12)
lines(matp$matp13)
