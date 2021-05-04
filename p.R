###################### Exemples ###############################
library(evd)
# RNG
#set.seed(1)
set.seed(Sys.time())

# Easy distributions
size <- 20
t <- seq.int(size)
mu = seq(0, 5, length.out = size) 
X_gev = rgev(size, loc = 0, scale = 1, shape = 0)
Z_gev = rgev(size, loc = 2, scale = 1, shape = 0)
Z_gev_trend = rgev(size, loc = mu, scale = 1, shape = 0)

X_norm = rnorm(size,sd=2)
Z_norm = rnorm(size,0.4,2)

############# Example 1
# Gumbel - (set.seed(292))

# Simulations
size <-  250 * 4
rp <- 50
t <- seq.int(size)/size
mu = 2 + seq(0, 5, length.out = size) 
X = rgev(size * 5, loc = 0, scale = 1, shape = 0)
Z = rgev(size, loc = mu, scale = 1, shape = 0)
theta_theo <- 1 / exp(mu)
p12_theo <- 1 / (1 + theta_theo)

############ Example 2
# Frechet - (set.seed(304))
# Simulations
size <-  250 * 4
rp <- 50
t <- seq.int(size)/size
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
  tt <- (3/4)* (1-x^2)
  tt[-1>x] <- 0
  tt[x>1] <- 0
  return (tt)
} 
# Normal distribution
Kern_gauss <- function(x)exp(-x^2/2)

###################### p12 Non parametric computation ################
p12_NonPar <- function(X,Z,t,t_eval,h,kern= dEpan){
  # G_Z computation
  G_emp <- ecdf(X)
  G_Z <- G_emp(Z)
  # W computation
  Kij <- outer(t_eval,t,function(zz,z) kern((zz - z) / h))
  W <- Kij / rowSums(Kij)
  # p12_hat computation
  p12_hat <- W %*% G_Z
}

# plot p12_theo & p12_hat
p12_hat <- p12_NonPar(X_gev,Z_gev_trend,t,t,0.2,)
plot(t,p12_theo, type = "l", lwd = 4, xlab = "t", ylab = "p12")
lines(t,p12_hat, col = "red", lwd = 2)


######################### Computation of IC ###################################

# Computation of variance
G_emp <- ecdf(X)
G_Z <- G_emp(Z)

VarA <- size # or just size
t <- seq.int(size)/size
t_eval <- seq.int(size)/size

Khj <- outer(t_eval, t, function(zz,z) dEpan((zz - z) / 0.2))
Khi <- outer(t_eval, t, function(zz,z) dEpan((zz - z) / 0.2))
E <- mean(outer(G_Z,G_Z, pmin) - outer(G_Z,G_Z,"*"))

VarB_num <- as.numeric(rowSums(Khj) %*% rowSums(Khi))
VarB_denom <- length(t)*(rowSums(Khj))^2
VarB <- VarB_num*E/VarB_denom

sigma_m <- VarA + VarB

# Computation CI
s <- (sigma_m)^1/2
error <- qnorm(0.975)*s/sqrt(size)
left <- p12_hat-error
right <- p12_hat+error

left_ex <-p12_hat-0.05 #test pour voir si le plot marche bien
right_ex<-p12_hat+0.05

##### plot p12_hat with CI
library(ggplot2)
df <-data.frame(x=t,y=p12_hat,z1=left, z2=right)
p12plot<-ggplot(df,aes(x=t,y=p12_hat)) + geom_line()
p12plot <- p12plot + geom_ribbon(aes(ymin=left_ex, ymax=right_ex), linetype=2, alpha=0.1)
p12plot

dev.off()
