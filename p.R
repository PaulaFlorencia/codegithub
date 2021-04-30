###################### Exemples ###############################
library(evd)
# RNG
set.seed(1)

# Example 1
# Gumbel

# Simulations
size <-  250 * 4
rp <- 50
t <- seq.int(size)/size
mu = 2 + seq(0, 5, length.out = size) 
X = rgev(size * 5, loc = 0, scale = 1, shape = 0)
Z = rgev(size, loc = mu, scale = 1, shape = 0)
theta_theo <- 1 / exp(mu)
p12_theo <- 1 / (1 + theta_theo)

##########################################################################
# Example 2
# Frechet
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

#############################################################################
# Epanechnikov density distribution
depan <- function(x){d <- (3 / 4) * (1 - x^2)}

##########################################
p12_NonPar <- function(X,Z,t,h,m=512, k=dnorm){
  
  # Compute empirical cdf of NAT
  G_emp <- ecdf(X)
  G_Z <- G_emp(Z)
  
  # evaluation vector
  rg <- range(t)
  t_eval <- seq(rg[1] - 3 * h, rg[2] + 3 * h, length.out = m)
  kt <- numeric(m)
  
  # Numerator Kt
  const <- sqrt(2 * pi) * h * length(t)
  Kt <- sapply(t, function(tj) k((t_eval - tj)/h))
  
  # Ponderation
  W <- Kt / rowSums(Kt)
  
  p12 <- W %*% G_Z
  
  drop(p12)
  
}