################################################################################################
################################################################################################
###                               Working with CMIP data                                     ###
################################################################################################
################################################################################################

library(ismev)
library(dplyr)
library(ggplot2)
library(evd)

############################## Data visualisation ##############################################
################################################################################################
################################################################################################

#### Latitud and longitude database
df_lonlat <- readRDS("~/Documents/LSCE/datasets/df_lonlat.rds")
View(df_lonlat)

#### CMIP6 TMAX dataset
tmax_cmip6_yearmax <- readRDS("~/Documents/LSCE/datasets/tmax_cmip6_yearmax.rds")
table(tmax_cmip6_yearmax$experiment) # there are three possible scenarios
unique(subset(tmax_cmip6_yearmax, experiment == "rcp85")$year) # years for which the scenario rcp85 was calculated
unique(subset(tmax_cmip6_yearmax, experiment == "historical" & model == "BCC-CSM2-MR")$run) # runs donne with the model BCC-CSM2-MR for the scenario "historical"
# 1984-2014
#### CMIP56 prmax dataset
df_pr_cmip56_yearmax <- readRDS("~/Documents/LSCE/datasets/df_pr_cmip56_yearmax.rds")

#### Example plot 

# Selection of one run of the sceario rcp85 for 1 point in the space
# For 2 different models:between 2015 and 2100
df_tmax_rcp85_BCC_CSM2_MR <- tmax_cmip6_yearmax %>% select (institute, model, experiment, run, year, "1") %>% 
  filter(experiment == "rcp85" & model == "BCC-CSM2-MR")%>%rename(temp=6) %>% arrange(year)

df_tmax_rcp85_IPSL<-tmax_cmip6_yearmax %>% select (institute, model, experiment, run, year, "1") %>% 
  filter(experiment == "rcp85" & model == "IPSL-CM6A-LR" &  run=="r1i1p1f1")%>%rename(temp=6) %>% arrange(year) %>% slice(1:86)

# TMAX rcp85 sceario for BCC_CSM2_MR and IPSL models
v1 <- df_tmax_rcp85_BCC_CSM2_MR$year
v2 <-df_tmax_rcp85_BCC_CSM2_MR$temp
v3 <- df_tmax_rcp85_IPSL$temp
data <- data.frame(v1, v2, v3)
ggplot(data) + geom_line(aes(v1, v2), color = "darkred")+ geom_line(aes(v1, v3), color="steelblue")+
  theme_bw() + theme(legend.position = "left") +xlab(" years")+ylab("temperature")+ggtitle("Rcp85 scenario ")


############################## Computation of p1r for TMAX #############################################
################################################################################################
################################################################################################


##  1. We choose the data ######################################################################
################################################################################################

# Model: BCC-CSM2-MR
# years <- X: all historicalNat trajectory, Z1: 1984 - 2014  from historical, Z2: 2015-2045 from rcp85
# point 6 , lon:-152.5	lat:-87.5

# X
tmax_historicalNat<- tmax_cmip6_yearmax %>% select (institute, model, experiment, run, year, "1") %>% 
  filter(experiment == "historicalNat" & model == "BCC-CSM2-MR") %>% rename(temp=6) %>% arrange(year) %>%
  select(temp)
X <-as.vector(tmax_historicalNat[[1]])/10

# Z1
tmax_rcp85<- tmax_cmip6_yearmax %>% select (institute, model, experiment, run, year, "1") %>% 
  filter(experiment == "rcp85" & model == "BCC-CSM2-MR"& between(year, 2015, 2045)) %>% rename(temp=6) %>% arrange(year) %>%
  select(temp)
Z1 <-as.vector(tmax_rcp85[[1]])/10

# Z2
tmax_historical<- tmax_cmip6_yearmax %>% select (institute, model, experiment, run, year, "1") %>% 
  filter(experiment == "historical" & model == "BCC-CSM2-MR" & between(year, 1984, 2014))%>%rename(temp=6) %>% arrange(year) %>% 
  select(temp)
Z2 <-as.vector(tmax_historical[[1]])/10


tt1 <- seq.int(length(Z1))
tt2 <- seq.int(length(Z2))

##  2. We estimate the parameters of the distributions #########################################
################################################################################################

# 2.1 both stationary

fitX <- gev.fit(X)
ksiX <- fitX$mle[1] # 25.93
sigX <- fitX$mle[2] # 0.12
muX <- fitX$mle[3] # -0.15

fitZ1 <- gev.fit(Z1)
muZ1 <- fitX$mle[1] # 25.93
sigZ1 <- fitX$mle[2] # 0.12
ksiZ1 <- fitX$mle[3] # -0.15

fitZ2 <- gev.fit(Z2)
muZ2 <- fitX$mle[1] # 25.93
sigZ2 <- fitX$mle[2] # 0.12
ksiXZ2 <- fitX$mle[3] # -0.15


#### !!! Aparantement toutes ont les memes paralametres de GEV ce qui est bizarre 
#### C'est parce qu'en 30 ans ça ne change pas beaucoup? 

library(EnvStats)
dataX <- egevd(X)

# location = 25.9335925
# scale    =  0.1289788
# shape    =  0.1528549

dataZ1 <- egevd(Z1)
# location = 26.0828541
# scale    =  0.1038379
# shape    =  0.1009625

dataZ2 <- egevd(Z2)
# location = 26.0828541
# scale    =  0.1038379
# shape    =  0.1009625

#### Q: !!!! Avec cette autre méthode on observe un changement entre X et Z 
#### mais le parametre de forme est positif, ce qui est bizarre parce que pour
#### les temperatures les ksi sont censée d'être négatifs 

# 2.2 Z non statotionary

# trouver comment faire

##  3. Optimal bandwidth #######################################################################
################################################################################################

h1_seq<-seq(from=2, to=31,by=1)
error_vector1<-CV_error(X,Z1,tt1,tt1,h1_seq,kern= dEpan)
h1_opt <- h1_seq[which(error_vector1 == min(error_vector1))]; h1_opt 
# optimal badwidth is h= 30


h2_seq<-seq(from=10, to=31,by=1)
error_vector2<-CV_error(X,Z2,tt2,tt2,h2_seq,kern= dEpan)
h2_opt <- h2_seq[which(error_vector2 == min(error_vector2))]; h2_opt
# optimal badwidth is h= 30

### Q: !!!! Le bandwidth opt est de prendre tout la trajectoire
### cela veut dire que c'est plutot stationaire ? 


##  4. We compute p12t_hat & p13t_hat with kernel method #######################################
################################################################################################

G_emp<-ecdf(X)

# historalNat vs rcp85
GZ1<-G_emp(Z1)
plot(tt1,GZ1)
p12_hat_Z1<-p12_NonPar(X,Z1,tt1,tt1,31) 
p13_hat_Z1<-p13_NonPar(X,Z1,tt1,tt1,31)
par(1,2)
par(mfrow=c(1,2))
plot(tt1, p12_hat_Z1, type="l", col="blue")
plot(tt1, p13_hat_Z1, type="l", col="black")


# historalNat vs historical
GZ2<-G_emp(Z2)
plot(tt1,GZ2)
p12_hat_Z2<-p12_NonPar(X,Z2,tt1,tt1,31) # faut calculer l'optimal bandwidth
p13_hat_Z2<-p13_NonPar(X,Z2,tt1,tt1,31)
par(mfrow=c(1,2))
plot(tt1, p12_hat_Z2, type="l", col="blue")
plot(tt1, p13_hat_Z2, type="l", col="black")

#### !!!! Q: significat de ces résultats

##  4. Theoretical computation of lambda & k paramaters ########################################
################################################################################################

## Pour l'instant on va utiliser les paramètres issues de la deuxieme méthode mais en supposant 
## que les paramètres de forme sont négatifs

k1_theo <- 0.1009625/0.1528549
lam1_theo <- (k1_theo * (0.1289788/0.1289788))^(1/0.1009625) #(k1_theo * (sigZ/sigX))^(-1/ksiX)

k2_theo <- 0.1009625/0.1528549 # 0.660512
lam2_theo <- (k2_theo * (0.1289788/0.1289788))^(1/0.1009625) # 0.016

## Q: lambda est très petit et k constant 

##  4. estimation of lambda and k paramaters ###################################################
################################################################################################

## We will just computate for Z1 (because for both trajectories we estimated the same gev parameters)

# lambdahat_t and khat_t estimated with quasi-newton optimisation method 
thetha_opt<-P12_P13_estimation(X,Z1,tt1,tt1,31,kern=dEpan)
thetahat_opt<-weibullOptim (matp$matp12,matp$matp13,truevalues=NULL)
lamhat_opt <- thetahat_opt[[1]] # vers 0.26
khat_opt <- thetahat_opt[[2]] # vers 1

# lambdahat_t and khat_t estimated using GMM initializing from t-1 values
GmZ1<-matGZ_func(X,Z1)
param<-weibullGMM_NonStationaire(GmZ1, tt1, tt1, 31, kern=dEpan, truevalues=NULL)
lamhat_gmm <- param[[1]] # vers 0.26
khat_gmm <- param[[2]] # vers 1

# lambdahat_t and khat_t estimated using GMM with fixed initialization 
param_startval<-weibullGMM_NonStationaire_startval_1 (GmZ1, tt1, tt1, 31, kern=dEpan, truevalues=NULL)
lamhat_gmm_startval <- param_startval[[1]] # vers 0.26
khat_gmm_startval<- param_startval[[2]] # vers 1

## Q: On ne trouve pas du tout les mêmes paramètres que dans la théorique
## !! La bon chose est que les methodes nous données à peu près les memes résultats
## Peut-être c'est la méthode pour trouver les parametres de la gev la mauvaise??
## peut-être c'est pas bien de suposser que Z est stationnaire?

##  5. p1r Computation using the estimated lambda_t  and k_t ###################################
################################################################################################

r1=4
r2=10

## For r=4  ####################################################################################

# Parameters coming from quasi-newton optimisation method
p1rfar_optim<-p1rfarW_temps(lamhat_opt,khat_opt,matrix(r1,ncol=1,nrow=length(Z1))) 
plot(tt1,p1rfar_optim$p1r, type="l" )
# Parameters coming from gmm method with t-1 initialization
p1rfar_gmm<-p1rfarW_temps(lamhat_gmm,khat_gmm,matrix(r1,ncol=1,nrow=length(Z1))) 
lines(tt1,p1rfar_gmm$p1r,col="blue" )
# Parameters coming from gmm method with fixed initialization 
p1rfar_gmm_startval<-p1rfarW_temps(lamhat_gmm_startval,khat_gmm_startval,matrix(r1,ncol=1,nrow=length(Z1))) 
lines(tt1,p1rfar_gmm_startval$p1r,col="green" )
# Just kernels
p1r_hat_noyaux<-p1r_NonPar(X,Z1,r1,tt1,tt1,31,kern=dEpan)
lines(tt1,p1r_hat_noyaux,col="purple" )

# using theoretic patameters lambda_t , k_t
k1_theo_t <- as.matrix(rep(k1_theo,length(Z1)))
lam1_theo_t <- as.matrix(rep(lam1_theo,length(Z1)))
p1rfar_param_theo<-p1rfarW_temps(lam1_theo.vec,k1_theo.vec,matrix(r1,ncol=1,nrow=length(Z1))) 
lines(tt1,p1rfar_param_theo$p1r,col="gray")
#### Q: Les parameters de la gev sont mal estimés ?

## For r=10 ####################################################################################
# Parameters coming from quasi-newton optimisation method
p1rfar_optim<-p1rfarW_temps(lamhat_opt,khat_opt,matrix(r2,ncol=1,nrow=length(Z1))) 
plot(tt1,p1rfar_optim$p1r, type="l" )
# Parameters coming from gmm method with t-1 initialization
p1rfar_gmm<-p1rfarW_temps(lamhat_gmm,khat_gmm,matrix(r2,ncol=1,nrow=length(Z1))) 
lines(tt1,p1rfar_gmm$p1r,col="blue" )
# Parameters coming from gmm method with fixed initialization 
p1rfar_gmm_startval<-p1rfarW_temps(lamhat_gmm_startval,khat_gmm_startval,matrix(r2,ncol=1,nrow=length(Z1))) 
lines(tt1,p1rfar_gmm_startval$p1r,col="green" )
# Just kernels
p1r_hat_noyaux<-p1r_NonPar(X,Z1,r2,tt1,tt1,31,kern=dEpan)
lines(tt1,p1r_hat_noyaux,col="purple" )

## Q:!! La methode a noyaux est mauvais por r grand?

# theoretic parameters
#p1rfar_param_theo<-p1rfarW_temps(lam1_theo.vec,k1_theo.vec,matrix(r2,ncol=1,nrow=length(Z1))) 
#lines(tt1,p1rfar_param_theo$p1r,col="gray")

