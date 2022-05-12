source("p.R") 
source("Methods_for_p1r_functions.R")
tmax_cmip6_yearmax <- readRDS("/Users/pgonzale/Documents/datasets/tmax_cmip6_yearmax.rds")
pr_cmip6_yearmax <- readRDS("/Users/pgonzale/Documents/datasets/df_pr_cmip56_yearmax.rds")
##########################################################################################################################################
##########################################################################################################################################
##                                       Application of Ws method on CMIP simulations                                                   ##
##########################################################################################################################################

##########################################################################################################################################
##########################################################################################################################################
############################################# 0. Stationary and non-stationary pG  #######################################################

# Yearly maxima of daily maxima of temperature:
##########################################################################################################################################
##########################################################################################################################################

# Examples
##########################################################################################################################################
mat_X_Z <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(200), model.choice="CanESM5", run.choice="r1i1p1f1",var="tmax")
matx <- mat_X_Z$matX; matz <- mat_X_Z$matZ[220:250]; tt <- mat_X_Z$time.vec[220:250]

mat_X_Z <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(700), model.choice="CanESM5", run.choice="r1i1p1f1",var="tmax")
matx <- mat_X_Z$matX; matz <- mat_X_Z$matZ; tt <- mat_X_Z$time.vec

mat_X_Z <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(300), model.choice="IPSL-CM6A-LR", run.choice="r1i1p1f1",var="tmax")
matx <- mat_X_Z$matX; matz <- mat_X_Z$matZ; tt <- mat_X_Z$time.vec

mat_X_Z <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(1345), model.choice="BCC-CSM2-MR", run.choice="r1i1p1f1",var="tmax")
matx <- mat_X_Z$matX; matz <- mat_X_Z$matZ; tt <- mat_X_Z$time.vec

mat_X_Z <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(1345), model.choice="BCC-CSM2-MR", run.choice="r1i1p1f1",var="tmax")
matx <- mat_X_Z$matX; matz <- mat_X_Z$matZ[180:210]; tt <- mat_X_Z$time.vecZ[180:210]

# pG estimation
##########################################################################################################################################
I <- dim(matx)[1]
## Stationary pG
matGZ_ecdf <- GZestimation(matx,matz)
pG_hat<- p1r_NonP(matGZ_ecdf,I) #function stationaire
cat("With stationaroty assumtion, pG = ", pG_hat) 
## Non-stationary pG
hopt <- 35 
pGt_hat <- p1r_NonPar(matx,matz,I,tt,tt,hopt,kern= dEpan) #function non stationaire
plot(tt, pGt_hat, col="blue",xlab="years", ylab="pG", main="pG evolution")
min(pGt_hat);pGt_hat[150];pGt_hat[200];max(pGt_hat)
cat("With Non-stationary assumtion, pG in 1850 is  ", pGt_hat[1], " and in 2100 is ",pGt_hat[250])

# Yearly maxima of daily precipitation
##########################################################################################################################################
##########################################################################################################################################
mat_X_Zpr <- traj_from_data(variable.df=pr_cmip6_yearmax, grid_points=200, model.choice="CNRM-CM5", run.choice="r1i1p1",var="pr")
matxpr <- mat_X_Zpr$matX
matzpr <- mat_X_Zpr$matZ
ttpr <- mat_X_Zpr$time.vec
Ipr <- dim(matxpr)[1]
# Stationary pG
matGZ_ecdfpr <- GZestimation(matxpr,matzpr)
pG_hatpr <- p1r_NonP(matGZ_ecdfpr,Ipr) 
cat("With stationaroty assumtion, pG = ", pG_hatpr) 
# Non-stationary pG
hoptpr <- 35
pGt_hatpr <- p1r_NonPar(matxpr,matzpr,Ipr,ttpr,ttpr,hoptpr,kern= dEpan)
plot(ttpr, pGt_hatpr, col="blue",xlab="years", ylab="pG", main="pG evolution")
min(pGt_hatpr); pGt_hatpr[150]; pGt_hatpr[200]; max(pGt_hatpr)
cat("With Non-stationary assumtion, pG in 1850 is  ", pGt_hatpr[1], " and in 2100 is ",pGt_hatpr[250])

# Non stationary pG plot for both climate variables
##########################################################################################################################################
tyears <- tt + 1850 
plot(tyears, pGt_hat,type="l",xlab="years", ylab=bquote(p[G]), main=bquote(p[G]~ "evolution"), col="orange")
lines(tyears, pGt_hatpr,xlab="years", ylab="pG", col="steelblue")
legend("topleft", legend=c("Tmax", "Prmax"),col=c("orange", "steelblue"), lty=1, cex=1.1)

##########################################################################################################################################
##########################################################################################################################################
######################################## 1. p1rWs with STATIONARITY assumption on Z  #####################################################
r <- 100; rvec <- c(10:100) 

# Data
##########################################################################################################################################
##########################################################################################################################################

# Point 200 - CanESM5
XZ <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(200), model.choice="CanESM5", run.choice="r1i1p1f1",var="tmax")
X <- XZ$matX; Z <- XZ$matZ[1:30]; tt <- XZ$time.vec[1:30]

XZ <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(200), model.choice="CanESM5", run.choice="r1i1p1f1",var="tmax")
X <- XZ$matX; Z <- XZ$matZ[100:130]; tt <- XZ$time.vec 
XZ <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(200), model.choice="CanESM5", run.choice="r1i1p1f1",var="tmax")
X <- XZ$matX; Z <- XZ$matZ[180:210]; tt <- XZ$time.vec

XZ <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(200), model.choice="CanESM5", run.choice="r1i1p1f1",var="tmax")
X <- XZ$matX; Z <- XZ$matZ[210:240]; tt <- XZ$time.vec

# Point 300 - IPSL-CM6A-LR
XZ <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(300), model.choice="IPSL-CM6A-LR", run.choice="r1i1p1f1",var="tmax")
X <- XZ$matX; Z <- XZ$matZ[1:50]; tt <- XZ$time.vec 

XZ <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(300), model.choice="IPSL-CM6A-LR", run.choice="r1i1p1f1",var="tmax")
X <- XZ$matX; Z <- XZ$matZ[100:150]; tt <- XZ$time.vec

XZ <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(300), model.choice="IPSL-CM6A-LR", run.choice="r1i1p1f1",var="tmax")
X <- XZ$matX; Z <- XZ$matZ[170:200]; tt <- XZ$time.vec 

XZ <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(300), model.choice="IPSL-CM6A-LR", run.choice="r1i1p1f1",var="tmax")
X <- XZ$matX; Z <- XZ$matZ[220:250]; tt <- XZ$time.vec

# Point 1345 BCC-CSM2-MR
XZ <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(1345), model.choice="BCC-CSM2-MR", run.choice="r1i1p1f1",var="tmax")
X <- XZ$matX; Z <- XZ$matZ[200:230]; tt <- XZ$time.vec[200:230] 

XZ <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(1345), model.choice="BCC-CSM2-MR", run.choice="r1i1p1f1",var="tmax")
X <- XZ$matX; Z <- XZ$matZ[220:250]; tt <- XZ$time.vec 

XZ <- traj_from_data(variable.df=tmax_cmip6_yearmax, grid_points=c(700), model.choice="CanESM5", run.choice="r1i1p1f1",var="tmax")
X <- XZ$matX; Z <- XZ$matZ[190:220]; tt <- XZ$time.vec[190:220]

I <- length(X); J <- length(Z)

# Non-parametric estimation method (stationary)
########################################################################################################################################## 
matGZ_ecdf <- GZestimation(X,Z,methodGhatm="ecdf")
p1r_nonpar <- p1r_NonP(matGZ_ecdf,I)
cat("p1,100 Nonpar estimation = ", p1r_nonpar) 
p1r_nonpar.vec <- p1r_NonP(matGZ_ecdf,rvec)
cat("p1,10 Nonpar estimation = ", p1r_nonpar.vec[1]," p1,100 Nonpar estimation= ", p1r_nonpar.vec[91] ) 

### Estimation under Wclass assumtion (stationary)
########################################################################################################################################## 
Wclass_theta <- weibullGMMestim(matGZ_ecdf) # we use the ecdf computated before
lamW_hat <- Wclass_theta$lambdahat; 
kW_hat <- Wclass_theta$khat
p1rWclass_hat <- p1rfarW(lamW_hat,kW_hat,r)$p1r
cat("p(1,100)W estimation = ", p1rWclass_hat ,"\n") 
p1rWclass_hat.vec <- p1rfarW(rep(lamW_hat,length(rvec)),rep(kW_hat,length(rvec)),rvec)$p1r # for each r between 10 and 100. We will plot it later (in a comparision plot)
cat("p(1,10)W estimation = ", p1rWclass_hat.vec[1]," p(1,100)W estimation = ", p1rWclass_hat.vec[91] ) 

# Estimation using Ws method (stationary)
########################################################################################################################################## 
x_max <- max(X) # Here we no longer know the actual support su we estimate it as max(x)
Zs_hat <- Z[ Z < x_max]
Ns_hat <- length(Zs_hat)

# N_G_hat <-J-Ns_hat
# p_G_hat2 <- (N_G_hat)/J
# cat("p_G_hat = P(Z>x_max) = ", p_G_hat2 )

matGZ_ecdf <- GZestimation(X,Z,methodGhatm="ecdf")
p_G_hat <- p1r_NonP(matGZ_ecdf,I)
p_G_hat <- as.numeric(p_G_hat)
cat("p_G_hat = P(Z>x_max) = ", p_G_hat ) 
matGZs_ecdf <- GZestimation(X,Zs_hat, methodGhatm="ecdf") 
Wsclass_theta <- WsGMMestim(matGZs_ecdf,p_G_hat,truevalues=NULL) # modifier la fonction pour pouvoir estimer pour les colonnes où pG n'est pas 1 
                                                                 # (dans le cas où chaque colonne est un point different)
ksiFWs_hat <- Wsclass_theta$xiFhat;
ksiGWs_hat <- Wsclass_theta$xiGhat; 
lamWs_hat <- Wsclass_theta$lambdahat
cat("parametres estimés (méthode Ws): "," xiF = ",ksiFWs_hat," xiG = ",ksiGWs_hat, " lambda = ",lamWs_hat, " pG = ",p_G_hat,"\n" ) 
p1rWsclass_andstar_hat <- p1rWs_class(r, lamWs_hat,p_G_hat, ksiGWs_hat, ksiFWs_hat)
p1rWsclass_hat <- p1rWsclass_andstar_hat$p1rWs
p1rstar_hat <- p1rWsclass_andstar_hat$p1rstar
cat("p(1,100)Ws estimation = ", p1rWsclass_hat,"\n") 
p1rWsclass_andstar_hat.vec <- p1rWs_class(rvec, lamWs_hat, p_G_hat, ksiGWs_hat, ksiFWs_hat)
p1rWsclass_hat.vec <- p1rWsclass_andstar_hat.vec$p1rWs  
p1rstar_hat.vec <- p1rWsclass_andstar_hat.vec$p1rstar
p1rWs_prop2_hat <- p1rWs_gamma(r, lamWs_hat,ksiGWs_hat, ksiFWs_hat, p_G_hat)
p1rWs_prop2_hat.vec <- p1rWs_gamma(as.numeric(rvec), lamWs_hat,ksiGWs_hat, ksiFWs_hat, p_G_hat)

# Comparison of estimations of p1r
########################################################################################################################################## 
cat("for r = ",r,",   p1rW_hat: ",p1rWclass_hat," - p1rWs_hat: ",p1rWsclass_hat," - p1r_nonParam: ",p1r_nonpar,"\n",
    " and p_Ghat = ",p_G_hat,"\n") # using r = 100
cat("p1r*: ",p1rstar_hat,"   -     p1r_prop2: ",p1rWs_prop2_hat,"\n") 
plot(p1r_nonpar.vec,type="l",ylim=c(0,1),xlab="r",ylab="p1r",main = " p1r estim ")
lines(p1rWclass_hat.vec,col="blue")
lines(p1rWsclass_hat.vec,col="green")
legend("left",legend=c("p1rNonPar","p1rW_hat", "p1rWs_hat"), 
       col=c("black", "blue", "green"),lty=c(1,1,1,1),cex=0.8)

# Comparison between estimations p1r* and prop2_p1rstar 
########################################################################################################################################## 
plot(p1rstar_hat.vec,type="l",ylim=c(0,1),xlab="r",ylab="p1r",main = " p1r estim ") 
lines(p1rWs_prop2_hat.vec,col="green")
legend("topright",legend=c("p1rNonPar","p1rW_hat", "p1rWs_hat"), 
       col=c("black", "green"),lty=c(1,1),cex=0.8)

##########################################################################################################################################
##########################################################################################################################################
######################################## 2. p1rWs with Non-stationary assumption on Z  ###################################################
r <- 100; rvec <- c(10:100) 

p12_hatNonpar <- p1r_NonPar(matx,matz,2,tt,tt,hopt,kern= dEpan) # creer une function comme matp<-P12_P13_estimation(X,Z,tt,tt,0.11,kern=dEpan)
p13_hatNonpar <- p1r_NonPar(matx,matz,3,tt,tt,hopt,kern= dEpan)
p14_hatNonpar <- p1r_NonPar(matx,matz,4,tt,tt,hopt,kern= dEpan)
# ...

#
#
#
##########################################################################################################################################
##########################################################################################################################################
# Work in progress
# pG maps
##########################################################################################################################################
# map_X_Z <- map_from_data(variable.df=tmax_cmip6_yearmax, model.choice="CanESM5", run.choice="r1i1p1f1",var="tmax", savemat=FALSE)
# mapx <- map_X_Z$matX
# mapz <- map_X_Z$matZ
# tt <- map_X_Z$time.vec
# I <- dim(mapx)[1]
# J <- dim(mapz)[1]
# hopt <- 35 
# pGt_hat.mat <- p1r_NonPar(mapx,mapz,rep(I,J),tt,tt,hopt,kern= dEpan)
# # CA MARCHE PAS, JE NE COMPRENDS PAS POURQUOI
# pGvec <- pGt_hat.mat[200,]
# gridmat_pG <- matrix(data=pGvec,nrow=72,ncol=36)
# grid_lon<-seq(-177.5, 177.5, length.out = 72) 
# grid_lat<-seq(-87.5, 87.5, length.out = 36)
# image.plot(grid_lon,grid_lat,gridmat_pG)
# world(lwd = 1,add=TRUE) 
# 
# Map_pG(pGvec,2050)
# dev.off()

# pG maps functions
##########################################################################################################################################
##########################################################################################################################################
# map_from_data <-function(variable.df, model.choice, run.choice, varname, savemat=FALSE){
#   
#   stringYear<-toString(model.choice)
#   stringYear<-toString(run.choice)
#   
#   Z_historical <- variable.df %>% 
#     filter(experiment == "historical" & model == model.choice & run == run.choice & between(year, 1850, 2014))  %>% 
#     arrange(year) %>%
#     select(!c(institute,model,experiment,run, year))
#   Z_historical<-as.data.frame(Z_historical)
#   Z_rcp85<- variable.df %>% 
#     filter(experiment == "rcp85" & model == model.choice & run == run.choice & between(year, 2015, 2100))  %>%
#     arrange(year) %>% 
#     select(!c(institute,model,experiment,run, year))
#   Z_rcp85<-as.data.frame(Z_rcp85)
#   matz <- bind_rows(Z_historical,Z_rcp85)
#   matz <-as.matrix(matz)
#   
#   matx<- variable.df %>%
#     filter(experiment == "historicalNat" & model == model.choice & run == run.choice)  %>%
#     arrange(year) %>% select(!c(institute,model,experiment,run, year))
#   matx<-as.matrix(matx)
#   
#   tt<-c(1:dim(matz)[1])
# 
#   if (savemat==TRUE){
#     matx <- as.data.frame(matx)
#     matz <- as.data.frame(matz)
#     matx_file_name <- paste("matx_",varname,model.choice,"_",run.choice,".rds",sep ="")
#     matz_file_name <- paste("matz_",varname,model.choice,"_",run.choice,".rds",sep ="")
#     saveRDS(matx, file=matx_file_name)
#     saveRDS(matz, file=matz_file_name)
#   }
#   
#   return(list("matX"=matx,"matZ"=matz,"time.vec"=tt))
# }
# 
# 
# Map_pG <- function(pG.vec,year){
#   
#   gridmat_pG <- matrix(data=pG.vec,nrow=72,ncol=36)
#   grid_lon<-seq(-177.5, 177.5, length.out = 72) 
#   grid_lat<-seq(-87.5, 87.5, length.out = 36)
#   rfcol <- colorRampPalette(brewer.pal(9,'Reds'))
#   colsp <- rfcol(64)
#   r_breaks <- seq(0,1,length.out=65)
#   image.plot(grid_lon,grid_lat,gridmat_pG,col = colsp,breaks = r_breaks,
#              legend.args = list( text = expression(p["G"]),
#                                  cex = 1,
#                                  side = 3,
#                                  line = .5),ylab="latitude",xlab="longitude")
#   stringYear<-toString(year)
#   title(stringYear,cex=0.5,line = .8)
#   world(lwd = 1,add=TRUE) 
# }
