#install.packages("datasets.load")
library(datasets)
View(Nile) #view the NIle data
plot(Nile)

#exports=read.table(file = 'year_origin_destination_sitc_rev2.tsv', sep = '\t', header = TRUE)
#View(exports)

require(KFAS)
# Example of local level model for Nile series
# See Durbin and Koopman (2012)

model_Nile <- SSModel(Nile ~
                        SSMtrend(1, Q = list(matrix(NA))), H = matrix(NA))
model_Nile
model_Nile <- fitSSM(model_Nile, c(log(var(Nile)), log(var(Nile))),
                     method = "BFGS",hessian=TRUE)

#hessian matrix from LL fitting
model_Nile$optim.out$hessian

# Filtering and state smoothing
out_Nile <- KFS(model_Nile$model, filtering = "state", smoothing = c("state","disturbance"),simplify=F)
out_Nile
ts.plot(cbind(Nile, out_Nile$a[-1], out_Nile$att[-1],out_Nile$alphahat[-1]),col = c(1:4),
        ylab = "filtered and smoothed states", main = "River Nile")

# plot to compare with matlab outputs
plot(as.ts(out_Nile$a[-1]))
plot(as.ts(t(out_Nile$r)))
plot(as.ts(t(out_Nile$N)))
plot(as.ts(out_Nile$etahat))
plot(as.ts(out_Nile$epshat))
plot(as.ts(rstandard(out_Nile,type="state")))
plot(as.ts(rstandard(out_Nile,type="recursive")))


## shock testing
dt=matrix(NA,length(Nile),2)
dt2=matrix(NA,length(Nile),2)
for (tt in 1:length(Nile)){
  shock=rep(0,length(Nile))
  shock[tt]=1
  
  shock2=rep(0,length(Nile))
  shock2[tt:length(Nile)]=1
  
  ## dummy in observation equation
  model_Nile_shock <- SSModel(Nile ~ shock +
                          SSMtrend(1, Q = list(matrix(NA))), H = matrix(NA))
  model_Nile_shock
  model_Nile_shock <- fitSSM(model_Nile_shock, c(log(var(Nile)), log(var(Nile))),
                       method = "BFGS",hessian=TRUE)
  
  #hessian matrix from LL fitting
  model_Nile_shock$optim.out$hessian
  
  # Filtering and state smoothing
  out_Nil_shock<- KFS(model_Nile_shock$model, filtering = "state", smoothing = c("state","disturbance"),simplify=F)
  
  #calculation of d(t,j)
  dt[tt,1]=(log(model_Nile_shock$model['Q'][1,1,])-log(model_Nile$model['Q'][1,1,]))/sqrt(model_Nile$optim.out$hessian[1,1])
  dt[tt,2]=(log(model_Nile_shock$model['H'][1,1,])-log(model_Nile$model['H'][1,1,]))/sqrt(model_Nile$optim.out$hessian[2,2])
  
  
  #################################################################
  ## Second type of schock
  ##############################################################
  
  ## dummy in observation equation
  model_Nile_shock2 <- SSModel(Nile ~ shock2 +
                                SSMtrend(1, Q = list(matrix(NA))), H = matrix(NA))
  model_Nile_shock2
  model_Nile_shock2 <- fitSSM(model_Nile_shock2, c(log(var(Nile)), log(var(Nile))),
                             method = "BFGS",hessian=TRUE)
  
  #hessian matrix from LL fitting
  model_Nile_shock2$optim.out$hessian
  
  # Filtering and state smoothing
  out_Nil_shock2<- KFS(model_Nile_shock2$model, filtering = "state", smoothing = c("state","disturbance"),simplify=F)
  
  #calculation of d(t,j)
  dt2[tt,1]=(log(model_Nile_shock2$model['Q'][1,1,])-log(model_Nile$model['Q'][1,1,]))/sqrt(model_Nile$optim.out$hessian[1,1])
  dt2[tt,2]=(log(model_Nile_shock2$model['H'][1,1,])-log(model_Nile$model['H'][1,1,]))/sqrt(model_Nile$optim.out$hessian[2,2])
  
}  

plot(dt[,1],type='h')
plot(dt[,2],type='h')

###############################################################
# Same but with a simulated data
###############################################################
T=100;  # sample size
d = 0;      # intercept parameter in update equation
c = 0;     # autoregressive parameter in update equation
sigma_eta = 1469.1; # standard error of innovations in update equation Q
sigma_eps = 15099;  # standard error of innovations in observation equation H
mu1 = 1120; # define initial value for time series x

#generate innovations
eta =  sigma_eta*rnorm(T); #innovations for observations
epsilon = sigma_eps*rnorm(T);#innovations for observations

#define time series vector
mu = rep(0,T);
y = rep(0,T);

#initialization
mu[1] = mu1;

#generate the observations and states
for (t in 1:T){
  mu[t+1] = d + mu[t] + eta[t];
  y[t] = c + mu[t] + epsilon[t];
}

plot(as.ts(y))

###SIMULATE SHOCKS  
y[15] = y[43] - 0.3 * y[43];
#y[44:80] = y[44:80] * 1.3;

plot(as.ts(y))

######################################################
###*** Apply our algorithm ***########################
#####################################################

#firstly fit the basic model
model<- SSModel(y ~ SSMtrend(1, Q = list(matrix(NA))), H = matrix(NA))
model <- fitSSM(model, c(0, 0),
                     method = "BFGS",hessian=TRUE)

#hessian matrix from LL fitting
model$optim.out$hessian

# Filtering and state smoothing
out <- KFS(model$model, filtering = "state", smoothing = c("state","disturbance"),simplify=F)
out

## shock testing
dt=matrix(NA,length(y),2)
dt2=matrix(NA,length(y),2)
for (tt in 1:length(y)){
  shock=rep(0,length(y))
  shock[tt]=1
  
  shock2=rep(0,length(y))
  shock2[tt:length(y)]=1
  
  ## dummy in observation equation
  model_shock <- SSModel(y ~ Xt +
                                SSMtrend(1, Q = list(matrix(NA))), H = matrix(NA))
  model_shock
  model_shock <- fitSSM(model_shock, c(0, 0),
                             method = "L-BFGS-B",hessian=TRUE)
  
  #hessian matrix from LL fitting
  model_shock$optim.out$hessian
  
  # Filtering and state smoothing
  out_shock<- KFS(model_shock, filtering = "state", smoothing = c("state","disturbance"),simplify=T)
  
  #calculation of d(t,j)
  dt[tt,1]=(log(model_shock$model['Q'][1,1,])-log(model$model['Q'][1,1,]))/sqrt(1)
  dt[tt,2]=(log(model_shock$model['H'][1,1,])-log(model$model['H'][1,1,]))/sqrt(1)
  
}  

plot(dt[,1],type='h')
plot(dt[,2],type='h')

###********############################################################
# using package "stsm" to get second derivates and scores of likelihood function
###********####################################################
library(stsm)
#initialisation of variance
initpars <- c(var1 = 1, var2 = 1,P0=10000000,a0=0)
mairp <- stsm.model(model = "local-level",y = Nile, 
                pars = initpars, nopars = NULL, transPars = NULL)

# procedure 0: ML-TD StructTS()
res0 <- StructTS(Nile, type = "level")$coef

#mairp <- stsm.model(model = "local-level", y = Nile,
#                  transPars = "StructTS")

res1 <- maxlik.td.optim(mairp, KF.version = "KFKSDS",
                          KF.args = list(P0cov = TRUE), method = "L-BFGS-B", gr = "analytical")

mairpNile <- set.pars(mairp, pmax(res3$par, .Machine$double.eps))
round(get.pars(mairpNile), 6)
all.equal(get.pars(mairpNile), res0, tol = 1e-05, check.att = FALSE)

#analytical hessian computed via Newton-Raphson
res2 <- maxlik.fd.scoring(mairp, step = NULL,
                    information = "observed", 
                    control = list(maxit = 10000, tol = 0.001),
                    ls = list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1))

#this one gives the right values for parameters
res3 <- maxlik.td.scoring(mairp, 
                          KF.args = list(P0cov = FALSE), check.KF.args = TRUE,
                          ls = list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1),
                          control = list(maxit = 10000, tol = 0.001, trace = TRUE))
res2$Dmat  
res3$infomat

# Scoring algorithm (information matrix)
res4 <- maxlik.fd.scoring(m = m, step = NULL,
                          information = "expected", 
                          control = list(maxit = 10000, tol = 0.001),
                          ls = list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1))
res4$Dmat

print(res2)

plot(tsSmooth(StructTS(mairpNile@y, type = "level"))) #smoothed state
ssNile <- char2numeric(mairpNile) #produces parameters and properties of our ss model
kfNile <-KFKSDS:::KF(mairpNile@y, ssNile) #Kalman filtering
ksNile <- KFKSDS:::KS(mairpNile@y, ssNile, kfNile) #Kalman smoothing
KFKSDS 
plot(ksNile$ahat)

ts.plot(cbind(Nile, kfNile$a.upd,kfNile$a.pred,ksNile$ahat),col = c(1:4),
        ylab = "filtered and smoothed states", main = "River Nile")

#compare with KFAS output
ts.plot(cbind(Nile,out_Nile$att,out_Nile$a[-1],out_Nile$alphahat),col = c(1:4),
        ylab = "filtered and smoothed states", main = "River Nile")

#compute derivates of log-likelihood
kfNile.deriv = KFKSDS:::KF.deriv(mairpNile@y, ssNile)

#### Experiment with Nile data by stsm package (replication of Atkinson et.al)
out = matrix(0,length(Nile),2) #output matrix for two types
for (t in seq((1:length(Nile)))){
  
  #simulate additive schock at time t
  Xt = rep(0,length(y));
  Xt[t] = 10; # put the shock on the first time index t=1
  
  #call auxiliary KFS
  out[t,1] =kafs(Nile, ssNile, Xt)$d[1]
  out[t,2] =kafs(Nile, ssNile, Xt)$d[2] 
}

#plot scaled distances
plot(out[,1],type='h')
plot(out[,2],type='h')

###############################################################
# Simulate data with use of package tsoutliers
###############################################################
T=100;  # sample size
d = 0;      # intercept parameter in update equation
c = 0;     # autoregressive parameter in update equation
sigma_eta = 1469.163; # standard error of innovations in update equation Q
sigma_eps = 15098.65;  # standard error of innovations in observation equation H
mu1 = 1120; # define initial value for time series x

#generate innovations
eta =  sigma_eta*rnorm(T); #innovations for observations
epsilon = sigma_eps*rnorm(T);#innovations for observations

#define time series vector
mu = rep(0,T);
y = rep(0,T);

#initialization
mu[1] = mu1;

#generate the observations and states
for (t in 1:T){
  mu[t+1] = d + mu[t] + eta[t];
  y[t] = c + mu[t] + epsilon[t];
}

plot(as.ts(y))

#produce outliers (2 AO, 1 LS)
y[15] = -70000;
y[45] = 47002.25;
y[70:90] = y[70:90] + 60000.45;

library(tsoutliers)
y = as.ts(y)
res0 <- StructTS(y, type = "level")$coef
mairp <- stsm.model(model = "local-level", y = y,
                    transPars = "StructTS")
res1 <- maxlik.td.optim(mairp, KF.version = "KFKSDS",
                        KF.args = list(P0cov = TRUE), method = "L-BFGS-B", gr = "analytical")
mair<- set.pars(mairp, pmax(res1$par, .Machine$double.eps))

plot(tsSmooth(StructTS(mairp@y, type = "level"))) #smoothed state
ss <- char2numeric(mairp) #produces parameters and properties of our ss model
kf <-KFKSDS:::KF(mairp@y, ss) #Kalman filtering
ks <- KFKSDS:::KS(mairp@y, ss, kf) #Kalman smoothing
ds <- KFKSDS:::DS(mairp@y, ss, kf,ks) #Disturbance smoothing

resid = ds$epshat #additive shocks
resid2 = ds$etahat #level shifts
n <- length(resid)

#standart deviations of residuals
sigma = ds$vareps
sigma2 = ds$vareta

#detection of outliers via outer and inner loop
stage1 <- locate.outliers.oloop(y, res1, types = c("IO", "AO", "LS", "TC"))
stage1$outliers

#*********************************************#
# TEST FOR i simulations
#*********************************************#
for (i in seq_len(dim(llm)[2]))
  {
  y=llm[,i]
  y[15] = min(y)-mean(y);
  y[45] = max(y)+mean(y);
  y[70:90] = y[70:90] + min(y)+mean(y)-y[70:90]*0.12;
  mairp <- stsm.model(model = "local-level", y = y,
                      transPars = "StructTS")
  res1 <- maxlik.td.optim(mairp, KF.version = "KFKSDS",
                          KF.args = list(P0cov = TRUE), method = "L-BFGS-B", gr = "analytical")
  stage1 <- locate.outliers.oloop(y, res1, types = c("IO", "AO", "LS", "TC"))
  print(stage1$outliers)
}

#### same step by all-automated function
#EXPERIMENT I.
res=rep(0,dim(llm)[2])
fit=matrix(0,dim(llm)[2],2)
colnames(fit)=c('var1','var2')
for (i in seq_len(dim(llm)[2]))
  {
  y=llm[,i]
  y[15] = min(y)-mean(y)
  y[45] = max(y)+2*mean(y)
  y[70:90] = y[70:90] + min(y)+1.2*mean(y)-y[70:90]*0.12
  res[i] <- tso(y = y, types = c("AO", "LS", "TC"),
                tsmethod = "stsm", args.tsmodel = list(model = "local-level",
                                                       pars = c("var1" = log(var(y)), "var2" =  log(var(y)))
                                                       )
                )
  #plot(res[i])
  plot(tso(y = y, types = c("AO", "LS", "TC"),
            tsmethod = "stsm", args.tsmodel = list(model = "local-level",
                                                   pars = c("var1" = log(var(y)), "var2" =  log(var(y))))
  ))
  fit[i,] = tso(y = y, types = c("AO", "LS", "TC"),
               tsmethod = "stsm", args.tsmodel = list(model = "local-level"))$fit$model@pars[1:2]
}

print(res)
print(fit)

#for (i in seq_len(dim(llm)[2])){
#  plot(res[[i]])
#}

#EXPERIMENT II.
# generate a quarterly series from a local level plus seasonal model
require(stsm)
pars <- c(var1 = 15098.65, var2 = 1469.163,a01 = 1120)
m <- stsm.model(model = "local-level", y = ts(seq(100)), 
                pars = pars, nopars = NULL)
ss <- char2numeric(m)

iter = 10 #set the number of simulations
y = matrix(0,T,iter)
mu = matrix(0,T,iter)

res=rep(0,dim(llm)[2])
fit=matrix(0,dim(llm)[2],2)
colnames(fit)=c('var1','var2')

set.seed(123)
for (i in seq(iter)){
  y[,i] <- datagen.stsm(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                  n0 = 0, old.version = TRUE)$data
  mu[,i] <- datagen.stsm(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                  n0 = 0, old.version = TRUE)$components

  mu[15,i] =max(mu[15,i]) + mean(mu[15,i])

  #reflect the shock in observations
  eps1 <- rnorm(T, sd = sqrt(ss$H[1]))

  for (j in seq(T)) {
    y[j,i] <- ss$Z %*% mu[j,i] + eps1[j]
  }
  y[,i] <- ts(y[,i], frequency =1)
  mu <- ts(mu, frequency = 1)
  
  y_m=as.ts(y[,i])
  res[i] <- tso(y = y_m, types = c("AO", "LS", "TC"),
                tsmethod = "stsm", args.tsmodel = list(model = "local-level",
                                                       pars = c("var1" = log(var(y_m)), "var2" =  log(var(y_m)))
                )
  )
  plot(tso(y = y_m, types = c("AO", "LS", "TC"),
           tsmethod = "stsm", args.tsmodel = list(model = "local-level",
                                                  pars = c("var1" = log(var(y_m)), "var2" =  log(var(y_m))))
  ))
  fit[i,] = tso(y = y_m, types = c("AO", "LS", "TC"),
                tsmethod = "stsm", args.tsmodel = list(model = "local-level"))$fit$model@pars[1:2]
}
print(res)
print(fit) 

#############################################################################
#### Experiment with simulated data by stsm package (replication of Atkinson et.al)
########

#simulate your LLM
pars <- c(var1 = 15098.65, var2 = 1469.163,a01 = 1120)
m <- stsm.model(model = "local-level", y = ts(seq(100)), 
                pars = pars, nopars = NULL)
ss <- char2numeric(m)

iter = 10 #set the number of simulations
y = matrix(0,T,iter)
mu = matrix(0,T,iter)

res=rep(0,length(iter))
fit=matrix(0,length(iter),2)
colnames(fit)=c('var1','var2')

set.seed(123)
for (i in seq(iter)){
  y[,i] <- datagen.stsm(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                        n0 = 0, old.version = TRUE)$data
  #simulation of two AO and one LS
  y_m=as.ts(y[,i])
  y_m[15] = min(y_m)-0.3*mean(y_m)
  y_m[45] = max(y_m)+0.2*mean(y_m)
  y_m[70:90] = y_m[70:90] + min(y_m)+0.15*mean(y_m)-y_m[70:90]*0.12
  results=list()
  out = matrix(0,length(y_m),2) #output matrix for two types
  #out2 = matrix(0,length(y_m),2) #output matrix for two types
    for (t in seq((1:length(y_m)))){
      #simulate additive schock at time t
      Xt = rep(0,length(y_m));
      Xt[t] = 10; # put the shock on the first time index t=1 AO
      #Xt[t:(t+10)] =Xt[t:(t+10)] + 10; # put the shock on the first time index t=1 LS
  
      #call auxiliary KFS
      results[[t]] =kafs(y_m, Xt)
      #out2[t,1] =kafs(y_m, Xt)$dist.v2[1]
      #out[t,1] =results$dist[1]
      #out[t,2] =results$dist[2]
}

    #plot scaled distances
    plot(out[,1],type='h')
    plot(out[,2],type='h')
    #plot(out2[,1],type='h')
   # plot(out2[,2],type='h')
}

#############################################################################
#### Experiment with simulated data by gets package (indicator saturation approach)
########

library(gets)

#simulate your LLM
pars <- c(var1 = 15098.65, var2 = 1469.163,a01 = 1120)
m <- stsm.model(model = "local-level", y = ts(seq(100)), 
                pars = pars, nopars = NULL)
ss <- char2numeric(m)

T=100 #number of time series
iter = 10 #set the number of simulations
y = matrix(0,T,iter)
mu = matrix(0,T,iter)

res=rep(0,length(iter))
out = list() 

set.seed(123)
for (i in seq(iter)){
  y[,i] <- datagen.stsm(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                        n0 = 0, old.version = TRUE)$data
 
  #simulation of two AO and one LS
  y_m=as.ts(y[,i])
  y_m[15] = min(y_m) - 0.3 * mean(y_m)
  y_m[45] = max(y_m) + 0.2 * mean(y_m)
  y_m[70:90] = y_m[70:90] + min(y_m)+0.15*mean(y_m)-y_m[70:90]*0.12
  
  #measure computation time
  start_time <- Sys.time()
  
  #output list
  out[[i]]=isat(y_m,iis = TRUE, sis = TRUE, tis = FALSE,t.pval=1/T)
  plot(out[[i]])
  end_time <- Sys.time()
  print(end_time - start_time)
}

outliers = c(15,45,70:90)

#performace of indicator saturation $ i dont know about this one, you have to figure out
#retention rate
library(stringr)
retention=rep(0,iter)
retention.out=rep(0,iter)
for (t in seq(iter)){
  indexes=rep(0,length(out[[1]]$ISnames))
  indexes.sis1=rep(0,length(out[[1]]$ISnames))
  indexes.sis2=rep(0,length(out[[1]]$ISnames))
  for (j in 1:length(out[[1]]$ISnames)){
    indexes[j]= as.numeric(unlist(strsplit(out[[t]]$ISnames[j], split='iis', fixed=TRUE))[2])
    indexes.sis1[j] = as.numeric(unlist(strsplit(out[[t]]$ISnames[j], split='sis', fixed=TRUE))[2])[c(TRUE,FALSE)]
    indexes.sis2[j] = as.numeric(unlist(strsplit(out[[t]]$ISnames[j], split='sis', fixed=TRUE))[2])[c(FALSE,TRUE)]
    #y_m[indexes[j],t]#have to figure out this part
  }
  retention[t]=length(intersect(outliers,indexes))
  retention.out[t]=length(setdiff(outliers,indexes))
}

#compute retention rate
retention_rate=mean(retention)
retention_rate.out=mean(retention.out)

#recal/potency
potency = retention_rate/length(outliers)
#gauge
gauge = (1/(T-length(outliers)))*retention_rate
#########******************************************############################
#### Experiment with simulated data by gets package (indicator saturation approach)
#### My version ####################################################################
##############################################################v###################
library(gets)

#simulate your LLM
T = 100
pars <- c(var1 = 15099, var2 = 1469.1,a01 = 1120)
m <- stsm.model(model = "local-level", y = ts(seq(100)), 
                pars = pars, nopars = NULL)
ss <- char2numeric(m)

iter = 10 #set the number of simulations
y = matrix(0,T,iter)
mu = matrix(0,T,iter)

res=rep(0,length(iter))
out3 = list() 

set.seed(123)
for (i in seq(iter)){
  y[,i] <- datagen.stsm(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                        n0 = 0, old.version = TRUE)$data
  
  #simulation of two AO and one LS
  y_m=as.ts(y[1:100,i])
  y_m[15] = min(y_m) - 0.3 * mean(y_m)
  y_m[45] = max(y_m) + 0.2 * mean(y_m)
  y_m[70:90] = y_m[70:90] + min(y_m)+0.15*mean(y_m)-y_m[70:90]*0.12
  
  #measure computation time
  start_time <- Sys.time()
  
  #output list, this guy counts cca 4-5mins for each time serie in TURBO version, 
  #probably would be desirable to put parallelision on or test without hessian estimation
  #(is i wont implement the version of atkinson et al)
  out3[[i]]=isat_my(y_m, t.pval=1/T, model.ssmod = "local-level",turbo=TRUE,
                   vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                   optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
  end_time <- Sys.time()
  print(end_time - start_time)
  plot(out3[[i]])
}

isat_my(y_m, t.pval=0.001, model.ssmod = "local-level",
                 vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                 optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))

#########******************************************############################
####incorporation with seasonality term ####################################################################
##############################################################v###################
#simulate your LLM
T = 100
#variance setted to 1 otherwise it is very slow to estimate coefficients, although if we have a time varying seasonality in data then it would be probably
#beneficial to do so
pars <- c(var1 = 15099, var2 = 1469.1 ,var3 = 1, a01 = 1120)

#for now only a weekly frequency set
m <- stsm.model(model = "llm+seas", y = ts(seq(104),start=c(2016,01),frequency = 52), #here a multinomial frequency in ts would be nice to have
                pars = pars, nopars = NULL)
ss <- char2numeric(m)

iter = 10 #set the number of simulations
y = matrix(0,T,iter)
#mu = matrix(0,T,iter)

res=rep(0,length(iter))
out.seas = list() 

set.seed(123)
for (i in seq(iter)){
  y[,i] <- datagen.stsm(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                        n0 = 0, old.version = TRUE)$data
  
  #simulation of two AO and one LS
  y_m=ts(y[,i],start=c(2016,01),frequency = 52) #
  y_m[15] = min(y_m) - 0.3 * mean(y_m)
  y_m[45] = max(y_m) + 0.2 * mean(y_m)
  y_m[70:90] = y_m[70:90] + min(y_m)+0.15*mean(y_m)-y_m[70:90]*0.12
  
  #measure computation time
  start_time <- Sys.time()
  
  #output list, this guy counts cca 4-5mins for each time serie in TURBO version, 
  #probably would be desirable to put parallelision on or test without hessian estimation
  #(is i wont implement the version of atkinson et al)
  out.seas[[i]]=isat_my(y_m, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                    vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                    optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
  end_time <- Sys.time()
  print(end_time - start_time)
  out.seas.adj=out.seas[[i]]
  out.seas.adj$mean.fit = out.seas.adj$mean.fit[,"level"]
  #out.seas.adj$aux$mX=coredata(out.seas.adj$aux$mX)
  plot.isat(out.seas.adj)
}

#########******************************************######################
####BASIC STRUCTURAL MODEL ##############################################
##############################################################v###################
#simulate your LLM
T =7*10
T=100
#T=200
#variance setted to 1 otherwise it is very slow to estimate coefficients, although if we have a time varying seasonality in data then it would be probably
#beneficial to do so
#pars <- c(var1 = 15099, var2 = 1469.1 ,var3 = 5, var4 = 1, a01 = 1120)
pars <- c(var1 = 10, var2 = 6 ,var3 = 0.0001, var4 = 8, a01 = 25)

#scenario sTsS
pars <- c(var1 = 10, var2 = 6 ,var3 = 0.00008, var4 =  0.00005, a01 = 25)

#scenario uTsS
pars <- c(var1 = 10, var2 = 6 ,var3 = 0.8, var4 =  0.00005, a01 = 25)

#scenario sTuS
pars <- c(var1 = 10, var2 = 6 ,var3 = 0.00008, var4 = 0.5, a01 = 25)

#scenario uTuS
pars <- c(var1 = 10, var2 = 6 ,var3 = 0.8, var4 =  0.5, a01 = 25)

#for now only a weekly frequency set
m <- stsm::stsm.model(model = "BSM", y = ts(seq(T),freq=7), #here a multinomial frequency in ts would be nice to have
                pars = pars, nopars = NULL)
ss <- char2numeric(m)

iter = 10 #set the number of simulations
y = matrix(0,T,iter)
#mu = matrix(0,T,iter)

res=rep(0,length(iter))
out.bsc = list() 
out.bsc.eps= list()
out.bsc.eta= list()
out.bsc.zeta= list()
out.bsc.omega= list()
#states = array(0, dim = c(T,3, iter))

set.seed(123)
for (i in seq(iter)){
  #this one is the old simulation setting
 # y[,i] <- datagen.stsm(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
  #                      n0 = 70, old.version = TRUE)$data
 # y[,i]= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
 #                n0 = 70, old.version = TRUE,
 #                 AOplus=c(15),AOplusmag=6, AOminus=c(40),AOminusmag=10,
 #                 LSplus=c(50),LSplusmag=8,
 #                TCminus=c(10), TCminusmag=8,SAOminus=c(70), SAOminusmag=7
 #)$data
  y[,i]= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                         n0 = 70, old.version = TRUE,
                         AOplus=c(50),AOplusmag=10, AOminus=c(40),AOminusmag=0,
                         LSplus=c(50),LSplusmag=0,
                         TCminus=c(10), TCminusmag=0,SAOminus=c(70), SAOminusmag=0
  )$data
  #setting as in marczak
 #y[,i]= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
 #                       n0 = 70, old.version = TRUE,
   #                     AOplus=c(15),AOplusmag=7, AOminus=c(40),AOminusmag=7,
   #                     LSplus=c(50),LSplusmag=7,
    #                     TCminus=c(10), TCminusmag=7,SAOminus=c(70), SAOminusmag=7
 # )$data
  #simulation of two AO and one LS
  y_m=ts(y[,i],frequency =7) #
 # y_m[15] = min(y_m) - 0.3 * mean(y_m)
 # y_m[45] = max(y_m) + 0.25 * mean(y_m)
 # y_m[70:90] = y_m[70:90] + min(y_m)+0.15*mean(y_m)-y_m[70:90]*0.12
  
  #measure computation time
  start_time <- Sys.time()
  print(i) #show the iteration number
  #output list, this guy counts cca 4-5mins for each time serie in TURBO version, 
  #probably would be desirable to put parallelision on or test without hessian estimation
  #(is i wont implement the version of atkinson et al)
  out.bsc[[i]]=isat_my(y_m, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                        vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                        optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
  plot(out.bsc[[i]])
  out.bsc.eps[[i]]=isat(out.bsc[[i]]$pure$std.eps,sis=T,iis=F)
  plot(out.bsc.eps[[i]])
  out.bsc.eta[[i]]=isat(out.bsc[[i]]$pure$std.eta,sis=T,iis=F)#for the level
  plot(out.bsc.eta[[i]])
  out.bsc.zeta[[i]]=isat(out.bsc[[i]]$pure$std.zeta,sis=T,iis=F)#for the slope
  plot(out.bsc.zeta[[i]])
  out.bsc.omega[[i]]=isat(out.bsc[[i]]$pure$std.omega,sis=T,iis=F)#for seasonality
  plot(out.bsc.omega[[i]])
  end_time <- Sys.time()
  print(end_time - start_time)
 # out.bsc.adj= out.bsc[[i]]
  #out.bsc.adj$mean.fit = out.bsc.adj$mean.fit[,"level"]
  #out.seas.adj$aux$mX=coredata(out.seas.adj$aux$mX)
 # plot.isat(out.bsc.adj)
}
for (i in seq(iter)){
  # print(i)
   string=out.bsc.eps[[i]]$ISnames #3  #6.571429e+14
  # string=out.bsc.eta[[i]]$ISnames
  # string=out.bsc.zeta[[i]]$ISnames
   #  string=out.bsc.omega[[i]]$ISnames
print(as.numeric(gsub("\\D", "",string)) )
}
confusion.iis(out.bsc.eps,out.bsc.eta,out.bsc.zeta,out.bsc.omega,vec.outlier1)

plot(as.ts(y[,1],frequency=7),main='Series with two additive outliers (t=15, t=45) and one level shift (t=50)'
     ,ylab='y')
plot(as.ts(y[,2],frequency=7),main='Series with two additive outliers (t=15, t=45) and one level shift (t=50)'
     ,ylab='y')

par(mfrow = c(1,1),mai = c(1, 0.42, 0.1,0.1))
plot(as.ts(pokus,frequency=7),main=''
     ,ylab='',cex.axis = 0.5,bty="n",frame.plot=FALSE)
plot(as.ts(pokus1,frequency=7),main=''
     ,ylab='', cex.axis = 0.5,bty="n",frame.plot=FALSE)


################################################################
# AO scenarios ##
T=100

#SCENARIO 1
pars <- c(var1 = 10, var2 = 6 ,var3 =  8*0.0001, var4 =5*0.0001, a01 = 25)

#for now only a weekly frequency set
m <- stsm::stsm.model(model = "BSM", y = ts(seq(T),freq=7), #here a multinomial frequency in ts would be nice to have
                      pars = pars, nopars = NULL)
ss <- char2numeric(m)

#t=10
series1a= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                n0 = 70, old.version = TRUE,
                AOplus=c(10),AOplusmag=6)$data

#t=80
series1b= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                         n0 = 70, old.version = TRUE,
                         AOplus=c(80),AOplusmag=6)$data

#plot series
par(mfrow = c(1,1),mai = c(1, 0.42, 0.1,0.1))
plot(as.ts(series1a,frequency=7),main=''
     ,ylab='',cex.axis = 0.5,bty="n",frame.plot=FALSE)
plot(as.ts(series1b,frequency=7),main=''
     ,ylab='', cex.axis = 0.5,bty="n",frame.plot=FALSE)

#SCENARIO 2
pars <- c(var1 = 10, var2 = 6 ,var3 = 0.8, var4 =5*0.0001, a01 = 25)

#for now only a weekly frequency set
m <- stsm::stsm.model(model = "BSM", y = ts(seq(T),freq=7), #here a multinomial frequency in ts would be nice to have
                      pars = pars, nopars = NULL)
ss <- char2numeric(m)

#t=10
series1a= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          AOplus=c(10),AOplusmag=6)$data

#t=80
series1b= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          AOplus=c(80),AOplusmag=6)$data

#plot series
par(mfrow = c(1,1),mai = c(1, 0.42, 0.1,0.1))
plot(as.ts(series1a,frequency=7),main=''
     ,ylab='',cex.axis = 0.5,bty="n",frame.plot=FALSE)
plot(as.ts(series1b,frequency=7),main=''
     ,ylab='', cex.axis = 0.5,bty="n",frame.plot=FALSE)

#SCENARIO 3
pars <- c(var1 = 10, var2 = 6 ,var3 = 8*0.0001, var4 =0.5, a01 = 25)

#for now only a weekly frequency set
m <- stsm::stsm.model(model = "BSM", y = ts(seq(T),freq=7), #here a multinomial frequency in ts would be nice to have
                      pars = pars, nopars = NULL)
ss <- char2numeric(m)

#t=10
series1a= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          AOplus=c(10),AOplusmag=6)$data

#t=80
series1b= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          AOplus=c(80),AOplusmag=6)$data

#plot series
par(mfrow = c(1,1),mai = c(1, 0.42, 0.1,0.1))
plot(as.ts(series1a,frequency=7),main=''
     ,ylab='',cex.axis = 0.5,bty="n",frame.plot=FALSE)
plot(as.ts(series1b,frequency=7),main=''
     ,ylab='', cex.axis = 0.5,bty="n",frame.plot=FALSE)

#SCENARIO 4
pars <- c(var1 = 10, var2 = 6 ,var3 = 0.8, var4 =0.5, a01 = 25)

#for now only a weekly frequency set
m <- stsm::stsm.model(model = "BSM", y = ts(seq(T),freq=7), #here a multinomial frequency in ts would be nice to have
                      pars = pars, nopars = NULL)
ss <- char2numeric(m)

#t=10
series1a= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          AOplus=c(10),AOplusmag=6)$data

#t=80
series1b= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          AOplus=c(80),AOplusmag=6)$data

#plot series
par(mfrow = c(1,1),mai = c(1, 0.42, 0.1,0.1))
plot(as.ts(series1a,frequency=7),main=''
     ,ylab='',cex.axis = 0.5,bty="n",frame.plot=FALSE)
plot(as.ts(series1b,frequency=7),main=''
     ,ylab='', cex.axis = 0.5,bty="n",frame.plot=FALSE)

#BENCHMARK SCENARIO
pars <- c(var1 = 10, var2 = 6 ,var3 = 0.0001, var4 =8, a01 = 25)

#for now only a weekly frequency set
m <- stsm::stsm.model(model = "BSM", y = ts(seq(T),freq=7), #here a multinomial frequency in ts would be nice to have
                      pars = pars, nopars = NULL)
ss <- char2numeric(m)

#t=10
series1a= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          AOplus=c(15),AOplusmag=6, AOminus=c(40),AOminusmag=10,
                          LSplus=c(50),LSplusmag=8,
                          TCminus=c(10), TCminusmag=8,SAOminus=c(70), SAOminusmag=7
                          )$data

# SCENARIO POSITION 1
series1b= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          AOplus=c(10),AOplusmag=6,
                          LSminus=c(80),LSminusmag=8)$data
# SCENARIO POSITION 2
series1c= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          AOminus=c(80),AOminusmag=10,
                          LSplus=c(10),LSplusmag=8)$data
# SCENARIO POSITION 3
series1d= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          AOminus=c(10),AOminusmag=10,
                          LSplus=c(80),LSplusmag=8)$data
# SCENARIO POSITION 4
series1e= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          AOminus=c(80),AOminusmag=10,AOplus=c(10),AOplusmag=6,
                          LSminus=c(80),LSminusmag=8)$data

#plot series
par(mfrow = c(1,1),mai = c(1, 0.42, 0.1,0.1))
plot(as.ts(series1a,frequency=7),main=''
     ,ylab='',cex.axis = 0.5,bty="n",frame.plot=FALSE)
plot(as.ts(series1b,frequency=7),main=''
     ,ylab='', cex.axis = 0.5,bty="n",frame.plot=FALSE)
plot(as.ts(series1c,frequency=7),main=''
     ,ylab='',cex.axis = 0.5,bty="n",frame.plot=FALSE)
plot(as.ts(series1d,frequency=7),main=''
     ,ylab='', cex.axis = 0.5,bty="n",frame.plot=FALSE)
plot(as.ts(series1e,frequency=7),main=''
     ,ylab='', cex.axis = 0.5,bty="n",frame.plot=FALSE)

# DIFFERENT SIZES
#LS
series1a= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          LSplus=c(50),LSplusmag=2)$data
series1b= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          LSplus=c(50),LSplusmag=20)$data

#plot series
par(mfrow = c(1,1),mai = c(1, 0.42, 0.1,0.1))
plot(as.ts(series1a,frequency=7),main=''
     ,ylab='',cex.axis = 0.5,bty="n",frame.plot=FALSE)
plot(as.ts(series1b,frequency=7),main=''
     ,ylab='', cex.axis = 0.5,bty="n",frame.plot=FALSE)

#TC
series1a= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          TCplus=c(50),TCplusmag=2)$data
series1b= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          TCplus=c(50),TCplusmag=20)$data

#plot series
par(mfrow = c(1,1),mai = c(1, 0.42, 0.1,0.1))
plot(as.ts(series1a,frequency=7),main=''
     ,ylab='',cex.axis = 0.5,bty="n",frame.plot=FALSE)
plot(as.ts(series1b,frequency=7),main=''
     ,ylab='', cex.axis = 0.5,bty="n",frame.plot=FALSE)

#SC
series1a= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          SAOplus=c(50),SAOplusmag=2)$data
series1b= datagen.outlier(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                          n0 = 70, old.version = TRUE,
                          SAOplus=c(50),SAOplusmag=20)$data

#plot series
par(mfrow = c(1,1),mai = c(1, 0.42, 0.1,0.1))
plot(as.ts(series1a,frequency=7),main=''
     ,ylab='',cex.axis = 0.5,bty="n",frame.plot=FALSE)
plot(as.ts(series1b,frequency=7),main=''
     ,ylab='', cex.axis = 0.5,bty="n",frame.plot=FALSE)

####Setting as in Marczak & Proietti ####################################################################
##############################################################v###################
##****in progresss doesnt work now

#simulate your LLM
T = 100
pars <- c(var1 = 15099, var2 = 1469.1,a01 = 1120)
m <- stsm.model(model = "local-level", y = ts(seq(100)), 
                pars = pars, nopars = NULL)
ss <- char2numeric(m)

iter = 10 #set the number of simulations
y = matrix(0,T,iter)
mu = matrix(0,T,iter)

res=rep(0,length(iter))
out = list() 

set.seed(123)
for (i in seq(iter)){
  y[,i] <- datagen.stsm(n = T, model = list(Z = ss$Z, T = ss$T, H = ss$H, Q = ss$Q, a0 = ss$a0), 
                        n0 = 0, old.version = TRUE)$data
  
  #simulation of two AO and one LS
  y_m=as.ts(y[2:100,i])
  y_m[15] = 9*var(y_m)
  y_m[45] = -8*var(y_m)
  #output list
  out[[i]]=isat(y_m)
  plot(out[[i]])
}
###############################################################
# old code - to be deleted
###############################################################
# basic local level model with normally distributed disturbances
###############################################################
# structural time series 
model_structural=SSModel(Nile~ -1 +
                           SSMtrend(degree =1, Q = list(matrix(15099)))
                         , H = matrix(1469.1))
fit_structural=fitSSM(model_structural, inits = c(0, 0),
                      method = "BFGS",hessian=T)
fit_structural$model["Q"] # display Q matrices for disturbances for observation eq and trend eq 
#(we set the variance of these disturbances to zero)
fit_structural$model["H"] 

#state smoothing
out=KFS(model_structural,smoothing=c("state","mean","disturbance"),simplify = F)
#plot only the level of the state
ts.plot(Nile, out$a[-1,1], out$att[,1], out$alpha[,1], col = 1:3)

#disturbance smoothing
#out_dist=KFS(fit_structural$model,smoothing="disturbance",simplify = F)

###############################################
###            outlier detection
##############################################

# APPROACH from the DK book
par(mfrow=c(1,2))
ts.plot(log(kpi.zoo), out$a[-1,1], out$att[,1], out$alpha[,1], col = 1:3)
plot(out_dist$etahat[,1],type="b") #structural break (state; equation)
ts.plot(log(kpi.zoo), out$a[-1,1], out$att[,1], out$alpha[,1], col = 1:3)
plot(out_dist$epshat[,1],type="b") #additive outlier (observation eq)
plot(out$e,type="b")

plot(rstandard(out,type="recursive"),type='h')
plot(rstandard(out,type="recursive",standardization_type ="cholesky"),type='h')
plot(rstandard(out,type="pearson",standardization_type ="cholesky"),type='h')

plot(rstandard(out,type="state"),type='h')
plot(rstandard(out,type="state",standardization_type ="cholesky"),type='h') #same output as with marginal standardisation

#score marginal log-likelihood
out$logLik

#calculation of ut from disturbance smoothing 
u=out$v/t(out$F)-(out$r[-1])/(out$K[1,1,])
D=1/t(out$F)+(out$K[1,1,])*out$N[1,1,-1]*(out$K[1,1,])

#################################################
### estimating d(t,j) given by (15) from the paper
##################################################
dj1=matrix(NA,length(Nile),2)
dj2=matrix(NA,length(Nile),2)
theta_tilde1=matrix(NA,length(Nile),2)
theta_tilde2=matrix(NA,length(Nile),2)
for (t in 1:length(Nile)){
  shock=rep(0,length(Nile))
  shock[t]=1
  
  ### dummy in observation equation
  model_structural=SSModel(Nile ~-1 + shock+
                             SSMtrend(degree =1, Q = list(matrix(NA)))
                           , H = matrix(NA))

  fit_structural=fitSSM(model_structural, inits = c(0, 0),
                      method = "BFGS",hessian=T)
  #state smoothing
  out_j1=KFS(fit_structural$model,smoothing=c("state","mean","disturbance"),simplify = F)
  
  dj1[t,1]=(fit_structural$model["Q"][1,1,]-out$model["Q"][1,1,])/(fit_structural$optim.out$hessian[1,1])^(-0.5) #is this formula ok????
  dj1[t,2]=(fit_structural$model["H"][1,1,]-out$model["H"][1,1,])/(fit_structural$optim.out$hessian[2,2])^(-0.5)
  theta_tilde1[t,1]=fit_structural$model["Q"][1,1,]+(out_j1$logLik-out$logLik)/(fit_structural$optim.out$hessian[1,1])^(-0.5)
  theta_tilde1[t,2]=fit_structural$model["H"][1,1,]+(out_j1$logLik-out$logLik)/(fit_structural$optim.out$hessian[2,2])^(-0.5)
  
  ### dummy in state equation
  Zt <- matrix(c(1, shock),1,length(Nile)+1)
  Ht <- matrix(NA)
  Tt <- matrix(c(1, rep(0,100), rep(0,100), diag(length(Nile))), length(Nile)+1,length(Nile)+1)       
  Rt <- matrix(c(1, rep(0,100)), length(Nile)+1, 1)
  Qt <- matrix(NA)
  a1 <- matrix(c(1, rep(0,100)), length(Nile)+1, 1)
  P1 <- matrix(0, length(Nile)+1, length(Nile)+1)
  P1inf <- diag(length(Nile)+1) #for regression coef, initial value can be taken both as diffuse or fixed. 
                   #in this case diffuse
  model_structural=SSModel(Nile  ~ -1 +
                             SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1,
                              P1inf = P1inf),
                           H = Ht)
  
  fit_structural=fitSSM(model_structural, inits = c(0, 0),
                        method = "BFGS",hessian=T)
  #state smoothing
  out_j2=KFS(fit_structural$model)
  
  dj2[t,1]=(fit_structural$model["Q"][1,1,]-out$model["Q"][1,1,])/(fit_structural$optim.out$hessian[1,1])^(-0.5)
  dj2[t,2]=(fit_structural$model["H"][1,1,]-out$model["H"][1,1,])/(fit_structural$optim.out$hessian[2,2])^(-0.5)
  theta_tilde2[t,1]=fit_structural$model["Q"][1,1,]+(out_j2$logLik-out$logLik)/(fit_structural$optim.out$hessian[1,1])^(-0.5)
  theta_tilde2[t,2]=fit_structural$model["H"][1,1,]+(out_j2$logLik-out$logLik)/(fit_structural$optim.out$hessian[2,2])^(-0.5)
  
}

#plot the results
par(mfrow=c(3,2))
plot(1:length(Nile),dj1[,1],type='h')
plot(1:length(Nile),dj1[,2],type='h')
plot(1:length(Nile),dj1[,3],type='h')

plot(1:length(Nile),dj2[,1],type='h')
plot(1:length(Nile),dj2[,2],type='h')
plot(1:length(Nile),dj2[,3],type='h')

plot(residuals(out,type="recursive"))
plot(residuals(out,type="response")) #same as pearson standartised
plot(residuals(out,type="state"))
plot(out$etahat) #nonstandardized etahats