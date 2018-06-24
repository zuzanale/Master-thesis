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
out = list() 

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
  
  #output list
  out[[i]]=isat_my(y_m, t.pval=0.001, model.ssmod = "local-level",
                   vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                   optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
  end_time <- Sys.time()
  print(end_time - start_time)
  plot(out[[i]])
}

isat_my(y_m, t.pval=0.001, model.ssmod = "local-level",
                 vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                 optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
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