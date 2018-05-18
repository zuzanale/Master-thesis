#install.packages("datasets.load")
library(datasets)
View(Nile) #view the NIle data
plot(Nile)

#exports=read.table(file = 'year_origin_destination_sitc_rev2.tsv', sep = '\t', header = TRUE)
#View(exports)

require(KFAS)
###############################################################
# basic local level model with normally distributed disturbances
###############################################################
# structural time series 
model_structural=SSModel(Nile~ -1 +
                           SSMtrend(degree =1, Q = list(matrix(NA)))
                         , H = matrix(NA))
fit_structural=fitSSM(model_structural, inits = c(0, 0),
                      method = "BFGS",hessian=T)
fit_structural$model["Q"] # display Q matrices for disturbances for observation eq and trend eq 
#(we set the variance of these disturbances to zero)
fit_structural$model["H"] 

#state smoothing
out=KFS(fit_structural$model,smoothing=c("state","mean","disturbance"),simplify = F)
#plot only the level of the state
ts.plot(log(kpi.zoo), out$a[-1,1], out$att[,1], out$alpha[,1], col = 1:3)

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