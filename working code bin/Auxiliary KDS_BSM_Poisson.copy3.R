kafs=function(y, Xt) 
{
  require(dlm)
  require(stsm)
  y=as.ts(y)
  Xt = as.matrix(Xt)
  initpars <- c(var1 =0, var2 =0,var3 = 0,var4 = 0,a01=0,a02=0,
                a03=0,a04=0)
  
  #seas is taken from the seasonality setted in a particular time series
  mairp <- stsm.model(model ="BSM",y = y, 
                      pars = initpars, nopars = NULL, transPars = NULL)
  #res.seas <- maxlik.td.optim(mairp, KF.version = "KFKSDS",
 #                         KF.args = list(P0cov = FALSE), method = "L-BFGS-B")
  #mairp <- set.pars(mairp, res.seas$par)
  
  #maybe estimation of variance for seasonal component should be removed, cuz it takes a lot od time
  buildFun <- function(z) {
    #dlmModPoly(1, dV = exp(x[1])^2, dW = exp(x[2])^2) +dlmModSeas(frequency(y), dV = 1,dW = c(exp(x[2])^2, rep(0, frequency(y) - 2))) # dlmModSeas(frequency(y))
  m=dlmModPoly(2)+dlmModSeas(frequency(y))
  diag(W(m))[1:3]= exp(z[2:4])
  V(m)[1]= exp(z[1])
return(m)
  }
  

  #model with poisson distribution using KFAS
 
model <- SSModel(y ~ Xt +SSMtrend(degree = 2, Q = list(NA,NA))
                            +SSMseasonal(period=7,Q= matrix(1),
                                       sea.type = "dummy")
                            ,distribution = "negative binomial")#u is an alternative for H-you have to figure out

update <- function(pars, model) {
  #diag(model["Q"][,,1])[1] <- exp(pars[1])
  model["Q", etas = "level"] <- exp(pars[1])
  model["Q", etas = "slope"] <- exp(pars[2])
  #diag(model["Q"][,,1])[2] <- exp(pars[2])
  model["Q", etas = "seasonal"] <- exp(pars[3])
  #diag(model["Q"][,,1])[3] <- exp(pars[3])
#  model["u"][1] <- exp(pars[4])
  #model["P1", states = "custom"] <- exp(pars[2])
 model
  }
fit <- fitSSM(model, c(2.692846,2.692846,2.692846),update, method = "BFGS",hessian=TRUE)
fit$model["Q", etas = "level"]
fit$model["Q", etas = "slope"]
fit$model["Q"]
  #res2 <- dlmMLE(y, parm = c(0,0,0,0), build = buildFun, hessian=T, method ="L-BFGS-B")#computes numerical hessian, 
                              #maybe i can change this to compute a numerical one, if i will have enough time
  #res2$vcov
  #dlmModel<- buildFun(res2$par)
 # param=c(V(dlmModel),diag(W(dlmModel))[1:3])
  
  #save information matrix of second derivates of likelihood
  infomat =(fit$optim.out$hessian)^(-1)
 # avarLog =solve(res2$hessian)
  #infomat= diag(exp(res2$par))%*% avarLog %*% diag(exp(res2$par))
  infomat.std.err=sqrt(diag(infomat))
  
  result<- KFAS::KFS(fit$model,convtol= 1e-06)
  out$vcov = result$V[1:NCOL(x),1:NCOL(x),out$n]
  out$fit = result$att
  out$coefficients= result$alphahat[out$n,1:NCOL(x)]
  if (out$k==1){
    out$t.stats =out$coefficients/sqrt(result$V[1:NCOL(x),1:NCOL(x),out$n])
    out$std.error = sqrt(result$V[1:NCOL(x),1:NCOL(x),out$n])
  }
  else{
    out$t.stats =  out$coefficients/sqrt(diag(result$V[1:NCOL(x),1:NCOL(x),out$n]))
    out$std.error = sqrt(diag(result$V[1:NCOL(x),1:NCOL(x),out$n]))
  }
  
 #preparation of the output list
 list( a.upd =a.upd, Q = Q, H = H, f = f, v = v, epshat = epshat, vareps = vareps, etahat = etahat, vareta = vareta, 
      r = r, N = N,V_x=V_x, ahat = ahat, varahat = varahat, Et=Et, ut = ut, ut_star = ut_star,
      N_star = N_star, Rt = Rt, rt_star = rt_star, Dt_star = Dt_star, Dt = Dt,
      St = Sn, bt = -bn, dist = dist, score_changed = score_changed, logl = res2$loglik, pars.change = pars.change,
      old.param = param ,infomat= infomat,infomat.std.err=infomat.std.err
    )
 #score_changed is wrong for some reason but doesn not really matter that much now :(((())))
 }