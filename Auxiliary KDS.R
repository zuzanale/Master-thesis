kafs=function (y, Xt) 
{
  require(dlm)
  require(stsm)
  y=as.ts(y)
  initpars <- c(var1 =0, var2 =0)
  
  mairp <- stsm.model(model = "local-level",y = y, 
                      pars = initpars,lower = c(var1 =0, var2 =0), nopars = NULL, transPars = NULL)
  
 # mairp.change <- stsm.model(model = "local-level",y =y, 
 #                    pars = initpars,lower = c(var1 =0, var2 =0),xreg=Xt, nopars = NULL, transPars = NULL)
  
  #using fisher scoring algorithm for optimization of log-likelihood
  #res3 <- maxlik.td.scoring(mairp,
   #                         KF.args = list(P0cov = FALSE), check.KF.args = TRUE,
   #                         ls = list(type = "optimize", tol = .Machine$double.eps^0.25, cap = 1),
   #                         control = list(maxit = 10000, tol = 0.001, trace = FALSE),debug=TRUE)
 # res2 <-maxlik.td.optim(mairp,
  #                       KF.args = list(P0cov = FALSE), check.KF.args = TRUE,
  #                       method = "L-BFGS-B",
  #                       gr ="analytical",optim.control=list(maxit = 100000000))
  
  buildFun <- function(x) {
    dlmModPoly(1, dV = exp(x[1])^2, dW = exp(x[2])^2)
  }
  
  res2 <- dlmMLE(y, parm = c(0,0), build = buildFun, hessian=T, method ="L-BFGS-B")#computes numerical hessian, 
                              #maybe i can change this to compute a numerical one, if i will have enough time
  #res2$vcov
  dlmModel<- buildFun(res2$par)
  param=c(V(dlmModel),W(dlmModel))
  
 #res.change <-maxlik.td.optim(mairp.change,
  #                         KF.args = list(P0cov = FALSE), check.KF.args = TRUE,
  #                          method = "L-BFGS-B",
  #                          gr ="analytical",optim.control=list(maxit = 100000000))
 #lik.change=res.change$loglik-res2$loglik
 
 # res2 <-maxlik.td.optim(mairp,``
  #                       optim.control=list(maxit = 100000000))
  
  #save information matrix of second derivates of likelihood
  infomat=solve(res2$hessian)
  res2$std.errors<-sqrt(diag(solve(res2$hessian)))
  #infomat = matrix(c(0.0101,-0.0263,-0.0263,0.188),2,2)
  mairp<- set.pars(mairp, pmax(param, .Machine$double.eps))
  ss<- char2numeric(mairp) #produces parameters and properties of our ss model
  #Xt = rep(0,length(y));
  #Xt[10] = 10; # put the shock on the first time index t=1
  
  convergence = c(0.001, length(y))
  tol <- convergence[1]
  maxiter <- convergence[2]
  checkconv <- maxiter < length(y)
  #convit <- ifelse(checkconv,0,NULL)
 
  ss <- char2numeric(mairp)
  c = 0 #mean adjustments
  d = 0 #mean adjustments
  n <- length(y)
  p <- 1
  #ss$a0=0
  #ss$P0=10000000
  a0 <- ss$a0
  P0 <- ss$P0
  Z <- ss$Z
  mT <- ss$T
  H <- ss$H
  Q <- ss$Q
 # H = 4.81
 # Q = 3.64
  tZ <- t(Z)
  tmT <- t(mT)
  notconv <- TRUE
  counter <- 0
  a.pred <- matrix(nrow = n, ncol = length(a0))
  P.pred <- array(NA, dim = c(dim(P0), n))
  a.upd <- matrix(nrow = n + 1, ncol = length(a0))
  P.upd <- array(NA, dim = c(dim(P0), n + 1))
  
  K <- matrix(nrow = n, ncol = length(a0))
  L <- array(NA, dim = c(length(a0), length(a0), n))
  v <- rep(NA, n)
  f <- rep(NA, n)
  lik.contrib <- rep(0, n)
  
  #for auxialiry part
  bt = matrix(nrow = n + 1, ncol = length(a0))
  At <- matrix(nrow = n + 1, ncol = length(a0))
  st <- matrix(nrow = n + 1, ncol = length(a0))
  St <- array(NA, dim = c(dim(P0), n + 1))
  V_x <- rep(NA, n)
  V_aux = matrix(length(V_x),n)
  
  a.upd[1, ] <- a0
  P.upd[, , 1] <- P0
  
  At[1, ] = mT %*% a0 #At je mozno spatne definovane
  st[1, ] <- 0
  St[, , 1] <- 0
  
  t0 = 1
  ###################
  #filtering
  #######################
  for (i in seq_len(n)) {
    a.pred[i, ] <- mT %*% a.upd[i, ]
    if (notconv) {
      P.pred[, , i] <- mT %*% P.upd[, , i] %*% tmT + Q
    }
    else P.pred[, , i] <- P.pred[, , i - 1]
    if (!is.na(y[i])) {
      v[i] <- y[i] - Z %*% a.pred[i, ] - c 
      V_x[i] = -Xt[i] - c -  Z %*% At[i, ]
      #V_aux[i,i] = V_x[i] #really not sure what supposted to be here
      if (notconv) {
        Mt <- tcrossprod(P.pred[, , i], Z)
        f[i] <- Z %*% Mt + H
      }
      else f[i] <- f[i - 1]
      lik.contrib[i] <- log(f[i]) + v[i]^2/f[i]
      K[i, ] <- Mt/f[i]
      a.upd[i + 1, ] <- a.pred[i, ] + K[i, ] * v[i] + d
      At[i + 1, ] =  mT %*% At[i, ] +  K[i, ] *V_x[i] + d;
      if (notconv) {
        P.upd[, , i + 1] <- P.pred[, , i] - tcrossprod(Mt)/f[i]
        K[i, ] <- mT %*% K[i, ]
        L[, , i] <- mT - K[i, ] %*% Z
      }
      else {
        P.upd[, , i + 1] <- P.upd[, , i]
        K[i, ] <- K[i - 1, ]
        L[, , i] <- L[, , i - 1]
      }
      if (checkconv && notconv) {
        if (i == 1) {
          fprev <- f[i] + tol + 1
        }
        if (abs(f[i] - fprev) < tol) {
          if (convit == i - 1) {
            counter <- counter + 1
          }
          else counter <- 1
          convit <- i
        }
        fprev <- f[i]
        if (counter == maxiter) {
          notconv <- FALSE
          convit <- i
        }
      }
      At[i + 1, ] =  mT %*% At[i, ] +  K[i, ] *V_x[i] + d;
      if (i!=1){ 
        st[i, ] = st[i-1, ] +  (t(V_x[i])/f[i])%*%v[i]
        St[, , i]  = St[, , i-1]  + crossprod(V_x[i]/f[i], V_x[i])
        if (St[, , i-1]==0) {
         bt[i, ] = 0
        }else {bt[i, ] = st[i, ]/St[, , i]} 
      }
    } else {
      a.upd[i + 1, ] <- a.pred[i, ]
      P.upd[, , i + 1] <- P.pred[, , i]
      if (!notconv) {
        notconv <- TRUE
        counter <- 1
      }
    }
  }
  convit <- NULL
  v <- ts(v)
  f <- ts(f)
  a.upd <- ts(a.upd[-1, ])
  K <- ts(K)
  tsp(v) <- tsp(f) <- tsp(a.upd) <- tsp(K) <- tsp(y)
  mll <- 0.5 * (n - t0 + 1) * log(2 * pi) + 0.5 * sum(lik.contrib[seq.int(t0, 
                                                                          n)], na.rm = TRUE)
  kf=list(v = v, f = f, K = K, L = L, a.upd = a.upd, P.upd = P.upd[, 
                                                                , -1], a.pred = a.pred, P.pred = P.pred, mll = mll, convit = convit)
  
  
  R <- ss$R
  Q <- ss$Q
  #Q=3.64
  v <- kf$v
  a.pred <- kf$a.pred
  P.pred <- kf$P.pred
  f <- kf$f
  K <- kf$K
  L <- kf$L
  r <- matrix(nrow = n + 1, ncol = ncol(a.pred))
  rt_star <- matrix(nrow = n + 1, ncol = ncol(a.pred))
  N <- array(dim = c(ncol(a.pred), ncol(a.pred), n + 1))
  N_star <- array(dim = c(ncol(a.pred), ncol(a.pred), n + 1))
  
  ahat <- matrix(nrow = n, ncol = ncol(a.pred))
  varahat <- array(NA, dim = c(rep(ncol(a.pred), 2), n))
  r[n + 1, ] <- 0
  rt_star[n + 1, ] <- 0
  N[, , n + 1] <- 0
  N_star[, , n + 1] <- 0
  epshat <- rep(NA, n)
  vareps <- matrix(nrow = n, ncol = 1)
  ut <- rep(NA, n)
  Dt <- matrix(nrow = n, ncol = 1)
  ut_star <- rep(NA, n)
  Dt_star <- matrix(nrow = n, ncol = 1)
  Et <- rep(NA, n)
  Rt <- matrix(nrow = n + 1, ncol = ncol(a.pred))
  Rt[n + 1, ] <- 0
  
  
  etahat <- matrix(nrow = n, ncol = ncol(ss$Q))
  vareta <- array(dim = c(n, ncol(ss$Q), ncol(ss$Q)))
  for (i in seq(n, 1)) {
    ip1 <- i + 1
    r[i, ] <- tZ * v[i]/f[i] + t(L[, , i]) %*% r[ip1, ]
    
    N[, , i] <- crossprod(Z)/f[i] + crossprod(L[, , i], N[, 
                                                          , ip1]) %*% L[, , i]
    
    ahat[i, ] <- a.pred[i, ] + P.pred[, , i] %*% r[i, ]
    varahat[, , i] <- P.pred[, , i] - P.pred[, , i] %*% N[, 
                                                          , i] %*% P.pred[, , i]
    epshat[i] <- H * (v[i]/f[i] - crossprod(cbind(K[i, ]), 
                                            r[ip1, ]))
    vareps[i, ] <- H - H^2 * (1/f[i] + crossprod(cbind(K[i, 
                                                         ]), N[, , ip1]) %*% cbind(K[i, ]))
    tmp <- tcrossprod(Q, R)
    etahat[i, ] <- tmp %*% cbind(r[ip1, ])
    vareta[i, , ] <- Q - tmp %*% N[, , ip1] %*% R %*% Q
    
    ut[i] = v[i]/f[i] - crossprod(cbind(K[i, ]), r[ip1, ])
    Dt[i, ] = 1/f[i] + crossprod(K[i, ], N[, , i]) %*% K[i, ]
    
    ##auxiliary filter part
    #first one should be correct
    #Et[i] = V_aux[,i]/f[i] - crossprod(cbind(K[i, ]), Rt[ip1, ]) #might have to be rewritten to V[i,i]
    Et[i] = V_x[i]/f[i] - crossprod(cbind(K[i, ]), Rt[ip1, ]) #might have to be rewritten to V[i,i]
   
     #auxiliary corrections
    Rt[i, ] = tZ * Et[i] + tmT %*% Rt[ip1, ]
    rt_star[i, ] = r[i, ] + Rt[i, ] %*% bt[n,];
    
   # r[i, ] <- tZ * v[i]/f[i] + t(L[, , i]) %*% r[ip1, ]
    
   # N[, , i] <- crossprod(Z)/f[i] + crossprod(L[, , i], N[,, ip1]) %*% L[, , i]
    
    
    N_star[, , i] <- N[, , i] - tcrossprod(Rt[i, ]/St[, , n], Rt[i, ])
    
   # ut[i] = v[i]/f[i] - crossprod(cbind(K[i, ]), r[ip1, ])
   # Dt[i, ] = 1/f[i] + crossprod(K[i, ], N[, , i]) %*% K[i, ]
    
    ut_star[i] = ut[i] + Et[i] %*% bt[n,]; # bt(T,1) in the paper????
    Dt_star[i, ] =  Dt[i, ] - tcrossprod(Et[i]/St[, , n],Et[i])
  }
  
  r <- ts(matrix(r[-1, ], nrow = nrow(r) - 1, ncol = ncol(r)), 
          frequency = frequency(y), end = end(y))
  rt_star <- ts(matrix(rt_star[-1, ], nrow = nrow(rt_star) - 1, ncol = ncol(rt_star)), 
          frequency = frequency(y), end = end(y))
  
  Rt <- ts(matrix(Rt[-1, ], nrow = nrow(Rt) - 1, ncol = ncol(Rt)), 
          frequency = frequency(y), end = end(y))
  
  N <- array(N[, , -1], dim = dim(N) - c(0, 0, 1))
  N_star <- array(N_star[, , -1], dim = dim(N_star) - c(0, 0, 1))
  ahat <- ts(ahat)
  tsp(ahat) <- tsp(y)
  epshat <- ts(epshat)
  ut =  ts(ut)
  Et =  ts(Et)
  ut_star = ts(ut_star)
  ut_star[is.na(ut_star)] <- 0
  Dt_star[1] <- 0
  etahat <- ts(etahat)
  tsp(epshat) <- tsp(etahat) <- tsp(y)
  
  
  ###################################
  # calculate d #####################
  ################################
  pars.change =rep(0,2)
  old_score = rep(0,2)
  new_score = rep(0,2)
  dist = rep(0,2)
  dist.v2 = rep(0,2)
  #score_changed = rep(0,2)
  
  #old_score[1]=1/2*sum(tcrossprod(ut^2 - Dt,2%*%H))
  #old_score[2]=1/2*sum(tcrossprod(r^2 - N[1,1,],2%*% Q)) #this works just for a model with one state equation
  old_score[1]=1/2*sum(tcrossprod(ut^2 - Dt,1))
  old_score[2]=1/2*sum(tcrossprod(r^2 - N[1,1,],1)) #this works just for a model with one state equation
  
  #new_score[1]=1/2*sum(tcrossprod(ut_star^2 - Dt_star,2%*%H))
  new_score[1]=1/2*sum(tcrossprod(ut_star^2 - Dt_star,1))
  #score_changed[2]=1/2*sum((tcrossprod(results$rt_star,results$rt_star) - results$Nt_star));
  #new_score[2]=1/2*sum(tcrossprod(rt_star^2 - N_star[1,1,], 2%*%Q)) #this works just for a model with one state equation
  new_score[2]=1/2*sum(tcrossprod(rt_star^2 - N_star[1,1,], 1)) #this works just for a model with one state equation
  
  #score_changed[1] = (bt(T-1)-St(T-1)^(-1))/2 * sum(Et.^2 * 1) + bt(T-1) * sum(ut.* Et)
  score_changed = new_score-old_score
  #std.errors <- sqrt(diag(solve(res$hessian)))
  #change in var1(H)
  dist[1] = (score_changed[1]*((infomat[1,1])^(-1)))/((diag(res2$std.errors)[1,1])^(-1))
  #d[1] = (new_score[1]/infomat[1,1])/(infomat[1,1]^(-0.5))
  #change in var2(Q)
  dist[2] =  (score_changed[2]*((infomat[2,2])^(-1)))/((diag(res2$std.errors)[2,2])^(-1))
  #d[2] = (new_score[2]/infomat[2,2])/(infomat[2,2]^(-0.5))
  pars.change[1] = ss$H + (infomat[1,1])^(-1) * score_changed[1]
  pars.change[2] = ss$Q + (infomat[2,2])^(-1) * score_changed[2]
 # dist.v2[1] = pars.change[1]/(diag(res2$std.errors)[1,1])^(1)
  #dist.v2[2] = pars.change[2]/(diag(res2$std.errors)[2,2])^(1)
  
  #compare
  #dist==dist.v2
  #lik.change==score_changed
  
 #preparation of the output list
 list( a.upd =a.upd, Q = Q, H = H, f = f, v = v, epshat = epshat, vareps = vareps, etahat = etahat, vareta = vareta, 
      r = r, N = N, ahat = ahat, varahat = varahat, Et=Et, ut = ut, ut_star = ut_star,
      N_star = N_star, Rt = Rt, rt_star = rt_star, Dt_star = Dt_star, Dt = Dt,
      St = St[,,n], bt = bt[n], dist = dist, score_changed = score_changed, logl = res2$loglik, pars.change = pars.change
      )
}