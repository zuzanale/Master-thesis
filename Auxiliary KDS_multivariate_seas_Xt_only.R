kafs.Xt=function(y, Xt,old.param,infomat, infomat.std.err) 
{
  require(dlm)
  require(stsm)
  y=as.ts(y)
  Xt = as.matrix(Xt)
  #seas is taken from the seasonality setted in a particular time series
  mairp <- stsm.model(model ="BSM",y = y, 
                      pars = initpars, nopars = NULL, transPars = NULL)
  #res.seas <- maxlik.td.optim(mairp, KF.version = "KFKSDS",
 #                         KF.args = list(P0cov = FALSE), method = "L-BFGS-B")
  #mairp <- set.pars(mairp, res.seas$par)
  param=old.param
  mairp<- set.pars(mairp, pmax(param, .Machine$double.eps))
  #mairp@ss$a0=rep(0,frequency(y)-1+2)
  ss<- char2numeric(mairp) #produces parameters and properties of our ss model
  ss$a0 = c(0,rep(0,frequency(y)-1+1))
  #Xt = rep(0,length(y));
  #Xt[10] = 10; # put the shock on the first time index t=1
  
  convergence = c(0.001, length(y))
  tol <- convergence[1]
  maxiter <- convergence[2]
  checkconv <- maxiter < length(y)
  if(checkconv){
    convit=0
  }else{convit=NULL}
  #convit <- ifelse(checkconv,0,NULL)
 

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
  Q <- ss$Q #this one was changed according to pdf
  R <- ss$R
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
  a.upd[1, ] <- a0
  P.upd[, , 1] <- P0
  
  #for auxialiry part
  bt = matrix(nrow = n + 1, ncol = NCOL(Xt))
  At = array(0, dim = c(n + 1, length(a0), NCOL(Xt)))
  #At = matrix(nrow = n + 1, ncol = length(a0)+NCOL(Xt))
  #At <- matrix(nrow = n + 1, ncol = length(a0))
  st <- matrix(nrow = n + 1, ncol = NCOL(Xt))
  St <- array(0, dim = c(NCOL(Xt), NCOL(Xt), n + 1))
  V_x <- matrix(nrow = NROW(Xt), ncol = NCOL(Xt))
  #V_aux = matrix(length(V_x),n)
  
  At[1 ,,] = mT %*% a0 #At je mozno spatne definovane
  st[1, ] <- 0
  St[, , 1] <-  rep(0,NCOL(Xt))
  
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
      v[i] <- y[i] - Z %*% a.pred[i, ]
        Mt <- tcrossprod(P.pred[, , i], Z)
        f[i] <- Z %*% Mt + H
      lik.contrib[i] <- log(f[i]) + v[i]^2/f[i]
      K[i, ] <- Mt/f[i]
      a.upd[i + 1, ] <- a.pred[i, ] + K[i, ] * v[i]
        P.upd[, , i + 1] <- P.pred[, , i] - tcrossprod(Mt)/f[i]
        K[i, ] <- mT %*% K[i, ]
        L[, , i] <- mT - K[i, ] %*% Z
        
        #auxiliary part 
        V_x[i,] = -Xt[i,] - c -  Z%*% At[i,,]
        At[i+1,,] =  mT %*% At[i,,] + crossprod(t(K[i, ]), (V_x[i,]))+ d
        
        #this is a part for bt estimation
        if (i!=1){ 
          st[i, ] = st[i-1, ] +  t(crossprod(t(V_x[i,]),1/f[i])%*%v[i])
          #st[i, ] = st[i-1, ] +  (t(V_x[i,])/rep(f[i],NCOL(Xt)))%*%rep(v[i],NCOL(Xt))
          St[, , i]  = St[, , i-1]  + t(crossprod(t(V_x[i,]),1/f[i])%*%V_x[i,])
          if (St[1,1 , i]==0) {
            bt[i, ] = 0
          }else {
            #bt[i, ] = (St[, , i]^(-1))%*%st[i, ]
            if (NCOL(Xt)!=1){
              bt[i, ] = st[i,]/diag(St[,,i])
            }else{
              bt[i, ] = st[i, ]/St[, , i]
            }
          }
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
    }
    else {
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
  kf=list(v = v, f = f, K = K, L = L, a.upd = a.upd, P.upd = P.upd[,  , -1],
          bt=bt,st=st,St=St,a.pred = a.pred, P.pred = P.pred, mll = mll, convit = convit)
  
  #smoothing
  if (is.null(kf$convit)) {
    convit <- n + 1
    nmconvit <- n
  }else {
    convit <- kf$convit
    nmconvit <- n - convit
  }
  Z <- ss$Z
  mT <- ss$T
  H <- ss$H
  Q <- ss$Q
  tZ <- t(Z)
  tmT <- t(mT)
  v <- kf$v
  a.pred <- kf$a.pred
  P.pred <- kf$P.pred
  f <- kf$f
  L <- kf$L
  r <- matrix(nrow = n + 1, ncol = NCOL(a.pred))
  N <- array(dim = c(NCOL(a.pred), NCOL(a.pred), n + 1))
  ahat <- matrix(nrow = n, ncol = NCOL(a.pred))
  varahat <- array(NA, dim = c(rep(NCOL(a.pred), 2), n))
  r[n + 1, ] <- 0
  N[, , n + 1] <- 0
  
  for (i in seq.int(n, 1)) {
    ip1 <- i + 1
    if (!is.na(y[i])) {
      r[i, ] <- tZ * v[i]/f[i] + t(L[, , i]) %*% r[ip1, 
                                                   ]
      ahat[i, ] <- a.pred[i, ] + P.pred[, , i] %*% r[i, 
                                                     ]
      if (i < convit || i > nmconvit) {
        N[, , i] <- crossprod(Z)/f[i] + crossprod(L[, 
                                                    , i], N[, , ip1]) %*% L[, , i]
        varahat[, , i] <- P.pred[, , i] - P.pred[, , 
                                                 i] %*% N[, , i] %*% P.pred[, , i]
      }
      else {
        N[, , i] <- N[, , ip1]
        varahat[, , i] <- varahat[, , ip1]
      }
    }
    else {
      r[i, ] <- r[ip1, ]
      ahat[i, ] <- a.pred[i, ] + P.pred[, , i] %*% r[i, 
                                                     ]
      N[, , i] <- N[, , ip1]
      varahat[, , i] <- P.pred[, , i] - P.pred[, , i] %*% 
        N[, , i] %*% P.pred[, , i]
    }
    ip2 <- i
    if (ip2 + 1 <= n && is.na(y[ip2])) {
      ip3 <- ip2 + 1
      r[ip3, ] <- NA
      ahat[ip2, ] <- NA
      N[, , ip3] <- NA
      varahat[, , ip2] <- NA
    }
  }
  r <- matrix(r[-1, ], nrow = nrow(r) - 1, ncol = ncol(r))
  N <- array(N[, , -1], dim = dim(N) - c(0, 0, 1))
  if (is.ts(y)) {
    r <- ts(r)
    ahat <- ts(ahat)
    tsp(r) <- tsp(ahat) <- tsp(y)
  }

  ###disturbance smoothing
  if (is.null(kf$convit)) {
    convit <- n + 1
    nmconvit <- n
  }else {
    convit <- kf$convit
    nmconvit <- n - convit
  }
  Z <- ss$Z
  mT <- ss$T
  H <- ss$H
  Q <- ss$Q
  R <- ss$R
  V <- ss$V
  tZ <- t(Z)
  tmT <- t(mT)
  VtR <- tcrossprod(V, R)
  #v <- kf$v
 # f <- kf$f
 # K <- kf$K
 # r <- ks$r
 # N <- ks$N
  epshat <- rep(NA, n)
  vareps <- rep(NA, n)
  etahat <- matrix(nrow = n, ncol = ncol(ss$V))
  vareta <- array(dim = c(n, ncol(ss$V), ncol(ss$V)))
  for (i in seq(n, 1)) {
    cKi <- cbind(K[i, ])
    if (!is.na(y[i])) {
      epshat[i] <- H * (v[i]/f[i] - crossprod(cKi, r[i, 
                                                     ]))
      etahat[i, ] <- VtR %*% cbind(r[i, ])
      if (i < convit || i > nmconvit) {
        vareps[i] <- H - H^2 * (1/f[i] + crossprod(cKi, 
                                                   N[, , i]) %*% cKi)
        vareta[i, , ] <- V - VtR %*% N[, , i] %*% R %*% 
          V
      }
      else {
        vareps[i] <- vareps[i + 1]
        vareta[i, , ] <- vareta[i + 1, , ]
      }
    }
    else {
    }
  }
  if (is.ts(y)) {
    epshat <- ts(epshat)
    etahat <- ts(etahat)
    tsp(epshat) <- tsp(etahat) <- tsp(y)
  }
  
  ###auxiliary smoother
  N_star <- array(dim = c(ncol(a.pred), ncol(a.pred), n + 1))
  N_star[, , n + 1] <-0
  ut <- rep(NA, n)
  Dt <- matrix(nrow = n, ncol = 1)
  ut_star <- rep(NA, n)
  Dt_star <- array(dim = c(NCOL(Xt), NCOL(Xt), n + 1))
  Et <-matrix(nrow = n , ncol = NCOL(Xt))
  Rt=array(dim = c(n + 1, NCOL(a.pred), NCOL(Xt)))
  #Rt <- matrix(nrow = n + 1, ncol = NCOL(Xt))
  rt_star <- matrix(nrow = n + 1, ncol = NCOL(a.pred))
  rt_star[n + 1, ] <- 0
  Rt[n + 1, ,] <- 0
  bn=as.matrix(bt[n,])
  Sn=as.matrix(St[, , n])
  for (i in seq(n, 1)) {
    ip1 <- i + 1
    cKi <- cbind(K[i, ])
    ut[i] = v[i]/f[i] - crossprod(cKi, r[i, 
                                          ])
    Dt[i, ] = 1/f[i] + crossprod(cKi, 
                                 N[, , i]) %*% cKi
    
    ##auxiliary filter part
    #Et[i] = V_aux[,i]/f[i] - crossprod(cbind(K[i, ]), Rt[ip1, ]) #might have to be rewritten to V[i,i]
    Et[i,] =(f[i]^(-1))%*%V_x[i,] - crossprod((K[i, ]), Rt[ip1, ,]) #have to be rewritten to V[i,i]
    Rt[i,, ] = tZ%*% Et[i,] + tmT %*% Rt[ip1, ,]
    ut_star[i] = ut[i] + Et[i,] %*% bn
    rt_star[i, ] = r[i, ] + Rt[i, ,] %*% bn
    N_star[, , i] <- N[, , i] - (Rt[i, ,]%*%(Sn)^(-1))%*% t(Rt[i, ,])
    Dt_star[,,i] =  Dt[i, ] - crossprod(Et[i,]%*%((Sn)^(-1)),t(Et[i,]))
   # r[i, ] <- tZ * v[i]/f[i] + t(L[, , i]) %*% r[ip1, ]
    
   # N[, , i] <- crossprod(Z)/f[i] + crossprod(L[, , i], N[,, ip1]) %*% L[, , i]
  }
  
  rt_star <- ts(matrix(rt_star[-1, ], nrow = nrow(rt_star) - 1, ncol = ncol(rt_star)), 
          frequency = frequency(y), end = end(y))
  
  Rt <- array(Rt[-1, ,], dim=c(NROW(Rt) - 1, NCOL(Rt),NCOL(Xt)))
  #        frequency = frequency(y), end = end(y))
  
  N_star <- array(N_star[, , -1], dim = dim(N_star) - c(0, 0, 1))
  ut =  ts(ut)
  Et =  ts(Et)
  ut_star = ts(ut_star)
  ut_star[is.na(ut_star)] <- 0
  Dt_star[,,101] <- 0
  etahat <- ts(etahat)
  tsp(epshat) <- tsp(etahat) <- tsp(y)
  

  ###################################
  # calculate d #####################
  ################################
  pars.change = matrix(0, 4,NCOL(Xt))
  old_score = rep(0,(length(a0)+1))
  new_score = matrix(0,(length(a0)+1),NCOL(Xt))
  dist = array(0, dim=c(4, 4,NCOL(Xt)))
  score_changed =matrix(0,(length(a0)+1),NCOL(Xt))
  
  #old_score[1]=1/2*sum(tcrossprod(ut^2 - Dt,2%*%H))
  #old_score[2]=1/2*sum(tcrossprod(r^2 - N[1,1,],2%*% Q)) #this works just for a model with one state equation
  old_score[1]=1/2*sum(tcrossprod(ut^2 - Dt,1))
  
  #prepare calculations to pass to score
  prod=matrix(0,n,(length(a0)))
  prod.star=matrix(0,n,(length(a0)))
  for (j in 1:(length(a0))){
    prod[,j]=(r^2)[,j]-N[j,j,]
    prod.star[,j]=(rt_star^2)[-dim(Dt_star)[3],j]-N_star[j,j,-dim(Dt_star)[3]]
  }
  
  #derivations of Q according to states
  Q.der=array(0,dim=c(length(a0),length(a0),length(a0)))
  for (j in 1:(length(a0))){
    Q.der[j,j,j]=1
    old_score[j]=1/2*sum(diag(prod%*% Q.der[j,,])) # the product is not asquare matrix,
    #probably, should be 100x100, it could be a problem with it????
    for (k in 1:NCOL(Xt)){
      new_score[1,k]=1/2*sum((ut_star^2 - Dt_star[k,k,-dim(Dt_star)[3]]))
      new_score[1+j,k]=1/2*sum(diag(prod.star %*% Q.der[j,,]))
    }
  }
  for (k in 1:NCOL(Xt)){
    score_changed[,k] = new_score[,k]-old_score
    #here you have to check what to do if one of elements of infomat.std.err is NaN
    #thus informat was negative
    dist[,,k] =  (((infomat)^(-1))%*%score_changed[1:4,k])%*%((infomat.std.err^(-1))^(-1))
    
    pars.change[1,k] = ss$H + (infomat[1,1])^(-1) %*% score_changed[1,k]
    #pars.change[1,2,k] = ss$H + (infomat[1,2])^(-1) * score_changed[1,k]
    
    pars.change[(2:4),k] =diag(ss$Q[1:3,1:3]) + ((infomat[2:4,2:4])^(-1))%*%score_changed[2:4,k]
    #pars.change[2,1] = ss$Q[1,1] + (infomat[2,1])^(-1) * score_changed[2]
   # pars.change[2,2] = ss$Q + (infomat[2,2])^(-1) * score_changed[2]
  }
  #bug fix
  if(is.na(bn)){
    bn=0.0000001
  }
 #preparation of the output list
 list( a.upd =a.upd, Q = Q, H = H, f = f, v = v, epshat = epshat, vareps = vareps, etahat = etahat, vareta = vareta, 
      r = r, N = N, ahat = ahat, varahat = varahat, Et=Et, ut = ut, ut_star = ut_star,
      N_star = N_star, V_x=V_x, Rt = Rt, rt_star = rt_star, Dt_star = Dt_star, Dt = Dt,
      St = Sn, bt = -bn, dist = dist, score_changed = score_changed, logl = res2$loglik, pars.change = pars.change
    )
 #score_changed is wrong for some reason but doesn not really matter that much now :(((())))
 }