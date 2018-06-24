kfs=function(y, x, tol = 1e-07, LAPACK = FALSE, method = 1) 
{   out <- list()
    out$n <- length(y)
    if (is.null(x)) {
      out$k <- 0
    }else {
      out$k <- NCOL(x)
    }
    ##############################
    #here my part starts
    ###################################
    out$df <- out$n - out$k
    if (out$k > 0) {
      qx <- qr(x, tol, LAPACK = LAPACK)
      out <- c(out, qx)
      out$xtxinv <- chol2inv(qx$qr, LINPACK = FALSE)
      #this part of if its what to do if we have less regresors then observations (shouldnt be same in our case???=> double check with DK book)
      result = kafs(y, Xt=x)
      out$coefficients =  result$bt
      out$t_stats = result$bt/sqrt(solve(result$St))
      out$eps.std = result$ut/sqrt(result$Dt)
      #out$eps.std.star = result$ut_star/sqrt(result$Dt_star)
      out$eta.std= result$r/sqrt(result$N[1,1,])#this has to be fixed for alpha more than 1 dimensions
      out$t_stats_eps = result$epshat/sqrt(result$Dt) #t-statistics for observations
      #out$t_stats_eps_star =((sum(crossprod(result$ut,result$Et),na.rm = TRUE)/sum(crossprod(result$Et,result$Et),na.rm = TRUE)))/result$St^(-1) #t-statistics for observation 
      #out$t_stats_eta = result$etahat/sqrt(result$N[1,1,]) #t-statistics for state
      #out$t_stats_eta_star =((sum(crossprod(result$r, result$Rt),na.rm = TRUE) /sum(crossprod(result$Rt,result$Rt),na.rm = TRUE)))/result$St^(-1) #t-statistics for state 
      
      #out$xtxinv <- chol2inv(qx$qr, LINPACK = FALSE)
     out$vcov = solve(result$St)#is this cov matrix btw????
      #out$vcov =  result$Dt
      out$fit <- result$a.upd #fitted values of y
      
      #this are just contenporary - erase
      out$Dt= result$Dt
      out$ut= result$ut
      out$r= result$r
      out$N= result$N
    }else{
      out$fit <- rep(0, out$n)
    }
    if (method==1){
      out$residuals <-result$epshat 
    }else{
      out$residuals <-result$etahat 
    }
    ######this part needs to be replaces
   # out$residuals2 <- out$residuals^2
   # out$rss <- sum(out$residuals2)
    out$sigma2 <- solve(result$St) # this is a matrix tho(for each state)
   # if (out$k > 0) {
      #out$vcov <- out$sigma2 * out$xtxinv
    #}
    out$logl <- -out$n * log(2 * pi)/2 - (1/2) * sum(log(result$f) + result$v^2/result$f,na.rm = TRUE)
    return(out)
}

############################################
#auxiliary version to replace arx funnction
############################################
kfs.arx=function(y, mc = FALSE, mxreg = NULL, 
                    asym = NULL,  vxreg = NULL, 
                  zero.adj = 0.1,  vcov.type = c("ordinary", 
                                                               "white", "newey-west"), qstat.options = NULL, user.estimator = NULL, 
                  user.diagnostics = NULL, tol = 1e-07, LAPACK = FALSE, plot = NULL) 
{
  
  #debug part
 # zero.adj = 0.1
# mxreg = mXis
# vcov.type = vcov.type
# qstat.options = qstat.options 
#user.diagnostics = user.diagnostics
#vxreg = NULL
#tol = tol
#LAPACK = LAPACK
#      plot = FALSE
  ###
  
  vcov.type <- match.arg(vcov.type)
  y.name <- deparse(substitute(y))
  if (is.zoo(y)) {
    y <- cbind(y)
  }else {
    y <- as.zoo(cbind(y))
  }
  y <- cbind(y)
  if (NCOL(y) > 1) 
    stop("Dependent variable not 1-dimensional")
  if (is.null(y.name)) {
    y.name <- colnames(y)[1]
  }
  if (y.name[1] == "") {
    y.name <- "y"
  }
  y.n <- NROW(y)
  y.index <- index(y)
  y <- coredata(y)
  mX <- NULL
  mXnames <- NULL
 if (identical(as.numeric(mc), 1)) {
    mX <- cbind(rep(1, y.n))
    mXnames <- "mconst"
  }
  tmp <- zoo(cbind(y, mX), order.by = y.index)
  tmp <- na.trim(tmp, sides = "both", is.na = "any")
  y <- tmp[, 1]
  y.n <- NROW(y)
  y.index <- index(y)
  t1 <- y.index[1]
  t2 <- y.index[y.n]
  y <- coredata(y)
  if (!is.null(mX)) {
    mX <- tmp[, 2:NCOL(tmp)]
    mX <- coredata(mX)
    mX <- cbind(mX)
    colnames(mX) <- NULL
  }
  if (!is.null(mxreg) && !identical(as.numeric(mxreg), 0)) {
    mxreg <- as.zoo(cbind(mxreg))
    mxreg.names <- colnames(mxreg)
    if (is.null(mxreg.names)) {
      mxreg.names <- paste("mxreg", 1:NCOL(mxreg), sep = "")
    }
    if (any(mxreg.names == "")) {
      missing.colnames <- which(mxreg.names == "")
      for (i in 1:length(missing.colnames)) {
        mxreg.names[i] <- paste("mxreg", i, sep = "")
      }
    }
    mXnames <- c(mXnames, mxreg.names)
    mxreg <- window(mxreg, start = t1, end = t2)
    mxreg <- cbind(coredata(mxreg))
    mX <- cbind(mX, mxreg)
    tmp <- zoo(cbind(y, mX), order.by = y.index)
    tmp <- na.trim(tmp, sides = "both", is.na = "any")
    y <- tmp[, 1]
    y.n <- NROW(y)
    y.index <- index(y)
    t1 <- y.index[1]
    t2 <- y.index[y.n]
    y <- coredata(y)
    mX <- tmp[, 2:NCOL(tmp)]
    mX <- coredata(mX)
    mX <- cbind(mX)
    colnames(mX) <- NULL
  }
#  if (is.null(qstat.options)) {
#    if (is.null(ar)) {
#      ar.lag <- 1
#    } else {
#      ar.lag <- max(ar) + 1
#    }
#    if (is.null(arch)) {
#      arch.lag <- 1
#    }else {
#      arch.lag <- max(arch) + 1
#    }
#    qstat.options <- c(ar.lag, arch.lag)
  aux <- list()
  aux$y <- y
  aux$y.index <- y.index
  aux$y.name <- y.name
  aux$y.n <- y.n
  if (!is.null(mX)) {
    colnames(mX) <- NULL
    aux$mX <- mX
    aux$mXnames <- mXnames
    aux$mXncol <- NCOL(mX)
  }
# aux$vc <- vc
 aux$zero.adj <- zero.adj
#  aux$vc.adj <- vc.adj
  aux$vcov.type <- vcov.type
  aux$qstat.options <- qstat.options
  aux$user.estimator <- user.estimator
  aux$user.diagnostics <- user.diagnostics
  aux$tol <- tol
  aux$LAPACK <- LAPACK
  sysCall <- sys.call()
  vcov.var <- NULL
  variance.results <- NULL
  meanSpec <- !is.null(mX)
    estMethod <- which(vcov.type == c("none", "none", "ordinary", 
                                      "white", "newey-west"))
    
    estim=list()
    out=NULL
    #this should be rewriteen into tidyverse world
    out$coefficients=t(rep(0,ncol(mX)))
    out$stats=t(rep(0,ncol(mX)))
    out$vcov=t(rep(0,ncol(mX)))
    out$df=t(rep(0,ncol(mX)))
    out$logl=t(rep(0,ncol(mX)))
    out$k=t(rep(0,ncol(mX)))
    out$ut=matrix(0,y.n,ncol(mX))
    out$Dt=matrix(0,y.n,ncol(mX))
    out$r=matrix(0,y.n,ncol(mX))
    out$N=matrix(0,y.n,ncol(mX))
    out$residuals=list()
    colnames(out$coefficients)=colnames(mX)
    colnames(out$stats)=colnames(mX)
    colnames(out$vcov)=colnames(mX)
    colnames(out$df)=colnames(mX)
    colnames(out$logl)=colnames(mX)
    colnames(out$k)=colnames(mX)
    colnames(out$residuals)=colnames(mX)
    colnames(out$ut)=colnames(mX)
    colnames(out$Dt)=colnames(mX)
    colnames(out$r)=colnames(mX)
    colnames(out$N)=colnames(mX)
  
    for (t in 1:ncol(mX)) {#should probably put here a bug fix (if i it is bigger than number of observations, put an error or sth)
      Xt= mX[,t]
      estim[[t]] <- kfs(y,Xt, tol = tol, LAPACK = LAPACK, method = 1)
      out$coefficients[1,t]= estim[[t]]$coefficients
      out$stats[1,t]= estim[[t]]$t_stats
      out$vcov[1,t]= estim[[t]]$vcov
      out$df[1,t]= estim[[t]]$df
      out$logl[1,t]= estim[[t]]$logl
      out$k[1,t]= estim[[t]]$k
      
      #to be erased
      out$ut[,t]= estim[[t]]$ut
      out$Dt[,t]= estim[[t]]$Dt
      out$r[,t]= estim[[t]]$r
      out$N[,t]= estim[[t]]$N[1,1,]

    }
    out$residuals= estim[[ncol(mX)]]$eps.std
    out$k = unique(t(out$k))
    out$logl = unique(t(out$logl))
    out$qr <- NULL
    out$rank <- NULL
    out$qraux <- NULL
    out$pivot <- NULL
    out$xtxinv <- NULL
    out$residuals2 <- NULL
    out$var.fit <- rep(out$vcov, aux$y.n)
    #out$std.residuals <- out$residuals/sqrt(out$vcov)
    out$std.residuals <-as.ts(out$residuals)
     aux$loge2.n <- aux$y.n
  if (meanSpec) {
    if (length(mXnames)!=1){
      out$vcov=diag(as.vector(t(out$vcov)))
    }
    stderrs <- sqrt(diag(out$vcov)) #vcov is Sn, might need to be changed for Dt or Nt
    
   # print(mXnames)
  #print(str(mXnames))
   # print(out$vcov)
  #print(str(out$vcov))
    
    colnames(out$vcov) <- mXnames
    rownames(out$vcov) <- mXnames
    t.stat=out$stats
   # t.stat <- out$coefficients/stderrs
    p.val <- pt(abs(t.stat), out$df, lower.tail = FALSE) * 
      2
    out$mean.results <- as.data.frame(cbind(t(out$coefficients), 
                                            (stderrs), t(t.stat), t(p.val)))
    colnames(out$mean.results) <- c("coef", "std.error", 
                                    "t-stat", "p-value")
    rownames(out$mean.results) <- mXnames
  }
    outNames <- names(out)
    whereIs <- which(outNames == "vcov")
    if (length(whereIs) > 0) {
      names(out)[whereIs] <- "vcov.mean"
    }
  out$diagnostics <- diagnostics(out, normality.JarqueB = 0, 
                                 user.fun = user.diagnostics, verbose = TRUE)
 
 # out$mean.fit <- zoo(out$mean.fit, order.by = y.index)
  out$residuals <- zoo(out$residuals, order.by = y.index)
  if (!is.null(out$var.fit)) {
    out$var.fit <- zoo(out$var.fit, order.by = y.index)
  }
  if (!is.null(out$ustar.residuals)) {
    out$ustar.residuals <- zoo(out$ustar.residuals, order.by = y.index)
  }
  if (!is.null(out$std.residuals)) {
    out$std.residuals <- zoo(out$std.residuals, order.by = y.index)
  }
  out <- c(list(call = sysCall, date = date(), aux = aux), 
           out)
  class(out) <- "arx"
  return(out)
}
