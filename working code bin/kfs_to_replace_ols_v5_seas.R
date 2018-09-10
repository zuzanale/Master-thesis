kfs=function(y, x, tol = 1e-07, LAPACK = FALSE, method = 1) 
{  
out <- list()
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
  
  #KFAS package with combination of dlm(for square root filtering for stability used)
  #instead of auxialiary filter from atkinson et al.
  require(dlm)
  require(stsm)
  require(KFAS)
  require(gets)
  y=as.ts(y)

  # Filtering and state smoothing
  result<- kafs(y,x)
  out$c= result$bt #a constant to be included in observation level
  out$infomat.std.err=result$infomat.std.err
  out$infomat=result$infomat
  out$old.param=result$old.param
  out$pars.change=result$pars.change
  out$score_changed=result$score_changed
  #out$residuals = result$v/sqrt(result$f)
  
  #this is nonstandartized version
  if (NCOL(result$V_x)>1){
    out$residuals = (result$v+rowSums(result$V_x,na.rm = TRUE))/sqrt(result$f) #plus or minus here?
  }else{
    out$residuals = (result$v+result$V_x)/sqrt(result$f) #plus or minus here?
  }
  #auxiliary residuals
  out$std.eps = result$epshat/sqrt(result$vareps)
  #out$std.eps =  result$Dt^(-1/2)*result$ut #result$ut/sqrt(result$Dt)
  out$std.state=matrix(0,out$n,3)
  out$std.state[,1] =result$etahat[,1]/sqrt(result$vareta[,1,1])
  out$std.state[,2] =result$etahat[,2]/sqrt(result$vareta[,2,2])
  out$std.state[,3] =result$etahat[,3]/sqrt(result$vareta[,3,3])
  #if ((x==rep(1,out$n))&&(NCOL(x)==1)){
    #out$fit = rowSums(result$a.upd[,c(1,3)],na.rm = TRUE)+out$c 
   # out$mean.fit <-rowSums(result$ahat[,c(1,3)],na.rm = TRUE)+out$c
  #}else{
    out$fit = rowSums(result$a.upd[,c(1,3)],na.rm = TRUE)
    out$mean.fit <-rowSums(result$ahat[,c(1,3)],na.rm = TRUE)
  #}
  out$vcov=(result$St)^(-1)
  if (out$k==1){
    out$t.stats =out$c/sqrt((result$St)^(-1))
    out$std.error = sqrt((result$St)^(-1))
  }
  else{
    out$t.stats =out$c/sqrt(diag((result$St)^(-1)))
    out$std.error = sqrt(diag((result$St)^(-1)))
  }
 #if and intercept is included in observation eq
  #rowSums(result$alphahat
  #out$mean.fit <- result$alphahat[,"level"] #fitted v
}else{
  # without regressor
  result<- kafs.without(y)
  out$c= 0
  out$infomat.std.err=result$infomat.std.err
  out$infomat=result$infomat
  out$old.param=result$old.param
  out$pars.change=result$pars.change
  out$score_changed=result$score_changed
  #out$residuals = result$v/sqrt(result$f)
  
  #nonstandartized version
  out$residuals= result$v/sqrt(result$f)
  #out$residuals = result$v
  #auxiliary residuals
  out$std.eps = result$epshat/sqrt(result$vareps)
  #out$std.eps =  result$Dt^(-1/2)*result$ut #result$ut/sqrt(result$Dt)
  out$std.state=matrix(0,out$n,3)
  out$std.state[,1] =result$etahat[,1]/sqrt(result$vareta[,1,1])
  out$std.state[,2] =result$etahat[,2]/sqrt(result$vareta[,2,2])
  out$std.state[,3] =result$etahat[,3]/sqrt(result$vareta[,3,3])
  out$fit = rowSums(result$a.upd[,c(1,3)],na.rm = TRUE) #if and intercept is included in observation eq
  
  out$mean.fit <-rowSums(result$ahat[,c(1,3)],na.rm = TRUE) #if and intercept is included in observation eq
}
out$residuals[is.na(out$residuals)] <-0
  out$std.state[is.na(out$std.state)] <-0
  out$std.eps[is.na(out$std.eps)] <-0
out$logl=result$logl #version from KFAS
#out$logl <- -out$n * log(2 * pi)/2 - (1/2) * sum(log(t(result$F)) + result$v^2/t(result$F),na.rm = TRUE)
return(out)
}

###########################################################################
#################version without likelihood estimation ###################
###########################################################################
kfs.XT=function(y, x, in.old.param, infomat,infomat.std.err,tol = 1e-07, LAPACK = FALSE, method = 1) 
{  
  out <- list()
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
    
    #KFAS package with combination of dlm(for square root filtering for stability used)
    #instead of auxialiary filter from atkinson et al.
    require(dlm)
    require(stsm)
    require(KFAS)
    require(gets)
    y=as.ts(y)
    
    # Filtering and state smoothing
    result<- kafs.Xt(y,x, in.old.param,infomat,infomat.std.err)
    out$c= result$bt #a constant to be included in observation level
    out$coefficients=result$bt
    out$infomat.std.err=result$infomat.std.err
    out$infomat=result$infomat
    out$old.param=result$old.param
    out$pars.change=result$pars.change
    
    #out$resid = result$v/sqrt(result$f)
    
    #nonstandartized version
    if (NCOL(result$V_x)>1){
      out$residuals = (result$v+rowSums(result$V_x,na.rm = TRUE))/sqrt(result$f) #plus or minus here?
    }else{
      out$residuals = (result$v+result$V_x)/sqrt(result$f) #plus or minus here?
    }
    
    out$logl= result$logl
    #auxiliary residuals
    out$std.eps = result$epshat/sqrt(result$vareps)
    #out$std.eps =  result$Dt^(-1/2)*result$ut #result$ut/sqrt(result$Dt)
    out$std.state=matrix(0,out$n,3)
    out$std.state[,1] =result$etahat[,1]/sqrt(result$vareta[,1,1])
    out$std.state[,2] =result$etahat[,2]/sqrt(result$vareta[,2,2])
    out$std.state[,3] =result$etahat[,3]/sqrt(result$vareta[,3,3])
    out$fit = rowSums(result$a.upd[,c(1,3)],na.rm = TRUE) #if and intercept is included in observation eq
    out$vcov=(result$St)^(-1)
    if (out$k==1){
      out$t.stats =out$c/sqrt((result$St)^(-1))
      out$std.error = sqrt((result$St)^(-1))
    }
    else{
      out$t.stats =out$c/sqrt(diag((result$St)^(-1)))
      out$std.error = sqrt(diag((result$St)^(-1)))
    }
    out$mean.fit <-rowSums(result$ahat[,c(1,3)],na.rm = TRUE)#if and intercept is included in observation eq
    #rowSums(result$alphahat
    out$logl <- -out$n * log(2 * pi)/2 - (1/2) * sum(log((result$f)) + result$v^2/(result$f),na.rm = TRUE)
  }else{
    # without regressor
    result<- kafs.without(y)
    out$c= 0
    out$infomat.std.err=result$infomat.std.err
    out$infomat=result$infomat
    out$old.param=result$old.param
    out$pars.change=result$pars.change
    out$residuals = result$v/sqrt(result$f)
    # out$residuals = result$v
    #auxiliary residuals
    out$std.eps = result$epshat/sqrt(result$vareps)
    #out$std.eps =  result$Dt^(-1/2)*result$ut #result$ut/sqrt(result$Dt)
    out$std.state=matrix(0,out$n,3)
    out$std.state[,1] =result$etahat[,1]/sqrt(result$vareta[,1,1])
    out$std.state[,2] =result$etahat[,2]/sqrt(result$vareta[,2,2])
    out$std.state[,3] =result$etahat[,3]/sqrt(result$vareta[,3,3])
    out$fit = rowSums(result$a.upd[,c(1,3)],na.rm = TRUE)#if and intercept is included in observation eq
    
    out$mean.fit <-rowSums(result$ahat[,c(1,3)],na.rm = TRUE) #if and intercept is included in observation eq
  }
  out$residuals[is.na(out$residuals)] <-0
  out$std.state[is.na(out$std.state)] <-0
  out$std.eps[is.na(out$std.eps)] <-0
  out$logl=result$logl #version from KFAS
  out$logl <- -out$n * log(2 * pi)/2 - (1/2) * sum(log((result$f)) + result$v^2/(result$f),na.rm = TRUE)
  return(out)
}

############################################
#auxiliary version to replace arx funnction
############################################
kfs.arx=function(y, mc = FALSE, mxreg = NULL, 
                 asym = NULL,  vxreg = NULL, 
                 zero.adj = 0.1,  vcov.type = c("ordinary", 
                                                "white", "newey-west"), qstat.options = NULL, user.estimator = NULL, 
                 user.diagnostics = NULL, tol = 1e-07, LAPACK = FALSE, plot = NULL,method.kfs=TRUE) 
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
  # mxreg = as.matrix(mxreg)
  vcov.type <- match.arg(vcov.type)
  y.name <- deparse(substitute(y))
  freq=frequency(y)
  if (is.zoo(y)) {
    y <- cbind(y)
  }else {
    y <- zoo(cbind(y),frequency=freq)
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
 # y <- coredata(y)
  mX <- NULL
  mXnames <- NULL
  if (identical(as.numeric(mc), 1)) {
    mX <- cbind(rep(1, y.n))
    mXnames <- "mconst"
  } 
  if (is.null(mX)){
    tmp <- zoo(y, order.by = y.index,frequency = freq)
  }else{
    tmp <- zoo(cbind(y, mX), order.by = y.index,frequency = freq)
  }
  #tmp <- zoo(cbind(y, mX), order.by = y.index,frequency = freq)
  tmp <- na.trim(tmp, sides = "both", is.na = "any")
  y <- tmp[, 1]
  y.n <- NROW(y)
  y.index <- index(y)
  t1 <- y.index[1]
  t2 <- y.index[y.n]
 # y <- coredata(y)
  if (!is.null(mX)) {
    mX <- tmp[, 2:NCOL(tmp)]
    #mX <- coredata(mX)
    mX <- cbind(mX)
    colnames(mX) <- NULL
  }
  if (!is.null(mxreg) && !identical(as.numeric(mxreg), 0)) {
    mxreg <- as.zoo(cbind(mxreg),frequency = freq)
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
    mxreg <- cbind(mxreg)
    if (is.null(mX)){
      mX <- cbind(mxreg)
    }else{
      mX <- cbind(mX, mxreg)
    }
    tmp <- zoo(cbind(y, mX), order.by = y.index,frequency = freq)
    tmp <- na.trim(tmp, sides = "both", is.na = "any")
    y <- tmp[, 1]
    y.n <- NROW(y)
    y.index <- index(y)
    t1 <- y.index[1]
    t2 <- y.index[y.n]
    #y <- coredata(y)
    mX <- tmp[, 2:NCOL(tmp)]
    #mX <- coredata(mX)
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
  
  out <- kfs(y,mX, tol = tol, LAPACK = LAPACK, method = 1)
  out$qr <- NULL
  out$rank <- NULL
  out$qraux <- NULL
  out$pivot <- NULL
  out$xtxinv <- NULL
 # out$residuals2 <- NULL
  out$var.fit <- rep(out$vcov, aux$y.n)
  #out$std.residuals <- out$residuals/sqrt(out$vcov)
    out$std.residuals <-as.ts(out$residuals)
    out$std.residuals.eps <-as.ts(out$std.eps)
    out$std.residuals.state <-as.ts(out$std.state)
    out$coefficients=out$c
    aux$loge2.n <- aux$y.n
  if (meanSpec) { #meanSpec is TRUE if there is some regressor included
    stderrs <-out$std.error  #vcov is Sn, might need to be changed for Dt or Nt
    
    # print(mXnames)
    #print(str(mXnames))
    # print(out$vcov)
    #print(str(out$vcov))
    out$vcov = as.matrix(out$vcov)
    colnames(out$vcov) <- mXnames
    rownames(out$vcov) <- mXnames
    stderrs = out$std.error 
    t.stat <- out$t.stats
    p.val <- pt(abs(t.stat), out$df, lower.tail = FALSE) * 
      2
    
    out$mean.results <- as.data.frame(cbind((out$coefficients), 
                                            (stderrs), (t.stat), (p.val)))
    colnames(out$mean.results) <- c("coef", "std.error", 
                                    "t-stat", "p-value")
    rownames(out$mean.results) <- mXnames
  }
  if (is.null(user.estimator)) {
    outNames <- names(out)
    whereIs <- which(outNames == "vcov")
    if (length(whereIs) > 0) {
      names(out)[whereIs] <- "vcov.mean"
    }
    whereIs <- which(outNames == "fit")
    names(out)[whereIs] <- "mean.fit"
  }
  out$diagnostics <- diagnostics(out, normality.JarqueB = 0, 
                                 user.fun = user.diagnostics, verbose = TRUE)
  
  out$mean.fit <- zoo(out$mean.fit, order.by = y.index)
  out$std.residuals<- zoo(out$std.residuals, order.by = y.index)
  out$std.residuals.eps <-zoo(out$std.residuals.eps, order.by = y.index)
  out$std.residuals.state <-zoo(out$std.residuals.state, order.by = y.index)
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
  if (is.null(plot)) {
    plot <- getOption("plot")
    if (is.null(plot)) {
      plot <- FALSE
    }
  }
  if (plot) {
    plot.arx(out)
  }
  return(out)
}

