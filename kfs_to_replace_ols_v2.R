kfs=function(y, x, tol = 1e-07, LAPACK = FALSE, method = 1) 
{  require(KFAS)
  require(stsm)
  require(dlm)
  
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
      y=as.ts(y)
        buildFun <- function(z) {
       m= dlmModPoly(1, dV = exp(z[1])^2, dW = exp(z[2])^2) #+ dlmModReg(x)
       return(m)
      }
      
      res2 <- dlmMLE(y, parm = c(0,0), build = buildFun, hessian=T, method ="L-BFGS-B")#computes numerical hessian, 
      #maybe i can change this to compute a numerical one, if i will have enough time
      #res2$vcov
      
      dlmModel<- buildFun(res2$par)
      
      model <- KFAS:::SSModel(y ~ x +
                               SSMtrend(1, Q = list(W(dlmModel)[1,1])), H = matrix(V(dlmModel)))
      result = KFAS:::KFS(model, filtering = "state", smoothing = c("state"),simplify=T)
      
      out$coefficients = result$alphahat[ out$n,1:ncol(x)]
      out$vcov =result$V[ 1:ncol(x), 1:ncol(x),out$n] #vcov matrix
      out$std.error = sqrt(diag(result$V[ 1:ncol(x), 1:ncol(x),out$n]))

      out$fit <- result$alphahat[out$n,(ncol(x)+1):(ncol(result$alphahat))] #smoothed values for unobserved componetns
      
    }else{
      out$fit <- rep(0, out$n)
    }
    if (method==1){
      out$residuals <-result$v
    }else{
      out$residuals <-result$v
      #ut$residuals <-result$epshat #be careful about this, if it is used somewhere else it should be mentioned
    }
    ######this part needs to be replaces
   # out$residuals2 <- out$residuals^2
   # out$rss <- sum(out$residuals2)
   # out$sigma2 <- solve(result$St) # this is a matrix tho(for each state)
   # if (out$k > 0) {
      #out$vcov <- out$sigma2 * out$xtxinv
    #}
    out$logl2 = result$logLik
    out$logl <- -out$n * log(2 * pi)/2 - (1/2) * sum(log(t(result$F)) + result$v^2/t(result$F),na.rm = TRUE)
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
    
    out=NULL
    out <- kfs(y,Xt, tol = tol, LAPACK = LAPACK, method = 1)
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
    stderrs <-out$std.error  #vcov is Sn, might need to be changed for Dt or Nt
    
   # print(mXnames)
  #print(str(mXnames))
   # print(out$vcov)
  #print(str(out$vcov))
    
    colnames(out$vcov) <- mXnames
    rownames(out$vcov) <- mXnames
    stderrs = out$std.error
    t.stat <- out$coefficients/stderrs
    p.val <- pt(abs(t.stat), out$df, lower.tail = FALSE) * 
      2
    out$mean.results <- as.data.frame(cbind((out$coefficients), 
                                            (stderrs), (t.stat), (p.val)))
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
