isat_my=function(y, mc = TRUE, ar = NULL, ewma = NULL, mxreg = NULL, 
               iis = TRUE, sis = TRUE, tis = FALSE, uis = FALSE, blocks = NULL, 
               ratio.threshold = 0.8, max.block.size = 30, t.pval = 0.001, 
               wald.pval = t.pval, vcov.type = c("ordinary", "white", "newey-west"), 
               do.pet = FALSE, ar.LjungB = NULL, arch.LjungB = NULL, normality.JarqueB = NULL, 
               user.diagnostics = NULL, info.method = c("sc", "aic", "hq"), 
               include.gum = NULL, include.1cut = FALSE, include.empty = FALSE, 
               max.paths = NULL, parallel.options = NULL, turbo =FALSE, 
               tol = 1e-07, LAPACK = FALSE, max.regs = NULL, print.searchinfo = TRUE, 
               plot = NULL, alarm = FALSE,
               #for my version of the code (for a local level model)
               model.ssmod = "local-level", pars.ssmod =  c(var1 = log(var(y)), var2 = log(var(y)),
                                                            P0=10000000,a0=y[1]),
               nopars.ssmod = NULL, cpar.ssmod = NULL, 
               xreg.ssmod = NULL, lower.ssmod = NULL, upper.ssmod  = NULL, transPars.ssmod = NULL,
               ssd.ssmod = FALSE, sgfc.ssmod = FALSE,
               kf.args = list(P0cov = FALSE), # if we would set this one to TRUE, we would say that
               # P0 for max likelihood with KF is not diagonal
               method.kfs ="L-BFGS-B",
               gr.kfs = c("numerical", "analytical"), 
               optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE) #this should be later rewritten 
               #into hessian = FALSE, can really slow down the whole process
               ) 
{
  #(y_m, t.pval=0.001, model.ssmod = "local-level",
 #  vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
  # optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
  #################################################
  
 
#  mc = TRUE
#  ar = NULL
#  ewma = NULL
#  mxreg = NULL 
#  iis = TRUE
#  sis = TRUE
#  tis = FALSE
#  uis = FALSE
#  blocks = NULL 
#  ratio.threshold = 0.8 
 # max.block.size = 30
 # t.pval = 0.001
 # wald.pval = t.pval
 # vcov.type = "ordinary"
 # do.pet = FALSE
 # ar.LjungB = NULL
 # arch.LjungB = NULL
 # normality.JarqueB = NULL
 # user.diagnostics = NULL
#  info.method ="aic"
#  include.gum = NULL
 # include.1cut = FALSE
 # include.empty = FALSE
 # max.paths = NULL
 # parallel.options = NULL
#  turbo = FALSE
#  tol = 1e-07
#  LAPACK = FALSE
#  max.regs = NULL
 # print.searchinfo = TRUE
#  plot = NULL
#  alarm = FALSE
 # model.ssmod = "local-level"
  #pars.ssmod <- c(var1 =var(y), var2 =var(y))
 # pars.ssmod =  c("var1" = var(y), "var2" = var(y),"P01"=10000000,"a01"=y[1])
 # nopars.ssmod = NULL
#  cpar.ssmod = NULL
 # xreg.ssmod = NULL
 # lower.ssmod = NULL
 # upper.ssmod  = NULL
#  transPars.ssmod = NULL
#  ssd.ssmod = FALSE
#  sgfc.ssmod = FALSE
#  kf.args = list(P0cov = FALSE) # if we would set this one to TRUE, we would say that
  # P0 for max likelihood with KF is not diagonal
#  method.kfs ="L-BFGS-B"
#  gr.kfs =  "analytical"
#  optim.kfs =list(maxit = 100000000) #this should be later rewritten 
  #into hessian = FALSE, can really slow down the whole process
  
  
  
  ###########################################
  require(zoo)
  require(stsm)
  require(dlm)
  isat.call <- sys.call()
  vcov.type <- match.arg(vcov.type)
  info.method <- match.arg(info.method)
  if (!is.null(include.gum)) {
    warning("The 'include.gum' argument is ignored (temporarily deprecated in isat)")
  }
  include.gum <- TRUE
  olsMethod <- switch(vcov.type, ordinary = 3, white = 4, `newey-west` = 5)
  if (!is.null(max.paths) && max.paths < 1) {
    stop("'max.paths' cannot be smaller than 1")
  }
  if (!is.null(parallel.options)) {
    if (is.numeric(parallel.options)) {
      clusterSpec <- parallel.options
    }
    OScores <- detectCores()
    if (parallel.options > OScores) {
      stop("parallel.options > number of cores/threads")
    }
  }
  #my own version for state space modelling
 # (y, mc = mc, ar = ar, ewma = ewma, mxreg = mxreg, 
  #  vcov.type = vcov.type, qstat.options = NULL, user.diagnostics = user.diagnostics, 
  #  tol = tol, LAPACK = LAPACK, plot = FALSE)
  y=as.ts(y)
  mod <- kfs.arx(y, mc = mc, mxreg = mxreg, 
             vcov.type = vcov.type, qstat.options = NULL, user.diagnostics = user.diagnostics, 
             tol = tol, LAPACK = LAPACK, plot = FALSE)
  #copy these for future analysis
  pure=NULL
  #if constant is not significant, we do not include it
  if (mod$mean.results$`p-value`>0.05)
  {
    mod$c=0
    param= mod$old.param
  }else{param= t(mod$pars.change)}
  pure$std.eps=mod$std.eps
  pure$std.eta=as.ts(mod$std.state[,1]) #for the level
  pure$std.zeta=as.ts(mod$std.state[,2]) #for the slope
  pure$std.omega=as.ts(mod$std.state[,3]) #for seasonality
 # mod <- arx(y, mc = mc, ar = ar, ewma = ewma, mxreg = mxreg, 
 #            vcov.type = vcov.type, qstat.options = NULL, user.diagnostics = user.diagnostics, 
 #            tol = tol, LAPACK = LAPACK, plot = FALSE)
  y <- mod$aux$y
  y.n <- mod$aux$y.n
  y.index <- mod$aux$y.index
  y.index.as.char <- as.character(y.index)
  y.name <- mod$aux$y.name
  mX <- mod$aux$mX
  mXnames <- mod$aux$mXnames
  colnames(mX) <- mXnames
  mXncol <- mod$aux$mXncol
  vcov.type <- mod$aux$vcov.type
  qstat.options <- mod$aux$qstat.options
  if (is.null(mX)){
    mxkeep <- NULL
  }else{
    mxkeep <- 1:mXncol
  }
  arLjungB <- NULL
  if (!is.null(ar.LjungB)) { #two item list with lag and pval for Ljung Box test for serial correlation in the standardised residuals.
    arLjungB <- c(NA, ar.LjungB$pval)
    if (is.null(ar.LjungB$lag)) {
      arLjungB[1] <- qstat.options[1]
    }else {
      arLjungB[1] <- ar.LjungB$lag
    }
  }
  archLjungB <- NULL
  if (!is.null(arch.LjungB)) {
    archLjungB <- c(NA, arch.LjungB$pval)
    if (is.null(arch.LjungB$lag)) {
      archLjungB[1] <- qstat.options[2]
    }else {
      archLjungB[1] <- arch.LjungB$lag
    }
  }
  ISmatrices <- list()
  if (iis) {
    mIIS <- matrix(0, y.n, y.n)
    diag(mIIS) <- 1
    colnames(mIIS) <- paste("iis", y.index.as.char, sep = "")
    ISmatrices <- c(ISmatrices, list(IIS = mIIS))
  }
  if (sis) {
    mSIS <- matrix(0, y.n, y.n)
    loop.indx <- 1:y.n
    tmp <- function(i) {
      mSIS[i, 1:i] <<- 1
    }
    tmp <- sapply(loop.indx, tmp)
    colnames(mSIS) <- paste("sis", y.index.as.char, sep = "")
    mSIS <- mSIS[, -1]
    ISmatrices <- c(ISmatrices, list(SIS = mSIS))
  }
  if (tis) {
    mTIS <- matrix(0, y.n, y.n)
    v1n <- seq(1, y.n)
    loop.indx <- 1:y.n
    tmp <- function(i) {
      mTIS[c(i:y.n), i] <<- v1n[1:c(y.n - i + 1)]
    }
    tmp <- sapply(loop.indx, tmp)
    colnames(mTIS) <- paste("tis", y.index.as.char, sep = "")
    mTIS <- mTIS[, -1]
    ISmatrices <- c(ISmatrices, list(TIS = mTIS))
  }
  if (!is.list(uis) && !identical(as.numeric(uis), 0)) {
    uis <- as.zoo(cbind(uis))
    uis.names <- colnames(uis)
    if (is.null(uis.names)) {
      uis.names <- paste("uisxreg", 1:NCOL(uis), sep = "")
    }
    if (any(uis.names == "")) {
      missing.colnames <- which(uis.names == "")
      for (i in 1:length(missing.colnames)) {
        uis.names[i] <- paste("uisxreg", missing.colnames[i], 
                              sep = "")
      }
    }
    uis <- na.trim(uis, sides = "both", is.na = "any")
    uis.index.as.char <- as.character(index(uis))
    t1 <- which(uis.index.as.char == y.index.as.char[1])
    t2 <- which(uis.index.as.char == y.index.as.char[length(y.index.as.char)])
    uis <- coredata(uis)
    uis <- window(uis, start = t1, end = t2)
    uis <- cbind(coredata(as.zoo(uis)))
    colnames(uis) <- uis.names
    if (nrow(uis) != y.n) 
      stop("nrow(uis) is unequal to no. of observations")
    ISmatrices <- c(ISmatrices, list(UIS = uis))
  }
  if (is.list(uis)) {
    for (i in 1:length(uis)) {
      uis[[i]] <- as.matrix(coredata(as.zoo(uis[[i]])))
      if (nrow(uis[[i]]) != y.n) {
        stop(paste("nrow(uis[[", i, "]]) is unequal to no. of observations", 
                   sep = ""))
      }
    }
    uis.names <- paste("UIS", 1:length(uis), sep = "")
    if (is.null(names(uis))) {
      names(uis) <- uis.names
    } else {
      for (i in 1:length(uis)) {
        if (names(uis)[i] == "") {
          names(uis)[i] <- uis.names[i]
        }else {
          names(uis)[i] <- paste(uis.names[i], ".", names(uis)[i], 
                                 sep = "")
        }
      }
    }
    ISmatrices <- c(ISmatrices, uis)
  }
  if (is.list(blocks)) {
    if (length(ISmatrices) != length(blocks)) {
      stop("No. of IS matrices is unequal to length(blocks)")
    }
    blocks.is.list <- TRUE
    ISblocks <- blocks
  }else {
    blocks.is.list <- FALSE
    ISblocks <- list()
  }
  ISfinalmodels <- list()
  for (i in 1:length(ISmatrices)) {
    if (!blocks.is.list) {
      ncol.adj <- NCOL(ISmatrices[[i]])
      if (is.null(blocks)) {
        blockratio.value <- ncol.adj/(ratio.threshold * 
                                        ncol.adj - mXncol)
        blocksize.value <- ncol.adj/min(y.n * ratio.threshold, 
                                        max.block.size)
        no.of.blocks <- max(2, blockratio.value, blocksize.value)
        no.of.blocks <- ceiling(no.of.blocks)
        no.of.blocks <- min(ncol.adj, no.of.blocks)
      }else{
        no.of.blocks <- blocks
      }
      blocksize <- ceiling(ncol.adj/no.of.blocks)
      partitions.t2 <- blocksize
      for (j in 1:no.of.blocks) {
        if (blocksize * j <= ncol.adj) {
          partitions.t2[j] <- blocksize * j
        }
      }
      if (partitions.t2[length(partitions.t2)] < ncol.adj) {
        partitions.t2 <- c(partitions.t2, ncol.adj)
      }
      blocksadj <- length(partitions.t2)
      partitions.t1 <- partitions.t2 + 1
      partitions.t1 <- c(1, partitions.t1[-blocksadj])
      tmp <- list()
      for (j in 1:blocksadj) {
        tmp[[j]] <- partitions.t1[j]:partitions.t2[j]
      }
      ISblocks[[i]] <- tmp
    }
    #there is some problem with this one
    ISblocksFun <- function(j, i, ISmatrices, ISblocks, mX, 
                            parallel.options, y, olsMethod, t.pval, wald.pval, 
                            do.pet, arLjungB, archLjungB, normality.JarqueB, 
                            user.diagnostics, info.method, mxkeep, include.gum, 
                            include.1cut, include.empty, max.paths, turbo, tol, 
                            LAPACK, max.regs, print.searchinfo) {
      if (length(ISblocks[[i]][[j]]) == 1) {
        tmp <- colnames(ISmatrices[[i]])[ISblocks[[i]][[j]]]
        mXis <- cbind(ISmatrices[[i]][, ISblocks[[i]][[j]]])
        colnames(mXis) <- tmp
        mXis <- cbind(mX, mXis)
      }else {
        mXis <- cbind(mX, ISmatrices[[i]][, ISblocks[[i]][[j]]])
      }
      mXis <- dropvar(mXis, tol = tol, LAPACK = LAPACK, 
                      silent = print.searchinfo)
      if (is.null(parallel.options)) {
        if (print.searchinfo) {
          message("\n", appendLF = FALSE)
          message(names(ISmatrices)[i], " block ", j, 
                  " of ", length(ISblocks[[i]]), ":", appendLF = TRUE)
        }
      }
      #this one performs backwards selection of the best model
      getsis <- getsFun(y, mXis, untransformed.residuals = NULL, 
                        user.estimator = list(name = "kfs.XT", tol = tol, 
                                              LAPACK = LAPACK, method = olsMethod,
                                              in.old.param=param, infomat=mod$infomat,
                                              infomat.std.err=mod$infomat.std.err), 
                        gum.result = NULL, 
                        t.pval = t.pval, wald.pval = wald.pval, do.pet = do.pet, 
                        ar.LjungB = arLjungB, arch.LjungB = archLjungB, 
                        normality.JarqueB = normality.JarqueB, user.diagnostics = user.diagnostics, 
                        gof.function = list(name = "infocrit", method = info.method), 
                        gof.method = "min", keep = mxkeep, include.gum = include.gum, 
                        include.1cut = include.1cut, include.empty = include.empty, 
                        max.paths = max.paths, turbo = turbo, tol = tol, 
                        LAPACK = LAPACK, max.regs = max.regs, print.searchinfo = print.searchinfo, 
                        alarm = FALSE)
      if (is.null(getsis$specific.spec)) {
        ISspecific.models <- NULL
      }else {
        ISspecific.models <- names(getsis$specific.spec)
      }
      return(ISspecific.models)
    }
    if (is.null(parallel.options)) {
      ISspecific.models <- lapply(1:length(ISblocks[[i]]), 
                                  ISblocksFun, i, ISmatrices, ISblocks, mX, parallel.options, 
                                  y, olsMethod, t.pval, wald.pval, do.pet, arLjungB, 
                                  archLjungB, normality.JarqueB, user.diagnostics, 
                                  info.method, mxkeep, include.gum, include.1cut, 
                                  include.empty, max.paths, turbo, tol, LAPACK, 
                                  max.regs, print.searchinfo)
    }
    if (!is.null(parallel.options)) {
      if (print.searchinfo) {
        message("\n", appendLF = FALSE)
        message("Preparing parallel computing...", appendLF = TRUE)
        message(names(ISmatrices)[i], " blocks to search in parallel: ", 
                length(ISblocks[[i]]), appendLF = TRUE)
        message("Searching...", appendLF = TRUE)
      }
      blocksClust <- makeCluster(clusterSpec, outfile = "")
      clusterExport(blocksClust, c("dropvar", "getsFun", 
                                   "ols", "infocrit", "diagnostics"), envir = .GlobalEnv)
      ISspecific.models <- parLapply(blocksClust, 1:length(ISblocks[[i]]), 
                                     ISblocksFun, i, ISmatrices, ISblocks, mX, parallel.options, 
                                     y, olsMethod, t.pval, wald.pval, do.pet, arLjungB, 
                                     archLjungB, normality.JarqueB, user.diagnostics, 
                                     info.method, mxkeep, include.gum, include.1cut, 
                                     include.empty, max.paths, turbo, tol, LAPACK, 
                                     max.regs, print.searchinfo)
      stopCluster(blocksClust)
    }
    if (print.searchinfo) {
      message("\n", appendLF = FALSE)
      message("GETS of union of retained ", names(ISmatrices)[i], 
              " variables... ", appendLF = TRUE)
    }
    if (length(ISspecific.models) == 0) {
      isNames <- NULL
      ISfinalmodels[[i]] <- NULL
    }
    if (length(ISspecific.models) > 0) {
      isNames <- NULL
      for (j in 1:length(ISspecific.models)) {
        if (!is.null(ISspecific.models[[j]])) {
          isNames <- union(isNames, ISspecific.models[[j]])
        }
      }
      isNames <- setdiff(isNames, mXnames)
      if (length(isNames) == 0) {
        ISfinalmodels[[i]] <- mXnames
      }
      else {
        mXisNames <- c(mXnames, isNames)
        mXis <- cbind(mX, ISmatrices[[i]][, isNames])
        colnames(mXis) <- mXisNames
        mXis <- dropvar(mXis, tol = tol, LAPACK = LAPACK, 
                        silent = print.searchinfo)
        getsis <- gets::getsFun(y, mXis, untransformed.residuals = NULL, 
                                list(name = "kfs.XT", tol = tol, 
                                     LAPACK = LAPACK, method = olsMethod,
                                     in.old.param=param, infomat=mod$infomat,
                                     infomat.std.err=mod$infomat.std.err), 
                                gum.result = NULL, 
                          t.pval = t.pval, wald.pval = wald.pval, do.pet = do.pet, 
                          ar.LjungB = arLjungB, arch.LjungB = archLjungB, 
                          normality.JarqueB = normality.JarqueB, user.diagnostics = user.diagnostics, 
                          gof.function = list(name = "infocrit", method = info.method), 
                          gof.method = "min", keep = mxkeep, include.gum = include.gum, 
                          include.1cut = include.1cut, include.empty = include.empty, 
                          max.paths = max.paths, turbo = turbo, tol = tol, 
                          LAPACK = LAPACK, max.regs = max.regs, print.searchinfo = print.searchinfo, 
                          alarm = FALSE)
        ISfinalmodels[[i]] <- names(getsis$specific.spec)
      }
    }
  }
  names(ISblocks) <- names(ISmatrices)
  if (print.searchinfo) {
    message("\n", appendLF = FALSE)
    message("GETS of union of ALL retained variables...", 
            appendLF = TRUE)
    message("\n", appendLF = FALSE)
  }
  if (length(ISfinalmodels) == 0) {
    ISfinalmodels <- NULL
    if (is.null(mX)) {
      mXis <- NULL
    }
    else {
      mXis <- zoo(cbind(mX), order.by = y.index)
      colnames(mXis) <- mXnames
    }
  }
  if (length(ISfinalmodels) > 0) {
    mIS <- NULL
    for (i in 1:length(ISfinalmodels)) {
      isNames <- NULL
      if (!is.null(ISfinalmodels[[i]])) {
        isNames <- setdiff(ISfinalmodels[[i]], mXnames)
      }
      if (length(isNames) > 0) {
        tmp <- cbind(ISmatrices[[i]][, isNames])
        colnames(tmp) <- isNames
        mIS <- cbind(mIS, tmp)
      }
    }
    if (is.null(mIS) && !is.null(mX)){
      mXis <- dropvar(cbind(mX), tol = tol, LAPACK = LAPACK, 
                      silent = print.searchinfo)
    }else if(!is.null(mIS) && is.null(mX)){
      mXis <- dropvar(cbind(mIS), tol = tol, LAPACK = LAPACK, 
                      silent = print.searchinfo)
    }else{ mXis <- dropvar(cbind(mX, mIS), tol = tol, LAPACK = LAPACK, 
                           silent = print.searchinfo)}
    mXis <- zoo(mXis, order.by = y.index)
  }
  #y <- as.ts(y, order.by = y.index)
  #y <- zoo(y, order.by = y.index)
  #mod <- arx(y, mxreg = mXis, vcov.type = vcov.type, qstat.options = qstat.options, 
  #           user.diagnostics = user.diagnostics, tol = tol, LAPACK = LAPACK, 
  #           plot = FALSE)
  #this part has to be done properly
  mod <-kfs.arx(y, mxreg = mXis, vcov.type = vcov.type, qstat.options = qstat.options, 
                user.diagnostics = user.diagnostics, tol = tol, LAPACK = LAPACK 
               )
  y <- zoo(y, order.by = y.index)
  getsis <- getsm_my(mod, keep = mxkeep, t.pval = t.pval, do.pet = do.pet, 
                  wald.pval = wald.pval, ar.LjungB = ar.LjungB, arch.LjungB = arch.LjungB, 
                  normality.JarqueB = normality.JarqueB, user.diagnostics = user.diagnostics, 
                  info.method = info.method, include.empty = include.empty, 
                  max.paths = max.paths, max.regs = max.regs, print.searchinfo = print.searchinfo, 
                  #added by me cuz of using kfs instead of 
                  vcov.type = vcov.type)
  ISnames <- setdiff(getsis$aux$mXnames, mXnames)
  if (length(ISnames) == 0) {
    ISnames <- NULL
  }
  colnames(getsis$aux$mX) <- getsis$aux$mXnames
  getsis$gets.type <- "isat"
  getsis$call <- isat.call
  getsis <- c(list(ISfinalmodels = ISfinalmodels, ISnames = ISnames), 
              getsis)
  getsis$aux$t.pval <- t.pval
  getsis$aux$mX=coredata(getsis$aux$mX)
  getsis$aux$y=coredata(getsis$aux$y)
  class(getsis) <- "isat"
  if (alarm) {
    alarm()
  }
  if (is.null(plot)) {
    plot <- getOption("plot")
    if (is.null(plot)) {
      plot <- FALSE
    }
  }
  if (plot) {
    plot.isat(getsis, coef.path = TRUE)
  }
  getsis$pure=pure
  return(getsis)
}