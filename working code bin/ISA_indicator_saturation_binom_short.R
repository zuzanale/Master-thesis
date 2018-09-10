isat_my_short=function(y, mc = TRUE, ar = NULL, ewma = NULL, mxreg = NULL, 
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
 # if (mod$mean.results$`p-value`>0.05)
  #{
  #  mod$c=0
  #  param= mod$old.param
#  }else{param= t(mod$pars.change)}
  pure$std.eps=mod$std.eps
  pure$std.eta=as.ts(mod$std.state[,1]) #for the level
  pure$std.zeta=as.ts(mod$std.state[,2]) #for the slope
  pure$std.omega=as.ts(mod$std.state[,3]) #for seasonality
  return(list=c(mod,pure))
}
