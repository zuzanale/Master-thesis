getsFun_my = function(y, x, untransformed.residuals = NULL, user.estimator = list(name = "kfs", 
                                                                                tol = 1e-07, LAPACK = FALSE, method = 1), gum.result = NULL, 
                    t.pval = 0.05, wald.pval = t.pval, do.pet = TRUE, ar.LjungB = NULL, 
                    arch.LjungB = NULL, normality.JarqueB = NULL, user.diagnostics = NULL, 
                    gof.function = list(name = "infocrit", method = "sc"), gof.method = c("min", 
                                                                                          "max"), keep = NULL, include.gum = FALSE, include.1cut = FALSE, 
                    include.empty = FALSE, max.paths = NULL, turbo = FALSE, tol = 1e-07, 
                    LAPACK = FALSE, max.regs = NULL, print.searchinfo = TRUE, 
                    alarm = FALSE) 
{########################################
  gof.method <- match.arg(gof.method)
  if (is.null(x) || NCOL(x) == 0) {
    stop("GUM regressor matrix is empty")
  }
  x <- cbind(x)
  aux <- list()
  aux$y.n <- NROW(y)
  aux$xNCOL <- NCOL(x)
  if (is.null(user.estimator$envir)) {
    user.estimator$envir <- .GlobalEnv
  }
  userEstArg <- user.estimator
  userEstArg$name <- NULL
  userEstArg$envir <- NULL
  if (length(userEstArg) == 0) {
    userEstArg <- NULL
  }
  if (is.null(gof.function$envir)) {
    gof.function$envir <- .GlobalEnv
  }
  if (gof.function$name == "infocrit" && is.null(gof.function$method)) {
    gof.function$method <- "sc"
  }
  gofFunArg <- gof.function
  gofFunArg$name <- NULL
  gofFunArg$envir <- NULL
  if (length(gofFunArg) == 0) {
    gofFunArg <- NULL
  }
  if (!is.null(max.paths) && max.paths < 1) {
    stop("'max.paths' cannot be smaller than 1")
  }
  if (!is.null(ar.LjungB) || !is.null(arch.LjungB) || !is.null(normality.JarqueB) || 
      !is.null(user.diagnostics)) {
    doDiagnostics <- TRUE
  }else {
    doDiagnostics <- FALSE
  }
  if (is.null(max.regs)) {
    max.regs <- 10 * aux$y.n
  }
  aux$mR <- matrix(0, aux$xNCOL, aux$xNCOL)
  diag(aux$mR) <- 1
  out <- list()
  out$time.started <- date()
  out$time.finished <- NA
  out$call <- sys.call()
  out$no.of.estimations <- 0
  out$messages <- NULL
  out$paths <- list()
  out$terminals <- list()
  out$terminals.results <- NULL
  row.labels <- NULL
  keep <- as.integer(keep)
  keep.n <- length(keep)
  gum <- 1:aux$xNCOL
  delete <- setdiff(gum, keep)
  delete.n <- length(delete)
  if (is.null(gum.result)) {
    #this part is for regression method used (eg ols), its in userEstArg in the "method" parameter
    #after performace if this part next estimator is tested
    #do.call was ols replaced by mine auxiliary filter
    est=NULL
    estim=list()
    #this should be rewriteen into tidyverse world
    est$coefficients=t(rep(0,ncol(x)))
    est$stats=t(rep(0,ncol(x)))
    est$vcov=t(rep(0,ncol(x)))
    est$df=t(rep(0,ncol(x)))
    est$logl=t(rep(0,ncol(x)))
    est$k=t(rep(0,ncol(x)))
    colnames(est$coefficients)=colnames(x)
    colnames(est$stats)=colnames(x)
    colnames(est$vcov)=colnames(x)
    colnames(est$df)=colnames(x)
    colnames(est$logl)=colnames(x)
    colnames(est$k)=colnames(x)
    for (t in 1:ncol(x)) {
      Xt= x[,t]
      estim[[t]] <- do.call(user.estimator$name, c(list(y = y, x = Xt), 
                                                   userEstArg), envir = user.estimator$envir)
      est$coefficients[1,t]= estim[[t]]$coefficients
      est$stats[1,t]= estim[[t]]$t_stats
      est$vcov[1,t]= estim[[t]]$vcov
      est$df[1,t]= estim[[t]]$df
      est$logl[1,t]= estim[[t]]$logl
      est$k[1,t]= estim[[t]]$k
    }
      out$no.of.estimations <- out$no.of.estimations + 1
  }else {
    est <- gum.result
  }
  if (doDiagnostics) {
    gumDiagnosticsOK <- diagnostics(est, ar.LjungB = ar.LjungB, 
                                    arch.LjungB = arch.LjungB, normality.JarqueB = normality.JarqueB, 
                                    verbose = FALSE, user.fun = user.diagnostics)
  }else {
    gumDiagnosticsOK <- TRUE
  }
  if (gumDiagnosticsOK) {
    gum.regs <- gum
    gum.coefs <- est$coefficients
    gum.varcovmat <- est$vcov #this i have to figure out cuz its needed for wald stats (i have to check on wiki how it is computed)
    stderrs <- sqrt(diag(est$vcov))
    #gum.tstat <- est$coefficients/stderrs
    gum.tstat =  est$stats #this is just for the observation level
    gum.pval <- pt(abs(gum.tstat), est$df, lower.tail = FALSE) *  2
    if (include.gum) {
      out$terminals[[1]] <- gum
      #function like information criteria for variables selection
      gofValue <- do.call(gof.function$name, c(list(x = est), 
                                               gofFunArg), envir = user.estimator$envir)
      out$terminals.results <- rbind(out$terminals.results, 
                                     c(gofValue, est$logl, est$n, est$k))
      row.labels <- c(row.labels, "spec 1 (gum):")
    }
  }else{
    out$messages <- paste(out$messages, "- MGUM does not pass one or more diagnostic checks", 
                          sep = "")
  }
  if (gumDiagnosticsOK && delete.n > 0 && include.1cut) {
    insig.regs <- setdiff(which(gum.pval > t.pval), keep)
    n.paths <- length(insig.regs)
    if (n.paths == 0) {
      out$messages <- paste(out$messages, "- 1-CUT not included (all non-keep regressors are significant)", 
                            sep = "")
    }
    if (n.paths > 0) {
      mXadj <- cbind(x[, -insig.regs])
      estim=list()
      #this should be rewriteen into tidyverse world
      est$coefficients=t(rep(0,ncol(mXadj)))
      est$stats=t(rep(0,ncol(mXadj)))
      est$vcov=t(rep(0,ncol(mXadj)))
      est$df=t(rep(0,ncol(mXadj)))
      est$logl=t(rep(0,ncol(mXadj)))
      est$k=t(rep(0,ncol(mXadj)))
      colnames(est$coefficients)=colnames(mXadj)
      colnames(est$stats)=colnames(mXadj)
      colnames(est$vcov)=colnames(mXadj)
      colnames(est$df)=colnames(mXadj)
      colnames(est$logl)=colnames(mXadj)
      colnames(est$k)=colnames(mXadj)
      for (t in 1:ncol(mXadj)) {
        Xt= mXadj[,t]
        estim[[t]] <- do.call(user.estimator$name, c(list(y = y, x = Xt), 
                                                     userEstArg), envir = user.estimator$envir)
        est$coefficients[1,t]= estim[[t]]$coefficients
        est$stats[1,t]= estim[[t]]$t_stats
        est$vcov[1,t]= estim[[t]]$vcov
        est$df[1,t]= estim[[t]]$df
        est$logl[1,t]= estim[[t]]$logl
        est$k[1,t]= estim[[t]]$k
      }
      
      out$no.of.estimations <- out$no.of.estimations + 
        1
      if (doDiagnostics) {
        diagnosticsOK <- diagnostics(est, ar.LjungB = ar.LjungB, 
                                     arch.LjungB = arch.LjungB, normality.JarqueB = normality.JarqueB, 
                                     verbose = FALSE, user.fun = user.diagnostics)
      }else {
        diagnosticsOK <- TRUE
      }
      if (diagnosticsOK) {
        if (do.pet) {
          mR <- rbind(aux$mR[insig.regs, ])
          mRestq <- mR %*% cbind(gum.coefs)
          wald.stat <- t(mRestq) %*% qr.solve(mR %*% 
                                                gum.varcovmat %*% t(mR), tol = tol) %*% mRestq
          petOK <- as.logical(wald.pval < pchisq(wald.stat, 
                                                 n.paths, lower.tail = FALSE))
        }else {
          petOK <- TRUE
        }
        if (petOK) {
          spec.1cut <- setdiff(gum, insig.regs)
          out$terminals[[length(out$terminals) + 1]] <- spec.1cut
          gofValue <- do.call(gof.function$name, c(list(x = est), 
                                                   gofFunArg), envir = user.estimator$envir)
          out$terminals.results <- rbind(out$terminals.results, 
                                         c(gofValue, est$logl, est$n, est$k))
          row.labels <- c(row.labels, paste("spec ", 
                                            length(out$terminals), " (1-cut):", sep = ""))
        }
      }
    }
  }
  if (gumDiagnosticsOK && delete.n > 0 && include.empty) {
    if (include.1cut && exists("spec.1cut")) {
      emptyEqualTo1cut <- identical(keep, spec.1cut)
    } else {
      emptyEqualTo1cut <- FALSE
    }
    if (emptyEqualTo1cut) {
      out$messages <- paste(out$messages, "- The empty model is equal to the 1-cut model", 
                            sep = "")
    } else {
      mXadj <- cbind(x[, keep])
      
      estim=list()
      #this should be rewriteen into tidyverse world
      est$coefficients=t(rep(0,ncol(mXadj)))
      est$stats=t(rep(0,ncol(mXadj)))
      est$vcov=t(rep(0,ncol(mXadj)))
      est$df=t(rep(0,ncol(mXadj)))
      est$logl=t(rep(0,ncol(mXadj)))
      est$k=t(rep(0,ncol(mXadj)))
      colnames(est$coefficients)=colnames(mXadj)
      colnames(est$stats)=colnames(mXadj)
      colnames(est$vcov)=colnames(mXadj)
      colnames(est$df)=colnames(mXadj)
      colnames(est$logl)=colnames(mXadj)
      colnames(est$k)=colnames(mXadj)
      for (t in 1:ncol(mXadj)) {
        Xt= mXadj[,t]
        estim[[t]] <- do.call(user.estimator$name, c(list(y = y, x = Xt), 
                                                     userEstArg), envir = user.estimator$envir)
        est$coefficients[1,t]= estim[[t]]$coefficients
        est$stats[1,t]= estim[[t]]$t_stats
        est$vcov[1,t]= estim[[t]]$vcov
        est$df[1,t]= estim[[t]]$df
        est$logl[1,t]= estim[[t]]$logl
        est$k[1,t]= estim[[t]]$k
      }
      
      out$no.of.estimations <- out$no.of.estimations + 
        1
      if (doDiagnostics) {
        diagnosticsOK <- diagnostics(est, ar.LjungB = ar.LjungB, 
                                     arch.LjungB = arch.LjungB, normality.JarqueB = normality.JarqueB, 
                                     verbose = FALSE, user.fun = user.diagnostics)
      } else {
        diagnosticsOK <- TRUE
      }
      if (diagnosticsOK) {
        out$terminals[[length(out$terminals) + 1]] <- keep
        gofValue <- do.call(gof.function$name, c(list(x = est), 
                                                 gofFunArg), envir = user.estimator$envir)
        out$terminals.results <- rbind(out$terminals.results, 
                                       c(gofValue, est$logl, est$n, est$k))
        row.labels <- c(row.labels, paste("spec ", length(out$terminals), 
                                          " (empty):", sep = ""))
      } else {
        out$messages <- paste(out$messages, "- Empty model not included (it does not pass one or more diagnostics)", 
                              sep = "")
      }
    }
  }
  insig.regs <- NULL
  pathsTerminals <- list()
  if (gumDiagnosticsOK && delete.n > 1) {
    insig.regs <- setdiff(which(gum.pval > t.pval), keep)
    if (!is.null(max.paths)) {
      if (max.paths < length(insig.regs)) {
        pvalRanksInv <- rank(1 - gum.pval[insig.regs])
        insig.regs <- insig.regs[pvalRanksInv <= max.paths]
      }
    }
    n.paths <- length(insig.regs)
    if (n.paths == 0) {
      out$messages <- paste(out$messages, "- All regressors significant in GUM mean equation", 
                            sep = "")
    }
    if (n.paths > 0) {
      if (print.searchinfo) {
        message(n.paths, " path(s) to search")
        message("Searching: ", appendLF = FALSE)
      }
      regsDeleteList <- list()
      regsKeepList <- list()
      regsMat <- NULL
      for (i in 1:n.paths) {
        if (print.searchinfo) {
          newLine <- ifelse(i == n.paths, TRUE, FALSE)
          message(i, " ", appendLF = newLine)
        }
        path <- insig.regs[i]
        delete.adj <- setdiff(delete, insig.regs[i])
        keep.adj <- keep
        for (j in 1:max.regs) {
          if (turbo && j > 1) {
            regsDeleteList.n <- length(regsDeleteList)
            if (regsDeleteList.n == 0 || i == 1) {
              counter <- regsDeleteList.n + 1
              regsDeleteList[[counter]] <- delete.adj
              regsKeepList[[counter]] <- keep.adj
              regsMat <- rbind(regsMat, c(i, length(path)))
            }else{
              whichOnesInDelete <- which(sapply(regsDeleteList, 
                                                setequal, delete.adj))
              if (length(whichOnesInDelete) == 0) {
                counter <- regsDeleteList.n + 1
                regsDeleteList[[counter]] <- delete.adj
                regsKeepList[[counter]] <- keep.adj
                regsMat <- rbind(regsMat, c(i, length(path)))
                regsDeleteAlreadyDone <- FALSE
              }else{
                regsDeleteAlreadyDone <- TRUE
              }
              if(regsDeleteAlreadyDone) {
                whichOnesInKeep <- which(sapply(regsKeepList, 
                                                setequal, keep.adj))
                whichOne <- intersect(whichOnesInDelete, 
                                      whichOnesInKeep)
                if (length(whichOne) == 1) {
                  regsKeepAlreadyDone <- TRUE
                }else{
                  counter <- regsDeleteList.n + 1
                  regsDeleteList[[counter]] <- delete.adj
                  regsKeepList[[counter]] <- keep.adj
                  regsMat <- rbind(regsMat, c(i, length(path)))
                  regsKeepAlreadyDone <- FALSE
                }
                if (regsKeepAlreadyDone) {
                  spec.adj <- pathsTerminals[[regsMat[whichOne, 
                                                      1]]]
                  pathtmp <- out$paths[[regsMat[whichOne, 
                                                1]]]
                  pathtmp <- pathtmp[-c(1:regsMat[whichOne, 
                                                  2])]
                  path <- c(path, pathtmp)
                  break
                }
              }
            }
          }
          mXadj <- cbind(x[, union(delete.adj, keep.adj)])
      
          estim=list()
          #this should be rewriteen into tidyverse world
          est$coefficients=t(rep(0,ncol(mXadj)))
          est$stats=t(rep(0,ncol(mXadj)))
          est$vcov=t(rep(0,ncol(mXadj)))
          est$df=t(rep(0,ncol(mXadj)))
          est$logl=t(rep(0,ncol(mXadj)))
          est$k=t(rep(0,ncol(mXadj)))
          colnames(est$coefficients)=colnames(mXadj)
          colnames(est$stats)=colnames(mXadj)
          colnames(est$vcov)=colnames(mXadj)
          colnames(est$df)=colnames(mXadj)
          colnames(est$logl)=colnames(mXadj)
          colnames(est$k)=colnames(mXadj)
          for (t in 1:ncol(mXadj)) {
            Xt= mXadj[,t]
            estim[[t]] <- do.call(user.estimator$name, c(list(y = y, x = Xt), 
                                                         userEstArg), envir = user.estimator$envir)
            est$coefficients[1,t]= estim[[t]]$coefficients
            est$stats[1,t]= estim[[t]]$t_stats
            est$vcov[1,t]= estim[[t]]$vcov
            est$df[1,t]= estim[[t]]$df
            est$logl[1,t]= estim[[t]]$logl
            est$k[1,t]= estim[[t]]$k
          }
         
          out$no.of.estimations <- out$no.of.estimations + 
            1
          if (doDiagnostics) {
            diagnosticsOK <- diagnostics(est, ar.LjungB = ar.LjungB, 
                                         arch.LjungB = arch.LjungB, normality.JarqueB = normality.JarqueB, 
                                         verbose = FALSE, user.fun = user.diagnostics)
          }else {
            diagnosticsOK <- TRUE
          }
          if (!diagnosticsOK) {
            path.n <- length(path)
            keep.adj <- union(path[path.n], keep.adj)
            path <- union(path, path[path.n] * c(-1))
            next
          }
          if (diagnosticsOK) {
            if (length(delete.adj) == 0) {
              spec.adj <- sort(keep.adj)
              break
            }
            stderrs <- sqrt(diag(est$vcov))
            t.stat <- est$stats
            p.val <- pt(abs(t.stat), est$df, lower.tail = FALSE) * 
              2 #this is zero for some unknown reason
            if (any(p.val[1:c(length(delete.adj))] > 
                    t.pval) > 0) {
              reg.no <- which.max(p.val[1:c(length(delete.adj))])
              if (do.pet) {
                deleted <- setdiff(delete, delete.adj[-reg.no])
                deleted <- sort(deleted)
                n.deleted <- length(deleted)
                mR <- rbind(aux$mR[deleted, ])
                mRestq <- mR %*% cbind(gum.coefs)
                wald.stat <- t(mRestq) %*% qr.solve(mR %*% 
                                                      gum.varcovmat %*% t(mR), tol = tol) %*% 
                  mRestq
                petOK <- as.logical(wald.pval < pchisq(wald.stat, 
                                                       n.deleted, lower.tail = FALSE))
              }else {
                petOK <- TRUE
              }
              if (petOK) {
                path <- union(path, delete.adj[reg.no])
                delete.adj <- delete.adj[-reg.no]
              }else {
                path <- union(path, delete.adj[reg.no] * 
                                c(-1))
                keep.adj <- union(delete.adj[reg.no], 
                                  keep.adj)
                delete.adj <- delete.adj[-reg.no]
              }
            }else {
              spec.adj <- sort(union(delete.adj, keep.adj))
              break
            }
          }
        }
        counter <- length(out$paths) + 1
        out$paths[[counter]] <- path
        pathsTerminals[[counter]] <- spec.adj
        if (length(out$terminals) == 0) {
          chk.spec <- FALSE
        }
        else {
          for (l in 1:length(out$terminals)) {
            chk.spec <- setequal(spec.adj, out$terminals[[l]])
            if (chk.spec == TRUE) {
              break
            }
          }
        }
        if (chk.spec == FALSE) {
          out$terminals[[length(out$terminals) + 1]] <- spec.adj
          gofValue <- do.call(gof.function$name, c(list(x = est), 
                                                   gofFunArg), envir = user.estimator$envir)
          out$terminals.results <- rbind(out$terminals.results, 
                                         c(gofValue, est$logl, est$n, est$k))
          row.labels <- c(row.labels, paste("spec ", 
                                            length(out$terminals), ":", sep = ""))
        }
      }
    }
  }
  if (!is.null(out$terminals.results)) {
    if (gof.method == "min") {
      out$best.terminal <- which.min(out$terminals.results[, 
                                                           1])
    }
    else {
      out$best.terminal <- which.max(out$terminals.results[, 
                                                           1])
    }
    if (length(out$best.terminal) > 1) {
      out$messages <- paste(out$messages, "- Several 'best' terminals, the first selected", 
                            sep = "")
    }
    out$best.terminal <- out$best.terminal[1]
    out$specific.spec <- out$terminals[[out$best.terminal]]
    if (length(out$specific.spec) == 0) {
      out$specific.spec <- NULL
    }
    else {
      out$specific.spec <- sort(out$specific.spec)
      names(out$specific.spec) <- colnames(x)[out$specific.spec]
    }
    if (gof.function$name == "infocrit") {
      col.labels <- c(paste("info(", gofFunArg$method, 
                            ")", sep = ""), "logl", "n", "k")
    }
    else {
      col.labels <- c("gof-value", "logl", "n", "k")
    }
    if (NCOL(out$terminals.results) != length(col.labels)) {
      col.labels <- c(col.labels[1], rep(NA, NCOL(out$terminals.results) - 
                                           1))
    }
    colnames(out$terminals.results) <- col.labels
    rownames(out$terminals.results) <- row.labels
    if (length(out$paths) == 0) {
      out$paths <- NULL
    }
  }
  out$time.finished <- date()
  if (alarm) {
    alarm()
  }
  return(out)
}