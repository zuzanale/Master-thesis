getsm_my=function (object, t.pval = 0.05, wald.pval = t.pval, vcov.type = NULL, 
                do.pet = TRUE, ar.LjungB = list(lag = NULL, pval = 0.025), 
                arch.LjungB = list(lag = NULL, pval = 0.025), normality.JarqueB = NULL, 
                user.diagnostics = NULL, info.method = c("sc", "aic", "hq"), 
                keep = NULL, include.gum = FALSE, include.empty = FALSE, 
                max.paths = NULL, max.regs = NULL, zero.adj = NULL, vc.adj = NULL, 
                verbose = TRUE, print.searchinfo = TRUE, estimate.specific = TRUE, 
                plot = NULL, alarm = FALSE) 
{
  info.method <- match.arg(info.method)
  if (!is.null(max.paths) && max.paths < 1) {
    stop("'max.paths' cannot be smaller than 1")
  }
  if (is.null(vcov.type)) {
    vcov.type <- object$aux$vcov.type
  }
  else {
    vcovTypes <- c("ordinary", "white", "newey-west")
    which.type <- charmatch(vcov.type, vcovTypes)
    vcov.type <- vcovTypes[which.type]
  }
  if (!is.null(ar.LjungB) && is.null(ar.LjungB$lag)) {
    ar.LjungB$lag <- object$aux$qstat.options[1]
  }
  ar.LjungB <- c(ar.LjungB$lag[1], ar.LjungB$pval[1])
  if (!is.null(arch.LjungB) && is.null(arch.LjungB$lag)) {
    arch.LjungB$lag <- object$aux$qstat.options[2]
  }
  arch.LjungB <- c(arch.LjungB$lag[1], arch.LjungB$pval[1])
  if (is.null(max.regs)) {
    max.regs <- 10 * object$aux$y.n
  }
  tol <- object$aux$tol
  LAPACK <- object$aux$LAPACK
  if (is.null(object$aux$mX)) {
    stop("Mean equation empty")
  }
  if (is.null(object$variance.results)) {
    var.spec.chk <- FALSE
  }
  else {
    var.spec.chk <- TRUE
  }
  if (!is.null(object$aux$user.estimator)) {
    if (estimate.specific) {
      estimate.specific <- FALSE
      message("  'estimate.specific' set to FALSE")
    }
    if (is.null(plot) || identical(plot, TRUE)) {
      plot <- FALSE
      message("  'plot' set to FALSE")
    }
  }
  out <- list()
  out$time.started <- date()
  out$time.finished <- NA
  out$call <- sys.call()
  messages <- NULL
  spec <- list()
  spec.results <- NULL
  keep.n <- length(keep)
  gum <- 1:object$aux$mXncol
  delete <- setdiff(gum, keep)
  delete.n <- length(delete)
  if (delete.n > 0) {
    mXdel <- cbind(object$aux$mX[, delete])
  }
  else {
    mXdel <- NULL
  }
  if (is.null(keep)) {
    mXndel <- NULL
  }
  else {
    mXndel <- cbind(object$aux$mX[, keep])
  }
  mXadj <- cbind(mXdel, mXndel)
  tmp <- rep(0, object$aux$mXncol)
  if (!is.null(keep)) {
    tmp[keep] <- 1
  }
  tmpdf <- cbind(tmp, object$mean.results)
  tmp <- 1:object$aux$mXncol
  tmpdf <- cbind(tmp, tmpdf)
  colnames(tmpdf)[1:2] <- c("reg.no", "keep")
  out$gum.mean <- tmpdf
  out$gum.variance <- object$variance.results
  out$gum.diagnostics <- object$diagnostics
  if (is.null(object$aux$user.estimator)) {
    estMethod <- which(vcov.type == c("none", "none", "ordinary", 
                                      "white", "newey-west"))
    est <- ols(object$aux$y, mXadj, tol = object$aux$tol, 
               LAPACK = object$aux$LAPACK, method = estMethod, user.fun = NULL, 
               user.options = NULL)
    est$std.residuals <- coredata(na.trim(object$std.residuals))
    est$logl <- object$logl
    if (!is.null(object$aux$loge2.n)) {
      est$n <- object$aux$loge2.n
    }
  }
  else {
    est <- do.call(object$aux$user.estimator$name, list(object$aux$y, 
                                                        mXadj), envir = .GlobalEnv)
    if (is.null(est$vcov) && !is.null(est$vcov.mean)) {
      est$vcov <- est$vcov.mean
    }
  }
  if (!is.null(est$residuals)) {
    gum.chk <- diagnostics(est, ar.LjungB = ar.LjungB, arch.LjungB = arch.LjungB, 
                           normality.JarqueB = normality.JarqueB, verbose = FALSE, 
                           user.fun = user.diagnostics)
  }
  else {
    gum.chk <- TRUE
  }
  if (gum.chk) {
    spec[[1]] <- spec.gum <- gum
    stderrs <- sqrt(diag(est$vcov))
    t.stat <- est$coefficients/stderrs
    p.val <- pt(abs(t.stat), est$df, lower.tail = FALSE) * 
      2
    info.results <- info.criterion(est$logl, est$n, est$k, 
                                   method = info.method)
    spec.results <- rbind(c(info.results$value, est$logl, 
                            info.results$n, info.results$k))
    col.labels <- c(paste("info(", info.method, ")", sep = ""), 
                    "logl", "n", "k")
    row.labels <- c("spec 1 (gum):")
    gum.regs <- c(delete, keep)
    gum.coefs <- object$mean.results[gum.regs, 1]
    gum.varcovmat <- est$vcov
  }
  else {
    messages <- paste(messages, "- MGUM does not pass one or more diagnostic checks", 
                      sep = "")
  }
  if (gum.chk && delete.n > 0 && include.empty) {
    if (is.null(object$aux$user.estimator)) {
      est <- ols(object$aux$y, mXndel, tol = object$aux$tol, 
                 LAPACK = object$aux$LAPACK, method = estMethod)
    }
    else {
      est <- do.call(object$aux$user.estimator$name, list(object$aux$y, 
                                                          mXndel), envir = .GlobalEnv)
      if (is.null(est$vcov) && !is.null(est$vcov.mean)) {
        est$vcov <- est$vcov.mean
      }
    }
   # if (var.spec.chk) {
   #   residsAdj <- zoo(est$residuals, order.by = object$aux$y.index)
   #   est.var <- arx(residsAdj, vc = object$aux$vc, arch = object$aux$arch, 
   #                  asym = object$aux$asym, log.ewma = object$aux$log.ewma, 
    #                 vxreg = object$aux$vxreg, zero.adj = object$aux$zero.adj, 
   #                  vc.adj = object$aux$vc.adj, tol = object$aux$tol, 
  #                   LAPACK = object$aux$LAPACK, plot = FALSE)
  #    est$std.residuals <- coredata(na.trim(est.var$std.residuals))
  #   est$var.fit <- coredata(na.trim(est.var$var.fit))
  #    est$logl <- est.var$logl
  #    est$residuals <- est$residuals[c(object$aux$y.n - 
  #                                       object$aux$loge2.n + 1):object$aux$y.n]
  #    est$n <- length(est$residuals)
   # }
    diagnostics.chk <- diagnostics(est, ar.LjungB = ar.LjungB, 
                                   arch.LjungB = arch.LjungB, normality.JarqueB = normality.JarqueB, 
                                   verbose = FALSE, user.fun = user.diagnostics)
    if (diagnostics.chk) {
      spec[[length(spec) + 1]] <- if (is.null(keep)) {
        0
      }
      else {
        keep
      }
      info.results <- info.criterion(est$logl, est$n, est$k, 
                                     method = info.method)
      spec.results <- rbind(spec.results, c(info.results$value, 
                                            est$logl, info.results$n, info.results$k))
      row.labels <- c(row.labels, paste("spec ", length(spec), 
                                        " (empty):", sep = ""))
    }
    else {
      messages <- paste(messages, "- Empty mean model does not pass one or more diagnostic checks", 
                        sep = "")
    }
  }
  insig.regs <- NULL
  paths <- list()
  if (gum.chk && delete.n > 1) {
    insig.regs <- delete[which(p.val[1:delete.n] > t.pval)]
    if (!is.null(max.paths)) {
      if (max.paths < length(insig.regs)) {
        pvalRanksInv <- rank(1 - p.val[insig.regs])
        insig.regs <- delete[pvalRanksInv <= max.paths]
      }
    }
    n.paths <- length(insig.regs)
    if (n.paths == 0) {
      messages <- paste(messages, "- All regressors significant in GUM mean equation", 
                        sep = "")
    }
    if (n.paths > 0) {
      if (print.searchinfo) {
        message(n.paths, " path(s) to search")
        message("Searching: ", appendLF = FALSE)
      }
      for (i in 1:n.paths) {
        if (print.searchinfo) {
          newLine <- ifelse(i == n.paths, TRUE, FALSE)
          message(i, " ", appendLF = newLine)
        }
        path <- insig.regs[i]
        delete.adj <- setdiff(delete, insig.regs[i])
        keep.adj <- as.numeric(keep)
        for (j in 1:max.regs) {
          mXdell <- if (length(delete.adj) == 0) {
            NULL
          }
          else {
            object$aux$mX[, delete.adj]
          }
          mXndell <- if (is.null(keep.adj)) {
            NULL
          }
          else {
            object$aux$mX[, keep.adj]
          }
          mXadj <- cbind(mXdell, mXndell)
          mXadj.k <- NCOL(mXadj)
          if (is.null(object$aux$user.estimator)) {
            est <- ols(object$aux$y, mXadj, tol = object$aux$tol, 
                       LAPACK = object$aux$LAPACK, method = estMethod, 
                       user.fun = NULL, user.options = NULL)
          }
          else {
            est <- do.call(object$aux$user.estimator$name, 
                           list(object$aux$y, mXadj), envir = .GlobalEnv)
            if (is.null(est$vcov) && !is.null(est$vcov.mean)) {
              est$vcov <- est$vcov.mean
            }
          }
         # if (var.spec.chk) {
        #    residsAdj <- zoo(est$residuals, order.by = object$aux$y.index)
         #   est.var <- arx(residsAdj, vc = object$aux$vc, 
        #                   arch = object$aux$arch, asym = object$aux$asym, 
       #                    log.ewma = object$aux$log.ewma, vxreg = object$aux$vxreg, 
        #                   zero.adj = object$aux$zero.adj, vc.adj = object$aux$vc.adj, 
       #                    tol = object$aux$tol, LAPACK = object$aux$LAPACK, 
       #                    plot = FALSE)
        #    est$std.residuals <- coredata(na.trim(est.var$std.residuals))
       #     est$var.fit <- coredata(na.trim(est.var$var.fit))
       #     est$logl <- est.var$logl
      #      est$residuals <- est$residuals[c(object$aux$y.n - 
      #                                         object$aux$loge2.n + 1):object$aux$y.n]
      #      est$n <- length(est$residuals)
       #   }
          diagnostics.chk <- diagnostics(est, ar.LjungB = ar.LjungB, 
                                         arch.LjungB = arch.LjungB, normality.JarqueB = normality.JarqueB, 
                                         verbose = FALSE, user.fun = user.diagnostics)
          if (!diagnostics.chk) {
            path.n <- length(path)
            keep.adj <- union(path[path.n], keep.adj)
            path <- union(path, path[path.n] * c(-1))
            next
          }
          if (diagnostics.chk) {
            if (length(delete.adj) == 0) {
              spec.adj <- keep.adj
              break
            }
            stderrs <- sqrt(diag(est$vcov))
            t.stat <- est$t.stats
            p.val <- pt(abs(t.stat), est$df, lower.tail = FALSE) * 
              2
            if (sum(p.val[1:c(length(delete.adj))] > 
                    t.pval) > 0) {
              reg.no <- which.max(p.val[1:I(length(delete.adj))])
              if (do.pet) {
                deleted <- setdiff(delete, delete.adj[-reg.no])
                n.deleted <- length(deleted)
                mR <- NULL
                for (k in 1:length(gum.regs)) {
                  if (gum.regs[k] %in% deleted) {
                    mR <- rbind(mR, c(rep(0, I(k - 1)), 
                                      1, rep(0, I(object$aux$mXncol - 
                                                    k))))
                  }
                }
                mRestq <- mR %*% cbind(gum.coefs)
                wald.stat <- t(mRestq) %*% qr.solve(mR %*% 
                                                      gum.varcovmat %*% t(mR), tol = object$aux$tol) %*% 
                  mRestq
                pet.chk <- as.logical(wald.pval < pchisq(wald.stat, 
                                                         n.deleted, lower.tail = FALSE))
              }
              else {
                pet.chk <- TRUE
              }
              if (pet.chk) {
                path <- union(path, delete.adj[reg.no])
                delete.adj <- delete.adj[-reg.no]
              }
              else {
                path <- union(path, delete.adj[reg.no] * 
                                I(-1))
                keep.adj <- union(delete.adj[reg.no], 
                                  keep.adj)
                delete.adj <- delete.adj[-reg.no]
              }
            }
            else {
              spec.adj <- union(delete.adj, keep.adj)
              break
            }
          }
        }
        paths[[length(paths) + 1]] <- path
        if (length(spec.adj) == 0) {
          spec.adj <- 0
        }
        for (l in 1:length(spec)) {
          chk.spec <- setequal(spec.adj, spec[[l]])
          if (chk.spec == TRUE) {
            break
          }
        }
        if (chk.spec == FALSE) {
          spec[[length(spec) + 1]] <- spec.adj
          if (spec.adj[1] == 0) {
            n.spec.adj <- 0
          }
          else {
            n.spec.adj <- length(spec.adj)
          }
          info.results <- info.criterion(est$logl, est$n, 
                                         n.spec.adj, method = info.method)
          spec.results <- rbind(spec.results, c(info.results$value, 
                                                est$logl, info.results$n, info.results$k))
          row.labels <- c(row.labels, paste("spec ", 
                                            length(spec), ":", sep = ""))
        }
      }
    }
  }
  if (!is.null(spec.results)) {
    J <- 1:NROW(spec.results)
    models <- cbind(J, spec.results)
    colnames(models) <- NULL
    if (include.gum) {
      min.value <- min(models[, 2])
      where <- which(min.value == models[, 2])
    }
    else {
      if (length(spec) == 1) {
        where <- 1
      }
      else {
        min.value <- min(models[-1, 2])
        where <- which(min.value == models[-1, 2]) + 
          1
      }
    }
    if (length(where) > 1) {
      messages <- paste(messages, "- Several terminal specifications attain the minimum information criterion", 
                        sep = "")
    }
    best.spec <- spec[[where[1]]]
  }
  out$keep <- keep
  if (is.null(spec.results)) {
    out$aux <- object$aux
    out$aux$vcov.type <- vcov.type
  }
  if (!is.null(spec.results)) {
    if (length(paths) == 0) {
      out$paths <- NULL
    }
    else {
      out$paths <- paths
    }
    out$terminals <- spec
    colnames(spec.results) <- col.labels
    where.empty <- which(spec.results[, "k"] == 0)
    if (include.empty == FALSE && length(where.empty) > 0) {
      row.labels[where.empty] <- paste("spec ", where.empty, 
                                       " (empty):", sep = "")
    }
    rownames(spec.results) <- row.labels
    out$terminals.results <- spec.results
    if (!estimate.specific) {
      if (best.spec == 0 || is.na(best.spec) || length(best.spec) == 
          0) {
        out$specific.spec <- NULL
      }
      else {
        specific <- sort(best.spec)
        names(specific) <- object$aux$mXnames[specific]
        out$specific.spec <- specific
      }
    }
    if (estimate.specific) {
      yadj <- zoo(cbind(object$aux$y), order.by = object$aux$y.index)
      colnames(yadj) <- object$aux$y.name
      specific <- sort(best.spec)
      if (specific[1] == 0) {
        mXadj <- NULL
      }
      else {
        mXadj <- cbind(object$aux$mX[, specific])
        colnames(mXadj) <- object$aux$mXnames[specific]
        mXadj <- zoo(mXadj, order.by = object$aux$y.index)
      }
      if (is.null(object$aux$vxreg)) {
        vxregAdj <- NULL
      }
      else {
        vxregAdj <- zoo(object$aux$vxreg, order.by = object$aux$y.index)
      }
      if (is.null(ar.LjungB)) {
        ar.LjungB <- object$aux$qstat.options[1]
      }
      if (is.null(arch.LjungB)) {
        arch.LjungB <- object$aux$qstat.options[2]
      }
      est <- kfs.arx(yadj, mxreg = mXadj,
                  asym = object$aux$asym, 
                vxreg = vxregAdj, 
                 zero.adj = object$aux$zero.adj, 
                 vcov.type = vcov.type, qstat.options = c(ar.LjungB[1], 
                                                          arch.LjungB[1]), user.diagnostics = user.diagnostics, 
                 tol = object$aux$tol, LAPACK = object$aux$LAPACK, 
                 plot = FALSE)
      est$call <- est$date <- NULL
      where.diagnostics <- which(names(est) == "diagnostics")
      if (length(where.diagnostics) > 0) {
        names(est)[where.diagnostics] <- "specific.diagnostics"
      }
      est$aux$y.name <- object$aux$y.name
      est <- unclass(est)
      names(specific) <- colnames(mXadj)
      out$specific.spec <- specific
      out <- c(out, est)
    }
  }
  if (!is.null(messages)) {
    out$messages <- messages
  }
  if (is.null(object$aux$user.estimator)) {
    out$aux$y <- object$aux$y
    out$aux$y.index <- object$aux$y.index
    out$aux$y.n <- object$aux$y.n
    out$aux$y.name <- object$aux$y.name
    out$aux$mXnames.gum <- object$aux$mXnames
    out$aux$call.gum <- object$call
    if (is.null(out$aux$vcov.type)) {
      out$aux$vcov.type <- vcov.type
    }
  }
  out <- c(list(date = date(), gets.type = "getsm"), out)
  out$time.finished <- date()
  class(out) <- "gets"
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
    plot.gets(out)
  }
  return(out)
}