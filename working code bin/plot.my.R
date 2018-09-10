plot.my=function (x, col = c("gray", "black"), lty = c(1, 2), 
                  lwd = c(1, 2), coef.path = TRUE, ...) 
{
  if (!is.null(x$mean.fit)) {
    if (length(lwd) == 1) {
      print("lwd needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lwd = rep(lwd, 2)
    }
    else if (length(lwd) > 2) {
      print("lwd needs two arguments, but more provided. First two used.")
      lwd = lwd[1:2]
    }
    if (length(lty) == 1) {
      print("lty needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lty = rep(lty, 2)
    }
    else if (length(lwd) > 2) {
      print("lty needs two arguments, but more provided. First two used.")
      lty = lty[1:2]
    }
    if (length(col) != 2) {
      randomcol <- function() {
        r.r <- runif(2)
        while (round(r.r[1], 1) == round(r.r[2], 1)) {
          r.r <- runif(2)
        }
        g.r <- runif(2)
        while (round(g.r[1], 1) == round(g.r[2], 1)) {
          g.r <- runif(2)
        }
        b.r <- runif(2)
        while (round(b.r[1], 1) == round(b.r[2], 1)) {
          b.r <- runif(2)
        }
        col <- rgb(runif(2), runif(2), runif(2))
        return(col)
      }
      clashcol <- function() {
        r.1 <- runif(1)
        g.1 <- runif(1)
        b.1 <- runif(1)
        b.2 <- min(r.1, min(g.1, b.1)) + max(r.1, max(g.1, 
                                                      b.1))
        col <- rgb(c(r.1, b.2 - r.1), c(g.1, b.2 - g.1), 
                   c(b.1, b.2 - b.1))
        return(col)
      }
      if (col[1] == "random") {
        col <- randomcol()
      }
      else if (col[1] == "awful.clash") {
        col <- clashcol()
      }
      else {
        print("Wrong number of colours specified; using random set of colours instead.")
        col <- randomcol()
      }
    }
    fitted <- x$mean.fit
    actual <- zoo(x$aux$y, order.by = x$aux$y.index)
    residuals <- x$std.residuals
    actual.name <- x$aux$y.name
    def.par <- par(no.readonly = TRUE)
    par(mar = c(2, 2, 0.5, 0.5))
    if ((x$gets.type == "isat" || coef.path == TRUE) && length(x$ISnames) != 
        0) {
      par(mfrow = c(3, 1))
      is.x <- cbind(x$aux$mX[, x$aux$mXnames %in% x$ISnames])
      is.coef.ests <- coef.isat(x)[x$ISnames]
      coef.path.0 <- zoo(is.x %*% is.coef.ests, order.by = x$aux$y.index)
    }
    else {
      par(mfrow = c(2, 1))
    }
    if (is.regular(actual)) {
      plot(actual, main = "", ylim = range(min(actual, 
                                               fitted, na.rm = TRUE), max(actual, fitted, na.rm = TRUE)), 
           type = "l", ylab = "", xlab = "", col = col[2])
      
    }
    else {
      plot(as.Date(index(actual)), coredata(actual), main = "", 
           ylim = range(min(actual, fitted, na.rm = TRUE), 
                        max(actual, fitted, na.rm = TRUE)), type = "l", 
           ylab = "", xlab = "", col = col[2])
    }
    if (is.regular(fitted)) {
      lines(fitted, col = col[1],lty=2,lwd=2)
    }
    else {
      lines(as.Date(index(fitted)), coredata(fitted), col = col[1],lty=2,lwd=2)
    }
    legend("topleft", lty = lty, lwd = lwd, ncol = 2, col = col[c(2, 
                                                                  1)], legend = c(actual.name, "fitted"), bty = "n")
    if (is.regular(residuals)) {
      plot(residuals, type = "h", col = col[1])
    }
    else {
      plot(as.Date(index(residuals)), coredata(residuals), 
           type = "h", col = col[1])
     
    }
    abline(0, 0)
    legend("topleft", lty = 1, col = col[1], legend = "standardised residuals", 
           bty = "n")
    if ((x$gets.type == "isat" | coef.path == TRUE) & length(x$ISnames) != 
        0) {
      if (!is.null(as.list(x$call)$tis) && as.list(x$call)$tis == 
          TRUE) {
        message("\n", appendLF = FALSE)
        message("NB: Because TIS selected, coefficient standard errors invalid hence not plotted", 
                appendLF = TRUE)
        ylim.values <- range(coef.path.0)
        if (is.regular(coef.path.0)) {
          ylim.values <- range(coef.path.0)
          plot(coef.path.0, type = "l", col = col[1], 
               ylim = ylim.values)
        }
        else {
          plot(as.Date(index(coef.path.0)), coredata(coef.path.0), 
               type = "l", col = col[1], ylim = ylim.values)
        }
      }
      else {
        coef.path.v <- isatvar(x)
        if (is.regular(coef.path.0)) {
          ylim.values <- range(min(coef.path.0 - qt(0.975, 
                                                    NROW(coef.path.0)) * coef.path.v$const.se), 
                               max(coef.path.0 + qt(0.975, NROW(coef.path.0)) * 
                                     coef.path.v$const.se))
          plot(coef.path.0, type = "l", col = col[1], 
               ylim = ylim.values)
          lines(coef.path.0 + qt(0.975, NROW(coef.path.0)) * 
                  coef.path.v$const.se, type = "l", col = col[1], 
                lty = 3)
          lines(coef.path.0 - qt(0.975, NROW(coef.path.0)) * 
                  coef.path.v$const.se, type = "l", col = col[1], 
                lty = 3)
        }
        else {
          plot(as.Date(index(coef.path.0)), coredata(coef.path.0), 
               type = "l", col = col[1], ylim = ylim.values)
          lines(as.Date(index(coef.path.0)), coredata(coef.path.0) + 
                  qt(0.975, NROW(coef.path.0)) * coef.path.v$const.se, 
                type = "l", col = col[1], lty = 3)
          lines(as.Date(index(coef.path.0)), coredata(coef.path.0) - 
                  qt(0.975, NROW(coef.path.0)) * coef.path.v$const.se, 
                type = "l", col = col[1], lty = 3)
        }
      }
      abline(0, 0, lty = 3)
      legend("topleft", lty = 1, col = col[1], legend = c(paste(actual.name, 
                                                                "Coefficient Path", sep = ": ")), bty = "n")
    }
    par(def.par)
  }
}

plot.my2=function (x, col = c("gray", "black"), lty = c(1, 2), 
                   lwd = c(1, 2), coef.path = TRUE, ...) 
{
  if (!is.null(x$mean.fit)) {
    if (length(lwd) == 1) {
      print("lwd needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lwd = rep(lwd, 2)
    }
    else if (length(lwd) > 2) {
      print("lwd needs two arguments, but more provided. First two used.")
      lwd = lwd[1:2]
    }
    if (length(lty) == 1) {
      print("lty needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lty = rep(lty, 2)
    }
    else if (length(lwd) > 2) {
      print("lty needs two arguments, but more provided. First two used.")
      lty = lty[1:2]
    }
    if (length(col) != 2) {
      randomcol <- function() {
        r.r <- runif(2)
        while (round(r.r[1], 1) == round(r.r[2], 1)) {
          r.r <- runif(2)
        }
        g.r <- runif(2)
        while (round(g.r[1], 1) == round(g.r[2], 1)) {
          g.r <- runif(2)
        }
        b.r <- runif(2)
        while (round(b.r[1], 1) == round(b.r[2], 1)) {
          b.r <- runif(2)
        }
        col <- rgb(runif(2), runif(2), runif(2))
        return(col)
      }
      clashcol <- function() {
        r.1 <- runif(1)
        g.1 <- runif(1)
        b.1 <- runif(1)
        b.2 <- min(r.1, min(g.1, b.1)) + max(r.1, max(g.1, 
                                                      b.1))
        col <- rgb(c(r.1, b.2 - r.1), c(g.1, b.2 - g.1), 
                   c(b.1, b.2 - b.1))
        return(col)
      }
      if (col[1] == "random") {
        col <- randomcol()
      }
      else if (col[1] == "awful.clash") {
        col <- clashcol()
      }
      else {
        print("Wrong number of colours specified; using random set of colours instead.")
        col <- randomcol()
      }
    }
    fitted <- x$mean.fit
    index(fitted)=index(fitted)*7
    actual <- zoo(x$aux$y, order.by = x$aux$y.index)
    residuals <- x$std.residuals
    actual.name <- x$aux$y.name
    index(actual)=index(actual)*7
    index(residuals)=index(residuals)*7
    def.par <- par(no.readonly = TRUE)
    par(mar = c(2, 2, 0.5, 0.5))
    if ((x$gets.type == "isat" || coef.path == TRUE) && length(x$ISnames) != 
        0) {
      par(mfrow = c(3, 1))
      is.x <- cbind(x$aux$mX[, x$aux$mXnames %in% x$ISnames])
      is.coef.ests <- coef.isat(x)[x$ISnames]
      coef.path.0 <- zoo(is.x %*% is.coef.ests, order.by = x$aux$y.index)
      index(coef.path.0)=index(coef.path.0)*7
    }
    else {
      par(mfrow = c(2, 1))
    }
   # time.index=(index(fitted)*7)
   # xticks=seq(time.index[1],time.index[length(time.index)],10)
  #  org.ticks<- seq(index(fitted)[1],index(fitted)[length(time.index)],10/7)
    
    if (is.regular(actual)) {
      plot(actual, main = "", ylim = range(min(actual, 
                                               fitted, na.rm = TRUE), max(actual, fitted, na.rm = TRUE)+0.5), 
           type = "l", ylab = "", xlab = "", col = col[2])
      
      #axis(1, at=org.ticks,labels = xticks)
    }
    else {
      plot(as.Date(index(actual)), coredata(actual), main = "", 
           ylim = range(min(actual, fitted, na.rm = TRUE), 
                        max(actual, fitted, na.rm = TRUE)+0.5), type = "l", 
           ylab = "", xlab = "", col = col[2])
      #axis(1,at=org.ticks,labels = xticks)
    }
    if (is.regular(fitted)) {
      lines(fitted, col = col[1],lty=2,lwd=2)
     # axis(1, at=org.ticks,labels = xticks)
    }
    else {
      lines(as.Date(index(fitted)), coredata(fitted), col = col[1],lty=2,lwd=2)
    #  axis(1, at=org.ticks,labels = xticks)
    }
    legend("topleft", lty = lty, lwd = lwd, ncol = 2, col = col[c(2, 
                                                                  1)], legend = c(actual.name, "fitted"), bty = "n")
    if (is.regular(residuals)) {
      plot(residuals, type = "h", col = col[1])
     # axis(1, at=org.ticks,labels = xticks)
    }
    else {
      plot(as.Date(index(residuals)), coredata(residuals), 
           type = "h", col = col[1])
    #  axis(1, at=org.ticks,labels = xticks)
    }
    abline(0, 0)
    legend("topleft", lty = 1, col = col[1], legend = "standardised residuals", 
           bty = "n")
    if ((x$gets.type == "isat" | coef.path == TRUE) & length(x$ISnames) != 
        0) {
      if (!is.null(as.list(x$call)$tis) && as.list(x$call)$tis == 
          TRUE) {
        message("\n", appendLF = FALSE)
        message("NB: Because TIS selected, coefficient standard errors invalid hence not plotted", 
                appendLF = TRUE)
        ylim.values <- range(coef.path.0)*7
        if (is.regular(coef.path.0)) {
          ylim.values <- range(coef.path.0)*7
          plot(coef.path.0, type = "l", col = col[1], 
               ylim = ylim.values)
        }
        else {
          plot(as.Date(index(coef.path.0)), coredata(coef.path.0), 
               type = "l", col = col[1], ylim = ylim.values)
        }
      }
      else {
        coef.path.v <- isatvar(x)
        if (is.regular(coef.path.0)) {
          ylim.values <- range(min(coef.path.0 - qt(0.975, 
                                                    NROW(coef.path.0)) * coef.path.v$const.se)+10, 
                               max(coef.path.0 + qt(0.975, NROW(coef.path.0)-10) * 
                                     coef.path.v$const.se))
          plot(coef.path.0, type = "l", col = col[1], 
               ylim = ylim.values)
          lines(coef.path.0 + qt(0.975, NROW(coef.path.0)) * 
                  coef.path.v$const.se, type = "l", col = col[1], 
                lty = 3)
          lines(coef.path.0 - qt(0.975, NROW(coef.path.0)) * 
                  coef.path.v$const.se, type = "l", col = col[1], 
                lty = 3)
        }
        else {
          plot(as.Date(index(coef.path.0)), coredata(coef.path.0), 
               type = "l", col = col[1], ylim = ylim.values)
          lines(as.Date(index(coef.path.0)), coredata(coef.path.0) + 
                  qt(0.975, NROW(coef.path.0)) * coef.path.v$const.se, 
                type = "l", col = col[1], lty = 3)
          lines(as.Date(index(coef.path.0)), coredata(coef.path.0) - 
                  qt(0.975, NROW(coef.path.0)) * coef.path.v$const.se, 
                type = "l", col = col[1], lty = 3)
        }
      }
      abline(0, 0, lty = 3)
      legend("topleft", lty = 1, col = col[1], legend = c(paste(actual.name, 
                                                                "Coefficient Path", sep = ": ")), bty = "n")
    }
    par(def.par)
  }
}

plot.my3=function (x, col = c("gray", "black"), lty = c(1, 2), 
                   lwd = c(1, 2), coef.path = TRUE, ...) 
{
  if (!is.null(x$mean.fit)) {
    if (length(lwd) == 1) {
      print("lwd needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lwd = rep(lwd, 2)
    }
    else if (length(lwd) > 2) {
      print("lwd needs two arguments, but more provided. First two used.")
      lwd = lwd[1:2]
    }
    if (length(lty) == 1) {
      print("lty needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lty = rep(lty, 2)
    }
    else if (length(lwd) > 2) {
      print("lty needs two arguments, but more provided. First two used.")
      lty = lty[1:2]
    }
    if (length(col) != 2) {
      randomcol <- function() {
        r.r <- runif(2)
        while (round(r.r[1], 1) == round(r.r[2], 1)) {
          r.r <- runif(2)
        }
        g.r <- runif(2)
        while (round(g.r[1], 1) == round(g.r[2], 1)) {
          g.r <- runif(2)
        }
        b.r <- runif(2)
        while (round(b.r[1], 1) == round(b.r[2], 1)) {
          b.r <- runif(2)
        }
        col <- rgb(runif(2), runif(2), runif(2))
        return(col)
      }
      clashcol <- function() {
        r.1 <- runif(1)
        g.1 <- runif(1)
        b.1 <- runif(1)
        b.2 <- min(r.1, min(g.1, b.1)) + max(r.1, max(g.1, 
                                                      b.1))
        col <- rgb(c(r.1, b.2 - r.1), c(g.1, b.2 - g.1), 
                   c(b.1, b.2 - b.1))
        return(col)
      }
      if (col[1] == "random") {
        col <- randomcol()
      }
      else if (col[1] == "awful.clash") {
        col <- clashcol()
      }
      else {
        print("Wrong number of colours specified; using random set of colours instead.")
        col <- randomcol()
      }
    }
    fitted <- x$mean.fit
    index(fitted)=index(fitted)*7
    actual <- zoo(x$aux$y, order.by = x$aux$y.index)
    residuals <- x$std.residuals
    actual.name <- x$aux$y.name
    index(actual)=index(actual)*7
    index(residuals)=index(residuals)*7
    def.par <- par(no.readonly = TRUE)
    par(mar = c(2, 2, 0.5, 0.5))
    if ((x$gets.type == "isat" || coef.path == TRUE) && length(x$ISnames) != 
        0) {
      par(mfrow = c(3, 1))
      is.x <- cbind(x$aux$mX[, x$aux$mXnames %in% x$ISnames])
      is.coef.ests <- coef.isat(x)[x$ISnames]
      coef.path.0 <- zoo(is.x %*% is.coef.ests, order.by = x$aux$y.index)
      index(coef.path.0)=index(coef.path.0)*7
    }
    else {
      par(mfrow = c(2, 1))
    }
    # time.index=(index(fitted)*7)
    # xticks=seq(time.index[1],time.index[length(time.index)],10)
    #  org.ticks<- seq(index(fitted)[1],index(fitted)[length(time.index)],10/7)
    
    if (is.regular(actual)) {
      plot(actual, main = "", ylim = range(min(actual, 
                                               fitted, na.rm = TRUE), max(actual, fitted, na.rm = TRUE)+0.5), 
           type = "l", ylab = "", xlab = "", col = col[2])
      
      #axis(1, at=org.ticks,labels = xticks)
    }
    else {
      plot(as.Date(index(actual)), coredata(actual), main = "", 
           ylim = range(min(actual, fitted, na.rm = TRUE), 
                        max(actual, fitted, na.rm = TRUE)+0.5), type = "l", 
           ylab = "", xlab = "", col = col[2])
      #axis(1,at=org.ticks,labels = xticks)
    }
    if(is.regular(fitted)) {
      lines(fitted, col = col[1],lwd=2)
      # axis(1, at=org.ticks,labels = xticks)
    }
    else{
      lines(as.Date(index(fitted)), coredata(fitted), col = col[1],lwd=2)
      #  axis(1, at=org.ticks,labels = xticks)
    }
    legend("topleft", lty = lty, lwd = lwd, ncol = 2, col = col[c(2, 
                                                                  1)], legend = c(actual.name, "fitted"), bty = "n")
    if (is.regular(residuals)) {
      plot(residuals, type = "h", col = col[1])
      # axis(1, at=org.ticks,labels = xticks)
    }
    else {
      plot(as.Date(index(residuals)), coredata(residuals), 
           type = "h", col = col[1])
      #  axis(1, at=org.ticks,labels = xticks)
    }
    abline(0, 0)
    legend("topleft", lty = 1, col = col[1], legend = "standardised residuals", 
           bty = "n")
    if ((x$gets.type == "isat" | coef.path == TRUE) & length(x$ISnames) != 
        0) {
      if (!is.null(as.list(x$call)$tis) && as.list(x$call)$tis == 
          TRUE) {
        message("\n", appendLF = FALSE)
        message("NB: Because TIS selected, coefficient standard errors invalid hence not plotted", 
                appendLF = TRUE)
        ylim.values <- range(coef.path.0)*7
        if (is.regular(coef.path.0)) {
          ylim.values <- range(coef.path.0)*7
          plot(coef.path.0, type = "l", col = col[1], 
               ylim = ylim.values)
        }
        else {
          plot(as.Date(index(coef.path.0)), coredata(coef.path.0), 
               type = "l", col = col[1], ylim = ylim.values)
        }
      }
      else {
        coef.path.v <- isatvar(x)
        if (is.regular(coef.path.0)) {
          ylim.values <- range(min(coef.path.0 - qt(0.975, 
                                                    NROW(coef.path.0)) * coef.path.v$const.se), 
                               max(coef.path.0 + qt(0.975, NROW(coef.path.0)-40) * 
                                     coef.path.v$const.se))
          plot(coef.path.0, type = "l", col = col[1], 
               ylim = ylim.values)
          lines(coef.path.0 + qt(0.975, NROW(coef.path.0)) * 
                  coef.path.v$const.se, type = "l", col = col[1], 
                lty = 3)
          lines(coef.path.0 - qt(0.975, NROW(coef.path.0)) * 
                  coef.path.v$const.se, type = "l", col = col[1], 
                lty = 3)
        }
        else {
          plot(as.Date(index(coef.path.0)), coredata(coef.path.0), 
               type = "l", col = col[1], ylim = ylim.values)
          lines(as.Date(index(coef.path.0)), coredata(coef.path.0) + 
                  qt(0.975, NROW(coef.path.0)) * coef.path.v$const.se, 
                type = "l", col = col[1], lty = 3)
          lines(as.Date(index(coef.path.0)), coredata(coef.path.0) - 
                  qt(0.975, NROW(coef.path.0)) * coef.path.v$const.se, 
                type = "l", col = col[1], lty = 3)
        }
      }
      abline(0, 0, lty = 3)
      legend("topleft", lty = 1, col = col[1], legend = c(paste(actual.name, 
                                                                "Coefficient Path", sep = ": ")), bty = "n")
    }
    par(def.par)
  }
}

plot.my4=function (x, col = c("gray", "black"), lty = c(1, 2), 
                   lwd = c(1, 2), coef.path = TRUE, ...) 
{
  if (!is.null(x$mean.fit)) {
    if (length(lwd) == 1) {
      print("lwd needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lwd = rep(lwd, 2)
    }
    else if (length(lwd) > 2) {
      print("lwd needs two arguments, but more provided. First two used.")
      lwd = lwd[1:2]
    }
    if (length(lty) == 1) {
      print("lty needs two arguments, but only one provided. Single argument applied to all lines plotted.")
      lty = rep(lty, 2)
    }
    else if (length(lwd) > 2) {
      print("lty needs two arguments, but more provided. First two used.")
      lty = lty[1:2]
    }
    if (length(col) != 2) {
      randomcol <- function() {
        r.r <- runif(2)
        while (round(r.r[1], 1) == round(r.r[2], 1)) {
          r.r <- runif(2)
        }
        g.r <- runif(2)
        while (round(g.r[1], 1) == round(g.r[2], 1)) {
          g.r <- runif(2)
        }
        b.r <- runif(2)
        while (round(b.r[1], 1) == round(b.r[2], 1)) {
          b.r <- runif(2)
        }
        col <- rgb(runif(2), runif(2), runif(2))
        return(col)
      }
      clashcol <- function() {
        r.1 <- runif(1)
        g.1 <- runif(1)
        b.1 <- runif(1)
        b.2 <- min(r.1, min(g.1, b.1)) + max(r.1, max(g.1, 
                                                      b.1))
        col <- rgb(c(r.1, b.2 - r.1), c(g.1, b.2 - g.1), 
                   c(b.1, b.2 - b.1))
        return(col)
      }
      if (col[1] == "random") {
        col <- randomcol()
      }
      else if (col[1] == "awful.clash") {
        col <- clashcol()
      }
      else {
        print("Wrong number of colours specified; using random set of colours instead.")
        col <- randomcol()
      }
    }
    fitted <- x$mean.fit
    index(fitted)=index(fitted)*7
    actual <- zoo(x$aux$y, order.by = x$aux$y.index)
    residuals <- x$std.residuals
    actual.name <- x$aux$y.name
    index(actual)=index(actual)*7
    index(residuals)=index(residuals)*7
    def.par <- par(no.readonly = TRUE)
    par(mar = c(2, 2, 0.5, 0.5))
    if ((x$gets.type == "isat" || coef.path == TRUE) && length(x$ISnames) != 
        0) {
      par(mfrow = c(3, 1))
      is.x <- cbind(x$aux$mX[, x$aux$mXnames %in% x$ISnames])
      is.coef.ests <- coef.isat(x)[x$ISnames]
      coef.path.0 <- zoo(is.x %*% is.coef.ests, order.by = x$aux$y.index)
      index(coef.path.0)=index(coef.path.0)*7
    }
    else {
      par(mfrow = c(2, 1))
    }
    # time.index=(index(fitted)*7)
    # xticks=seq(time.index[1],time.index[length(time.index)],10)
    #  org.ticks<- seq(index(fitted)[1],index(fitted)[length(time.index)],10/7)
    
    if (is.regular(actual)) {
      plot(actual, main = "", ylim = range(min(actual, 
                                               fitted, na.rm = TRUE), max(actual, fitted, na.rm = TRUE)+0.5), 
           type = "l", ylab = "", xlab = "", col = col[2])
      
      #axis(1, at=org.ticks,labels = xticks)
    }
    else {
      plot(as.Date(index(actual)), coredata(actual), main = "", 
           ylim = range(min(actual, fitted, na.rm = TRUE), 
                        max(actual, fitted, na.rm = TRUE)+0.5), type = "l", 
           ylab = "", xlab = "", col = col[2])
      #axis(1,at=org.ticks,labels = xticks)
    }
    if (is.regular(fitted)) {
      lines(fitted, col = col[1],lty=2,lwd=2)
      # axis(1, at=org.ticks,labels = xticks)
    }
    else {
      lines(as.Date(index(fitted)), coredata(fitted), col = col[1],lty=2,lwd=2)
      #  axis(1, at=org.ticks,labels = xticks)
    }
    legend("topleft", lty = lty, lwd = lwd, ncol = 2, col = col[c(2, 
                                                                  1)], legend = c(actual.name, "fitted"), bty = "n")
    outliers1=rep(NA,length(coredata(actual)))
    #index1=c(26,40,95,258)-6
    index1=c(15,22)-6
    outliers1[index1]=actual[index1]
    
    outliers2=rep(NA,length(coredata(actual)))
    index2=c(56,146)-7
    outliers2[index2]=actual[index2]
    
    outliers3=rep(NA,length(coredata(actual)))
    #index3=c(115,151,152,153,157,159,185,213)-7
    index3=c(185)-7
    outliers3[index3]=actual[index3]
    outliers4=rep(NA,length(coredata(actual)))
    index4=c(65)-7
    outliers4[index4]=actual[index4]
   lines(as.Date(index(fitted)),outliers1,col="red",type='o',lwd=2)
   # lines(as.Date(index(actual)),outliers2,col="red",type='o',lwd=2)
    #lines(as.Date(index(actual)),outliers3,col="red",type='o',lwd=2)
    # lines(as.Date(index(actual)),outliers4,col="red",type='o',lwd=2)
    
    text(x=as.Date(index(actual)), y=outliers1, pos=2, labels=c('AO', 'AO','AO','AO','AO','AO','AO'),col="red",font=2)
    #text(x=as.Date(index(actual)), y=outliers2, pos=3, labels=c('LS'),col="red",font=2)
   #text(x=as.Date(index(actual)), y=outliers3, pos=4, labels=c('TC','TC'),col="red",font=2)
   #text(x=as.Date(index(actual)), y=outliers4, pos=4, labels=c('SC','SC','SC'),col="red",font=2)
    #text(locator(), labels = c("AO","AO"))
    #with(outliers, text(sr~dpi, labels = c("AO","AO"), pos = 4))
    
    if (is.regular(residuals)) {
      plot(residuals, type = "h", col = col[1])
      # axis(1, at=org.ticks,labels = xticks)
    }
    else {
      plot(as.Date(index(residuals)), coredata(residuals), 
           type = "h", col = col[1])
      #  axis(1, at=org.ticks,labels = xticks)
    }
    abline(0, 0)
    legend("topleft", lty = 1, col = col[1], legend = "standardised residuals", 
           bty = "n")
    if ((x$gets.type == "isat" | coef.path == TRUE) & length(x$ISnames) != 
        0) {
      if (!is.null(as.list(x$call)$tis) && as.list(x$call)$tis == 
          TRUE) {
        message("\n", appendLF = FALSE)
        message("NB: Because TIS selected, coefficient standard errors invalid hence not plotted", 
                appendLF = TRUE)
        ylim.values <- range(coef.path.0)*7
        if (is.regular(coef.path.0)) {
          ylim.values <- range(coef.path.0)*7
          plot(coef.path.0, type = "l", col = col[1], 
               ylim = ylim.values)
        }
        else {
          plot(as.Date(index(coef.path.0)), coredata(coef.path.0), 
               type = "l", col = col[1], ylim = ylim.values)
        }
      }
      else {
        coef.path.v <- isatvar(x)
        if (is.regular(coef.path.0)) {
          ylim.values <- range(min(coef.path.0 - qt(0.975, 
                                                    NROW(coef.path.0)) * coef.path.v$const.se)+10, 
                               max(coef.path.0 + qt(0.975, NROW(coef.path.0)-10) * 
                                     coef.path.v$const.se))
          plot(coef.path.0, type = "l", col = col[1], 
               ylim = ylim.values)
          lines(coef.path.0 + qt(0.975, NROW(coef.path.0)) * 
                  coef.path.v$const.se, type = "l", col = col[1], 
                lty = 3)
          lines(coef.path.0 - qt(0.975, NROW(coef.path.0)) * 
                  coef.path.v$const.se, type = "l", col = col[1], 
                lty = 3)
        }
        else {
          plot(as.Date(index(coef.path.0)), coredata(coef.path.0), 
               type = "l", col = col[1], ylim = ylim.values)
          lines(as.Date(index(coef.path.0)), coredata(coef.path.0) + 
                  qt(0.975, NROW(coef.path.0)) * coef.path.v$const.se, 
                type = "l", col = col[1], lty = 3)
          lines(as.Date(index(coef.path.0)), coredata(coef.path.0) - 
                  qt(0.975, NROW(coef.path.0)) * coef.path.v$const.se, 
                type = "l", col = col[1], lty = 3)
        }
      }
      abline(0, 0, lty = 3)
      legend("topleft", lty = 1, col = col[1], legend = c(paste(actual.name, 
                                                                "Coefficient Path", sep = ": ")), bty = "n")
    }
    par(def.par)
  }
}
