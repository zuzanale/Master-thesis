plot.my.new=function (x,col = c("gray", "black"), lty = c(1, 2), 
                   lwd = c(1, 2), coef.path = TRUE, x.eps,x.eta, x.omega) 
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
    if ((x.eps$gets.type == "isat" || coef.path == TRUE) && length(x.eps$ISnames) != 
        0) {
     # par(mfrow = c(3, 1))
      is.x <- cbind(x.eps$aux$mX[, x.eps$aux$mXnames %in% x.eps$ISnames])
      is.coef.ests <- coef.isat(x.eps)[x.eps$ISnames]
      coef.path.0.eps <- zoo(is.x %*% is.coef.ests, order.by = x.eps$aux$y.index)
      index(coef.path.0.eps)=index(coef.path.0.eps)*7
    }
    if ((x.eta$gets.type == "isat" || coef.path == TRUE) && length(x.eta$ISnames) != 
        0) {
      # par(mfrow = c(3, 1))
      is.x <- cbind(x.eta$aux$mX[, x.eta$aux$mXnames %in% x.eta$ISnames])
      is.coef.ests <- coef.isat(x.eta)[x.eta$ISnames]
      coef.path.0.eta <- zoo(is.x %*% is.coef.ests, order.by = x.eta$aux$y.index)
      index(coef.path.0.eta)=index(coef.path.0.eta)*7
    }
    if ((x.omega$gets.type == "isat" || coef.path == TRUE) && length(x.omega$ISnames) != 
        0) {
      # par(mfrow = c(3, 1))
      is.x <- cbind(x.omega$aux$mX[, x.omega$aux$mXnames %in% x.omega$ISnames])
      is.coef.ests <- coef.isat(x.omega)[x.omega$ISnames]
      coef.path.0.omega <- zoo(is.x %*% is.coef.ests, order.by = x.omega$aux$y.index)
      index(coef.path.0.omega)=index(coef.path.0.omega)*7
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
      
    lines(outliers,col='red',type='o')
    
    
    ################# do we want to keep residuals????
   # if (is.regular(residuals)) {
   #   plot(residuals, type = "h", col = col[1])
      # axis(1, at=org.ticks,labels = xticks)
    #}
    #else {
     # plot(as.Date(index(residuals)), coredata(residuals), 
      #     type = "h", col = col[1])
      #  axis(1, at=org.ticks,labels = xticks)
    #}
    
  #  abline(0, 0)
   # legend("topleft", lty = 1, col = col[1], legend = "standardised residuals", 
    #       bty = "n")
   #       abline(0, 0, lty = 3)
   #
    }
    par(def.par)
}

outliers=rep(0,length(actual))
for(t in index(coef.path.0.eps)){
  print(t)
  if(coef.path.0.eps[t]>0){#tu potrebujeme dat looping cez index coeficientu nie actualu, 8 indexov posun for some reason
    sign.t=index(coef.path.0.eps[t])
    outliers[sign.t]=actual[sign.t]
  }
}

outliers=rep(NA,length(actual))
for(t in 1:length(coef.path.0.eps)){
  #print(t)
  if(coef.path.0.eps[t]!=0){#tu potrebujeme dat looping cez index coeficientu nie actualu, 8 indexov posun for some reason
    sign.t=t+9
    print(t)
    outliers[sign.t]=actual[sign.t]
  }
}
#index(outliers)=index(actual)
plot(actual)
lines(as.Date(index(actual)),outliers,col="red",type='o')

#plot.my.new(indicator1.res,x.eps=out.bsc.eps1,x.eta=out.bsc.eta1,x.omega=out.bsc.zeta1,out.bsc.omega1)
#plot.my.new(indicator2.res,x.eps=out.bsc.eps2,x.eta=out.bsc.eta2,x.omega=out.bsc.zeta2,out.bsc.omega2)
actual <- zoo(indicator1.res$aux$y, order.by = indicator1.res$aux$y.index)
index(actual)=index(actual)*7
x.eps=out.bsc.eps1
is.x <- cbind(x.eps$aux$mX[, x.eps$aux$mXnames %in% x.eps$ISnames])
is.coef.ests <- coef.isat(x.eps)[x.eps$ISnames]
coef.path.0.eps <- zoo(is.x %*% is.coef.ests, order.by = x.eps$aux$y.index)

plot(actual)
lines(as.Date(index(actual)),outliers,col='red',type='o')
index(coef.path.0.eps)=index(coef.path.0.eps)*7