library(MASS)
summary(geyser)
data1d = geyser[,1]
h1=4
data1d_kde = density(data1d, bw=h1, from=40, to=110)
## bw = h1: we set the smoothing banwdith = h1.
## from = 40, to = 110: the range we are evaluating the density is from 40 to 110.

##### make a plot
par(mar=c(4,4,2,1))
plot(data1d_kde, lwd=4, col="blue",ylim=c(0, 0.045), cex.axis=1.5,
     main="Kernel Density Estimator", cex.lab=2, cex.main=1.5, ylab="")
mtext("Density", side=2, line=2.2, cex=2)

#ks library to choose the bandwidth
#hns() for a multivariate case to choose the smoothing bandwidth according to the Silvermans rule
#hlscv() least square cross-validation approach
#hpi() applies the plug in method for bandwidth selection

library(ks)
h2 = hns(data1d)
h3 = hlscv(data1d)
h4 = hpi(data1d)
data1d_kde2 = density(data1d, bw=h2, from=40, to=110)
data1d_kde3 = density(data1d, bw=h3, from=40, to=110)
data1d_kde4 = density(data1d, bw=h4, from=40, to=110)
print(h2)
print(h3)
print(h4)

##### make a plot
par(mar=c(4,4,2,1))
plot(data1d_kde2, lwd=4, col="dodgerblue",ylim=c(0, 0.045), cex.axis=1.5,
     main="Kernel Density Estimator", cex.lab=2, cex.main=1.5, ylab="",
     xlab="")
lines(data1d_kde3, lwd=4, col="limegreen")
lines(data1d_kde4, lwd=4, col="brown")
mtext("Density", side=2, line=2.2, cex=2)
legend("topleft",c("hns()","hlscv()","hpi()"), col=c("dodgerblue","limegreen","brown"), lwd=5)

#confidence intervals construction with plug in approach or bootstrap approach
alpha0=0.05
n1 = length(data1d)
t0 = qnorm(1-alpha0)*sqrt(1/(2*sqrt(pi))/(n1*h1)*data1d_kde$y)

## data1d_kde$y: the density value evaluated at each point of data1d_kde$x.

##### make a plot
par(mar=c(4,4,2,1))
plot(data1d_kde, lwd=3, col="tan4", ylim=c(0, 0.045), cex.axis=1.5,
     main="Confidence Interval (PI)", xlab="", ylab="", cex.main=1.5)
mtext("Density", side=2, line=2.2, cex=2)
polygon(c(data1d_kde$x, rev(data1d_kde$x)),
        c(data1d_kde$y+t0, rev(data1d_kde$y-t0)),
        border="tan", col="tan")
abline(h=0, lwd=2)
lines(data1d_kde$x, data1d_kde$y, lwd=3, col="tan4")

#########################
# bootstrap approach ####
########################

n_BT = 1000

## number of bootstrap samples
kde_seq_BT_m = matrix(NA, nrow=n_BT, ncol= length(data1d_kde$y))
for(j in 1:n_BT){
  data1d_BT = data1d[sample(n1, n1, replace=T)]
  data1d_kde_BT = density(data1d_BT, bw=h1, from=40, to=110)
  kde_seq_BT_m[j,] = abs(data1d_kde_BT$y-data1d_kde$y)
}

## each row of 'kde_seq_BT_m' contains one bootstrap difference
t_pt = rep(0, length(data1d_kde$y))
for(l in 1:length(data1d_kde$y)){
  t_pt[l] = quantile(kde_seq_BT_m[,l], 1-alpha0)
}
## t_pt: the 1-\alpha quantile of deviation at each point

##### make a plot
par(mar=c(4,4,2,1))
plot(data1d_kde, lwd=3, col="blue", ylim=c(0, 0.045), cex.axis=1.5,
     main="Confidence Interval (Bootstrap)", xlab="", ylab="", cex.main=1.5)
mtext("Density", side=2, line=2.2, cex=2)
polygon(c(data1d_kde$x, rev(data1d_kde$x)),
        c(data1d_kde$y+t_pt, rev(data1d_kde$y-t_pt)),
        border="skyblue", col="skyblue")
abline(h=0, lwd=2)
lines(data1d_kde$x, data1d_kde$y, lwd=3, col="blue")

######################################################
### Confidence Band ##################################
######################################################

#Traditional boostrap approach

kde_seq_sup = rep(0, n_BT)
for(j in 1:n_BT){
  data1d_BT = data1d[sample(n1, n1, replace=T)]
  data1d_kde_BT = density(data1d_BT, bw=h1, from=40, to=110)
  ## bootstrap KDE
  kde_seq_sup[j] = max(abs(data1d_kde_BT$y-data1d_kde$y))
}
t_sup = quantile(kde_seq_sup, 1-alpha0)

##### make a plot
par(mar=c(4,4,2,1))
plot(data1d_kde$x, data1d_kde$y, lwd=3, col="purple", ylim=c(0, 0.045), cex.axis=1.5,
     main="Confidence Band", xlab="", ylab="", type="l", cex.main=1.5)
mtext("Density", side=2, line=2.2, cex=2)
polygon(c(data1d_kde$x, rev(data1d_kde$x)),c(data1d_kde$y+t_sup, rev(data1d_kde$y-t_sup)),
        border="plum1", col="plum1")
abline(h=0, lwd=2)
lines(data1d_kde$x, data1d_kde$y, lwd=3, col="purple")

#bootstrapping with debiased KDE approach
library(ks)

## this library provides density derivative estimation
X0_kdde = kdde(data1d, h=h1, eval.points = data1d_kde$x, deriv.order = 2)
## eval.points = data1d_kde$x: we evalute the result at data1d_kde$x
## deriv.order = 2: we are computing the second derivative

de_kde_seq = data1d_kde$y-0.5*h1^2*X0_kdde$estimate
## this is the debiased KDE, evaluated at points of data1d_kde$x

de_kde_seq_sup = rep(0, n_BT)
for(j in 1:n_BT){
  data1d_BT = data1d[sample(n1, n1, replace=T)]
  data1d_kde_BT = density(data1d_BT, bw=h1, from=40, to=110)
  data1d_kdde_BT = kdde(data1d_BT, h=h1, eval.points = data1d_kde$x, deriv.order = 2)
  de_kde_seq_BT = data1d_kde_BT$y-0.5*h1^2*data1d_kdde_BT$estimate
  ## the bootstrap debiased KDE
  de_kde_seq_sup[j] = max(abs(de_kde_seq-de_kde_seq_BT))
}
t_de_sup = quantile(de_kde_seq_sup, 1-alpha0)

##### make a plot
par(mar=c(4,4,2,1))
plot(data1d_kde$x, de_kde_seq, lwd=3, col="seagreen", ylim=c(0, 0.045), cex.axis=1.5,
     main="Confidence Band (Debiased)", xlab="", ylab="", type="l", cex.main=1.5)
mtext("Density", side=2, line=2.2, cex=2)
polygon(c(data1d_kde$x, rev(data1d_kde$x)),
        c(de_kde_seq+t_de_sup, rev(de_kde_seq-t_de_sup)),
        border="palegreen", col="palegreen")
abline(h=0, lwd=2)
lines(data1d_kde$x, de_kde_seq, lwd=3, col="seagreen")