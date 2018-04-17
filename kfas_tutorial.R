###EXAMple from KFAS vignette
### April 2018
### Author ZUzana Leova

#install.packages("KFAS")
require(KFAS)
library(ggplot2)
#local level with drift model
#mu_{t+1}=mu_{t}+nu+zeta_{t} where zeta_{t}~N(0,sigma^2_{zeta})

data("alcohol")
View("alcohol")
str(alcohol)

#?window

#"death at age 40-49"
deaths <- window(alcohol[,2], end = 2007) #window is a function from stats, extract the times 
                                          #between start and end 

# population by age 40-49 
population <- window(alcohol[, 6], end = 2007)
Zt <- matrix(c(1, 0), 1, 2)
Ht <- matrix(NA)
Tt <- matrix(c(1, 0, 1, 1), 2, 2)
Rt <- matrix(c(1, 0), 2, 1)
Qt <- matrix(NA)
a1 <- matrix(c(1, 0), 2, 1)
P1 <- matrix(0, 2, 2)
P1inf <-diag(2)
model_gaussian<-SSModel(deaths/population ~ -1 +SSMcustom(Z = Zt, T = Tt, 
                                       R = Rt, Q = Qt, a1 = a1, P1 = P1,
                                        P1inf = P1inf),H = Ht)

#estimation of unknown variance sigma_eps and sigma_eta
#BFGS is quasi-Newton methos
fit_gaussian <-fitSSM(model_gaussian,inits=c(0,0),method="BFGS")

#smoothing and filtering recursions
out_gaussian <- KFS(fit_gaussian$model) #slope estimate 0.8409 with st 0.34

ggplot(deaths/population, aes(x = Home.Value)) +
  geom_histogram()
plot(deaths/population,plot.type ="multiple")
plot(out_gaussian$a[1:39],col="red")
lines(out_gaussian$alphahat[1:39],col="blue")

######################################
## Exponential family distributions ##
######################################

#poisson distribution
model_poisson <- SSModel(deaths ~ -1 +
                  SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, P1inf = P1inf), 
                distribution="poisson", u = population)