## load zoo
require(zoo)
require(forecast)

####Number of new requests for NCSA  
ts1=score1[[1]][,c(4,5)]
ts1=ts1[1:180,]
ts1$start_at=as.Date(ts1$start_at)
## convert to a zoo object, with order given by the `datefield`
ts1.zoo=with(ts1, zoo(value, order.by = start_at))
ts1.ts=as.ts(ts1.zoo)
plot(ts1.zoo)

# weekly seasonality and deterministic trend adjustment
timeSeries = ts(ts1.zoo,frequency=7)#1067

####Summed scores of surveys for incidents resolved 30 days ago for APAC   
ts2=score2[[1]][,c(4,5)]
ts2$start_at=as.Date(ts2$start_at)
#ts2=ts2[1:180,]
## convert to a zoo object, with order given by the `datefield`
ts2.zoo=with(ts2, zoo(value, order.by = start_at))
ts2.ts=ts(ts2.zoo,frequency = 7)
plot(ts2.zoo,ylab='',xlab='Time')
plot(log(ts2.zoo),ylab='',xlab='Time')
# weekly seasonality and deterministic trend adjustment
timeSeries2 = ts(ts2.zoo,frequency=7)#275 we will take a time frame of 3 years
plot(timeSeries2)


####breakdown 2 type 2  #domain:local support
ts2.21=score2.breakdown21[[1]][,c(4,5)]
#ts2.11=ts2.11[1:180,]
ts2.21$start_at=as.Date(ts2.21$start_at)
## convert to a zoo object, with order given by the `datefield`
ts2.21.zoo=with(ts2.21, zoo(value, order.by = start_at))
ts2.21.ts=as.ts(ts2.21.zoo)
plot(ts2.21.zoo,ylab='',xlab='Time')
plot(log(ts2.21.zoo),ylab='',xlab='Time')
# weekly seasonality and deterministic trend adjustment
timeSeries2.21 = ts(ts2.21.zoo,frequency=7)

####breakdown 1 #department/studio: GNS
ts2.11=score2.breakdown11[[1]][,c(4,5)]
#ts2.11=ts2.11[1:180,]
ts2.11$start_at=as.Date(ts2.11$start_at)
## convert to a zoo object, with order given by the `datefield`
ts2.11.zoo=with(ts2.11, zoo(value, order.by = start_at))
ts2.11.ts=as.ts(ts2.11.zoo)
plot(ts2.11.zoo,ylab='',xlab='Time')
plot(log(ts2.11.zoo),ylab='',xlab='Time')

# weekly seasonality and deterministic trend adjustment
timeSeries2.11 = ts(ts2.11.zoo,frequency=7)

####breakdown 2 type 2  #domain:infrastructure
ts2.22=score2.breakdown22[[1]][,c(4,5)]
#ts2.11=ts2.11[1:180,]
ts2.22$start_at=as.Date(ts2.22$start_at)
## convert to a zoo object, with order given by the `datefield`
ts2.22.zoo=with(ts2.22, zoo(value, order.by = start_at))
ts2.22.ts=as.ts(ts2.22.zoo)
plot(ts2.22.zoo,ylab='',xlab='Time')
plot(log(ts2.22.zoo),ylab='',xlab='Time')

# weekly seasonality and deterministic trend adjustment
timeSeries2.22 = ts(ts2.22.zoo,frequency=7)

####breakdown 2 type 3  #domain:enterprise
ts2.23=score2.breakdown23[[1]][,c(4,5)]
#ts2.11=ts2.11[1:180,]
ts2.23$start_at=as.Date(ts2.23$start_at)
## convert to a zoo object, with order given by the `datefield`
ts2.23.zoo=with(ts2.23, zoo(value, order.by = start_at))
ts2.23.ts=as.ts(ts2.23.zoo)
plot(ts2.23.zoo,ylab='',xlab='Time')

# weekly seasonality and deterministic trend adjustment
timeSeries2.23 = ts(ts2.23.zoo,frequency=7)

indicator2=cbind(ts2.zoo,ts2.11.zoo,ts2.21.zoo,ts2.22.zoo)

#### Summed duration between Planned end and Closure 
ts3=score3[[1]][,c(4,5)]
#ts3=ts3[1:90,]
ts3$start_at=as.Date(ts3$start_at)
## convert to a zoo object, with order given by the `datefield`
ts3.zoo=with(ts3, zoo(value, order.by = start_at))
ts3.ts=as.ts(ts3.zoo)
plot(ts3.zoo,ylab='',xlab='Time')
lts3.zoo=log(ts3.zoo)
lts3.zoo[lts3.zoo==-Inf]=0
plot(lts3.zoo,ylab='',xlab='Time')

# weekly seasonality and deterministic trend adjustment
timeSeries3 = ts(ts3.zoo,frequency=7)

####breakdown1 
ts3.11= score3.breakdown11[[1]][,c(4,5)]
#ts3.11=ts3.11[(nrow(ts3)-90):nrow(ts3),]
ts3.11$start_at=as.Date(ts3.11$start_at)
## convert to a zoo object, with order given by the `datefield`
ts3.11.zoo=with(ts3.11, zoo(value, order.by = start_at))
ts3.11.ts=as.ts(ts3.11.zoo)
plot(ts3.11.zoo,ylab='',xlab='Time')
lts3.11.zoo=log(ts3.11.zoo)
lts3.11.zoo[lts3.11.zoo==-Inf]=0
plot(lts3.11.zoo,ylab='',xlab='Time')

# weekly seasonality and deterministic trend adjustment
timeSeries3.11 = ts(ts3.11.zoo,frequency=7)

####breakdown 2  
ts3.12= score3.breakdown12[[1]][,c(4,5)]
#ts3.11=ts3.11[(nrow(ts3)-90):nrow(ts3),]
ts3.12$start_at=as.Date(ts3.12$start_at)
## convert to a zoo object, with order given by the `datefield`
ts3.12.zoo=with(ts3.12, zoo(value, order.by = start_at))
ts3.12.ts=as.ts(ts3.12.zoo)
plot(ts3.12.zoo,ylab='',xlab='Time')
lts3.12.zoo=log(ts3.12.zoo)
lts3.12.zoo[lts3.12.zoo==-Inf]=0
plot(lts3.12.zoo,ylab='',xlab='Time')

# weekly seasonality and deterministic trend adjustment
timeSeries3.12 = ts(ts3.12.zoo,frequency=7)

####breakdown 3  
ts3.13= score3.breakdown13[[1]][,c(4,5)]
#ts3.11=ts3.11[(nrow(ts3)-90):nrow(ts3),]
ts3.13$start_at=as.Date(ts3.13$start_at)
## convert to a zoo object, with order given by the `datefield`
ts3.13.zoo=with(ts3.13, zoo(value, order.by = start_at))
ts3.13.ts=as.ts(ts3.13.zoo)
plot(ts3.13.zoo,ylab='',xlab='Time')
lts3.13.zoo=log(ts3.13.zoo)
lts3.13.zoo[lts3.13.zoo==-Inf]=0
plot(lts3.13.zoo,ylab='',xlab='Time')

# weekly seasonality and deterministic trend adjustment
timeSeries3.13 = ts(ts3.13.zoo,frequency=7)

indicator3=cbind(ts3.zoo,ts3.11.zoo,ts3.12.zoo,ts3.13.zoo)

####Summed scores of new surveys for GNS NCSA  Time Series (d = days)  
ts4=score4[[1]][,c(4,5)]
ts4$start_at=as.Date(ts4$start_at)
#ts2=ts2[1:180,]
## convert to a zoo object, with order given by the `datefield`
ts4.zoo=with(ts4, zoo(value, order.by = start_at))
ts4.ts=ts(ts4.zoo,frequency = 7)
plot(ts4.zoo,ylab='',xlab='Time') #this one is better for presentation
lts4.zoo=log(ts4.zoo)
lts4.zoo[lts4.zoo==-Inf]=0
plot(lts4.zoo,ylab='',xlab='Time')

# weekly seasonality and deterministic trend adjustment
timeSeries4 = ts(ts4.zoo,frequency=7)#275 we will take a time frame of 3 years
plot(timeSeries4) #244 obs

####Summed scores of new surveys for GNS NCSA  Time Series (d = days)  
ts5=score5[[1]][,c(4,5)]
ts5$start_at=as.Date(ts5$start_at)
#ts2=ts2[1:180,]
## convert to a zoo object, with order given by the `datefield`
ts5.zoo=with(ts5, zoo(value, order.by = start_at))
ts5.ts=ts(ts5.zoo,ylab='',xlab='Time')
plot(ts5.zoo) #this one is better for presentation
lts5.zoo=log(ts5.zoo)
lts5.zoo[lts5.zoo==-Inf]=0
plot(lts5.zoo,ylab='',xlab='Time')

# weekly seasonality and deterministic trend adjustment
timeSeries5 = ts(ts5.zoo,frequency=7)#275 we will take a time frame of 3 years
plot(timeSeries5) #244 obs

stargazer(as.data.frame(ts5.zoo),digits=2)
indicator3
ts4.zoo
ts5.zoo
data=do.call(merge,list(indicator2,indicator3,ts4.zoo,ts5.zoo))

hist(ts4.zoo, # histogram
     col="lightblue", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "",
     main = "")
ts2.21
hist(ts2.22.zoo, # histogram
     col="lightblue", # column color
     border="black",
     prob = TRUE, # show densities instead of frequencies
     xlab = "",
     main = "",
     ylab = ""
    #,ylim=c(0,0.08)
)
     lines(density(na.omit(ts2.22.zoo)), # density plot
      lwd = 2, # thickness of line
      col = "chocolate3",xlim=c(0.01,1500),
      ylab = "")
plot(density(na.omit(ts5.zoo)),main='')
par(mfrow = c(2,3))


####testing out algorithm for outlier detection
## indicator 1
lts2.zoo=log(ts2.zoo)
indicator1 = ts(lts2.zoo,frequency =7) 
indicator1.res=isat_my(indicator1, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
        vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
        optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
indicator1 = ts(ts2.zoo,frequency =7) 
indicator1.res.level=isat_my(indicator1, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                       vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                       optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))

plot.my2(indicator1.res)
 out.bsc.eps1=isat(indicator1.res$pure$std.eps,sis=FALSE,iis=TRUE,ar=2)
  plot.my2(out.bsc.eps1)
out.bsc.eta1=isat(diff(indicator1.res$pure$std.eta),sis=F,iis=T,ar=1)#for the level
plot.my(out.bsc.eta1)
 out.bsc.zeta1=isat(diff(indicator1.res$pure$std.zeta),sis=F,iis=T,ar=4)#for the slope
  plot.my(out.bsc.zeta1)
  out.bsc.omega1=isat(indicator1.res$pure$std.omega,sis=F,iis=T)#for seasonality
  plot.my(out.bsc.omega1)
  
  #SIS approach
  out.bsc.eps1.SIS=isat(indicator1.res$pure$std.eps,sis=TRUE,iis=FALSE)
  plot(out.bsc.eps1.SIS)
  out.bsc.eta1.SIS=isat(indicator1.res$pure$std.eta,sis=TRUE,iis=FALSE)#for the level
  plot(out.bsc.eta1.SIS)
  out.bsc.zeta1.SIS=isat(indicator1.res$pure$std.zeta,sis=TRUE,iis=FALSE)#for the slope
  plot(out.bsc.zeta1.SIS)
  out.bsc.omega1.SIS=isat(indicator1.res$pure$std.omega,sis=TRUE,iis=FALSE)#for seasonality
  plot(out.bsc.omega1.SIS)
  
  #alternative methods
 timeSeries1=indicator1 
  timeSeries1[is.na(timeSeries1)]=0
  #TS=decompose(timeSeries)
  TS1=decompose(timeSeries1,type="additive")
  y=TS1$x-TS1$seasonal-TS1$trend
  
  #testing for automatic order of ARIMA
  auto.arima(timeSeries1)
  auto.arima(y)
  
  ###IS with gets package
  indicator1.IS.level=isat(timeSeries1, ar=1,t.pval=1/T, turbo=TRUE,iis=T,sis=T,tis=T,
                     vcov.type = "ordinary",info.method = 'aic')
  indicator1.IS=isat(timeSeries1, ar=1,t.pval=1/T, turbo=TRUE,iis=T,sis=T,tis=T,
                           vcov.type = "ordinary",info.method = 'aic')
  
  indicator1.tso.level=tso(y = timeSeries1, types = c("AO", "LS", "TC","SLS"),
                     tsmethod = "stsm", args.tsmodel = list(model = "BSM")
  )
  
  indicator1.tso=tso(y = timeSeries1, types = c("AO", "LS", "TC","SLS"),
                           tsmethod = "stsm", args.tsmodel = list(model = "BSM")
  )
  
  #Sumanas approach 
  indicator1.WOM.level=outlier(indicator1)
  indicator1.WOM=outlier(indicator1)
  #### breakdown 1.11
  par(mfrow=c(1,1))
  lts2.11.zoo=log(ts2.11.zoo)
  indicator1.11 = ts(lts2.11.zoo,frequency =7) 
  indicator1.11.res=isat_my(indicator1.11, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                         vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                         optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
  indicator1.11 = ts(ts2.11.zoo,frequency =7) 
  indicator1.11.res.level=isat_my(indicator1.11, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                            vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                            optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
  
  plot(indicator1.11.res)
  out.bsc.eps1.11=isat(indicator1.11.res$pure$std.eps,sis=F,iis=T)
  plot(out.bsc.eps1.11)
  out.bsc.eta1.11=isat(diff(indicator1.11.res$pure$std.eta),sis=FALSE,iis=TRUE)#for the level
  plot(out.bsc.eta1.11)
  out.bsc.zeta1.11=isat(diff(indicator1.11.res$pure$std.zeta),sis=FALSE,iis=TRUE,ar=4)#for the slope
  plot(out.bsc.zeta1.11)
  out.bsc.omega1.11=isat(indicator1.11.res$pure$std.omega,sis=FALSE,iis=TRUE,ar=2)#for seasonality
  plot(out.bsc.omega1.11)
  
  #SIS approach
  out.bsc.eps1.11.SIS=isat(indicator1.11.res$pure$std.eps,sis=TRUE,iis=FALSE)
  plot(out.bsc.eps1.11.SIS)
  out.bsc.eta1.11.SIS=isat(diff(indicator1.11.res$pure$std.eta),sis=TRUE,iis=FALSE)#for the level
  plot(out.bsc.eta1.11.SIS)
  out.bsc.zeta1.11.SIS=isat(diff(indicator1.11.res$pure$std.zeta),sis=TRUE,iis=FALSE,ar=4)#for the slope
  plot(out.bsc.zeta1.11.SIS)
  out.bsc.omega1.11.SIS=isat(indicator1.11.res$pure$std.omega,sis=TRUE,iis=FALSE,ar=2)#for seasonality
  plot(out.bsc.omega1.11.SIS)
  
  #alternative methods
  timeSeries1.11=indicator1.11
  timeSeries1.11[is.na(timeSeries1.11)]=0
  #TS=decompose(timeSeries)
  TS1.11=decompose(timeSeries1.11,type="additive")
  y1.11=TS1.11$x-TS1.11$seasonal-TS1.11$trend
  
  #testing for automatic order of ARIMA
  auto.arima(timeSeries1.11)
  auto.arima(y1.11)
  
  ###IS with gets package
  indicator1.11.IS.level=isat(timeSeries1.11, ar=1,t.pval=1/T, turbo=TRUE,iis=T,sis=T,tis=T,
                     vcov.type = "ordinary",info.method = 'aic')
  
 
  indicator1.11.IS=isat(timeSeries1.11, ar=1,t.pval=1/T, turbo=TRUE,iis=T,sis=T,tis=T,
                              vcov.type = "ordinary",info.method = 'aic')
  
  indicator1.11.tso.level=tso(y = timeSeries1.11, types = c("AO", "LS", "TC","SLS"),
                     tsmethod = "stsm", args.tsmodel = list(model = "BSM")
  )
  
  indicator1.11.tso=tso(y = timeSeries1.11, types = c("AO", "LS", "TC","SLS"),
                              tsmethod = "stsm", args.tsmodel = list(model = "BSM")
  )
  
  indicator1.11.WOM.level=outlier(indicator1.11)
  indicator1.11.WOM=outlier(indicator1.11)
  
  #### breakdown 1.21
  par(mfrow=c(1,1))
  lts2.21.zoo=log(ts2.21.zoo)
  indicator1.21 = ts(lts2.21.zoo,frequency =7) 
  indicator1.21.res=isat_my(indicator1.21, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                            vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                            optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE),sis=FALSE,iis=TRUE)
  indicator1.21 = ts(ts2.21.zoo,frequency =7) 
  indicator1.21.res.level=isat_my(indicator1.21, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                            vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                            optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
  plot(indicator1.21.res)
  out.bsc.eps1.21=isat(indicator1.21.res$pure$std.eps,sis=FALSE,iis=TRUE)
  plot(out.bsc.eps1.21)
  out.bsc.eta1.21=isat(diff(indicator1.21.res$pure$std.eta),sis=FALSE,iis=TRUE)#for the level
  plot(out.bsc.eta1.21)
  out.bsc.zeta1.21=isat(diff(indicator1.21.res$pure$std.zeta),sis=FALSE,iis=TRUE,ar=2)#for the slope
  plot(out.bsc.zeta1.21)
  out.bsc.omega1.21=isat(indicator1.21.res$pure$std.omega,sis=FALSE,iis=TRUE,ar=4)#for seasonality
  plot(out.bsc.omega1.21)
  
  #SIS
  out.bsc.eps1.21.SIS=isat(indicator1.21.res$pure$std.eps,sis=TRUE,iis=FALSE)
  plot(out.bsc.eps1.21.SIS)
  out.bsc.eta1.21.SIS=isat(diff(indicator1.21.res$pure$std.eta),sis=TRUE,iis=FALSE)#for the level
  plot(out.bsc.eta1.21.SIS)
  out.bsc.zeta1.21.SIS=isat(diff(indicator1.21.res$pure$std.zeta),sis=TRUE,iis=FALSE,ar=2)#for the slope
  plot(out.bsc.zeta1.21.SIS)
  out.bsc.omega1.21.SIS=isat(indicator1.21.res$pure$std.omega,sis=FALSE,iis=FALSE,ar=4)#for seasonality
  plot(out.bsc.omega1.21.SIS)
  
  #alternative methods
  timeSeries1.21=indicator1.21
  timeSeries1.21[is.na(timeSeries1.21)]=0
  #TS=decompose(timeSeries)
  TS1.21=decompose(timeSeries1.21,type="additive")
  y1.21=TS1.21$x-TS1.21$seasonal-TS1.21$trend
  
  #testing for automatic order of ARIMA
  auto.arima(timeSeries1.21)
  auto.arima(y1.21)
  
  ###IS with gets package
  indicator1.21.IS.level=isat(timeSeries1.21, ar=1,t.pval=1/T, turbo=TRUE,iis=T,sis=T,tis=T,
                        vcov.type = "ordinary",info.method = 'aic')
  
  indicator1.21.tso.level=tso(y = timeSeries1.21, types = c("AO", "LS", "TC","SLS"),
                        tsmethod = "stsm", args.tsmodel = list(model = "BSM")
  )
  
  ###IS with gets package
  indicator1.21.IS=isat(timeSeries1.21, ar=1,t.pval=1/T, turbo=TRUE,iis=T,sis=T,tis=T,
                              vcov.type = "ordinary",info.method = 'aic')
  
  indicator1.21.tso=tso(y = timeSeries1.21, types = c("AO", "LS", "TC","SLS"),
                              tsmethod = "stsm", args.tsmodel = list(model = "BSM")
  )
  
  indicator1.21.WOM=outlier(indicator1.21)
  indicator1.21.WOM.level=outlier(indicator1.21)
  
  #### breakdown 1.22
  par(mfrow=c(1,1))
  lts2.22.zoo=log(ts2.22.zoo)
  indicator1.22 = ts(lts2.22.zoo,frequency =7) 
  indicator1.22.res=isat_my(indicator1.22, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                            vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                            optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
  plot.my(indicator1.22.res)
  
  #level
  indicator1.22 = ts(ts2.22.zoo,frequency =7) 
  indicator1.22.res.level=isat_my(indicator1.22, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                            vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                            optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
  
  out.bsc.eps1.22=isat(indicator1.22.res$pure$std.eps,sis=FALSE,iis=TRUE)
  plot(out.bsc.eps1.22)
  out.bsc.eta1.22=isat(diff(indicator1.22.res$pure$std.eta),ar=1,sis=FALSE,iis=TRUE)#for the level
  plot(out.bsc.eta1.22)
  out.bsc.zeta1.22=isat(diff(indicator1.22.res$pure$std.zeta),sis=FALSE,iis=TRUE,ar=1)#for the slope
  plot(out.bsc.zeta1.22)
  out.bsc.omega1.22=isat(indicator1.22.res$pure$std.omega,sis=FALSE,iis=TRUE,ar=2)#for seasonality
  plot.my(out.bsc.omega1.22)
  
  #SIS
  out.bsc.eps1.22.SIS=isat(indicator1.22.res$pure$std.eps,sis=TRUE,iis=FALSE)
  plot(out.bsc.eps1.22.SIS)
  out.bsc.eta1.22.SIS=isat(diff(indicator1.22.res$pure$std.eta),sis=TRUE,iis=FALSE)#for the level
  plot(out.bsc.eta1.22.SIS)
  out.bsc.zeta1.22.SIS=isat(diff(indicator1.22.res$pure$std.zeta),sis=TRUE,iis=FALSE,ar=1)#for the slope
  plot(out.bsc.zeta1.22.SIS)
  out.bsc.omega1.22.SIS=isat(indicator1.22.res$pure$std.omega,sis=TRUE,iis=FALSE)#for seasonality
  plot(out.bsc.omega1.22.SIS)
  
  #alternative methods
  timeSeries1.22=indicator1.22
  timeSeries1.22[is.na(timeSeries1.22)]=0
  #TS=decompose(timeSeries)
  TS1.22=decompose(timeSeries1.22,type="additive")
  y1.21=TS1.22$x-TS1.22$seasonal-TS1.22$trend
  
  #testing for automatic order of ARIMA
  auto.arima(timeSeries1.22)
  auto.arima(y1.22)
  
  ###IS with gets package
  indicator1.22.IS.level=isat(timeSeries1.22, ar=1,t.pval=1/T, turbo=TRUE,iis=TRUE,sis=TRUE,tis=TRUE,
                        vcov.type = "ordinary",info.method = 'aic')
  
  indicator1.22.tso.level=tso(y = timeSeries1.22, types = c("AO", "LS", "TC","SLS"),
                        tsmethod = "stsm", args.tsmodel = list(model = "BSM")
  )
  ###IS with gets package
  indicator1.22.IS=isat(timeSeries1.22, ar=1,t.pval=1/T, turbo=TRUE,iis=TRUE,sis=TRUE,tis=TRUE,
                              vcov.type = "ordinary",info.method = 'aic')
  
  indicator1.22.tso=tso(y = timeSeries1.22, types = c("AO", "LS", "TC","SLS"),
                              tsmethod = "stsm", args.tsmodel = list(model = "BSM")
  )
  
  indicator1.22.WOM.level=outlier(indicator1.22)
  indicator1.22.WOM=outlier(indicator1.22)
  
  ###################
  ## indicator 2
  ###################
  lts3.zoo=log(ts3.zoo)
  indicator2 = ts(lts3.zoo,frequency =7) 
  start_time <- Sys.time()
  indicator2.res=isat_my(indicator2, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                       vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                         optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
  end_time <- Sys.time()
  print(end_time - start_time)
  
  plot.my2(indicator2.res)
  
  #level
  indicator2 = ts(ts3.zoo,frequency =7) 
  indicator2.res.level=isat_my(indicator2, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                               vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                               optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
  plot.my(indicator2.res.level)
  #version without autocorrelated residuals
  out.bsc.eps2=isat((indicator2.res$pure$std.eps),sis=FALSE,iis=TRUE,ar=2)
  plot.my(out.bsc.eps2)
  out.bsc.eta2=isat(diff(indicator2.res$pure$std.eta),sis=FALSE,iis=TRUE,ar=0)#for the level
  plot.my(out.bsc.eta2)
  out.bsc.zeta2=isat(indicator2.res$pure$std.zeta,ar=1,sis=FALSE,iis=TRUE)#for the slope
  plot.my(out.bsc.zeta2)
  out.bsc.omega2=isat(indicator2.res$pure$std.omega,sis=FALSE,iis=TRUE,ar=1)#for seasonality
  plot.my(out.bsc.omega2)

timeSeries3=indicator2
timeSeries3[is.na(timeSeries3)]=0
#TS=decompose(timeSeries)
TS2=decompose(timeSeries3,type="additive")
y=TS2$x-TS2$seasonal-TS2$trend

#testing for automatic order of ARIMA
auto.arima(timeSeries3)
auto.arima(y)

###IS with gets package
indicator2.IS.level=isat(timeSeries3, ar=1,t.pval=1/T, turbo=TRUE,iis=T,sis=T,tis=T,
                         vcov.type = "ordinary",info.method = 'aic')

start_time <- Sys.time()
indicator2.IS=isat(timeSeries3, ar=1,t.pval=1/T, turbo=TRUE,iis=T,sis=T,tis=T,
  vcov.type = "ordinary",info.method = 'aic')
end_time <- Sys.time()
print(end_time - start_time)

indicator2.tso.level=tso(y = timeSeries3, types = c("AO", "LS", "TC","SLS"),
    tsmethod = "stsm", args.tsmodel = list(model = "BSM")
)
indicator2.WOM.level=outlier(indicator2)

start_time <- Sys.time()
indicator2.tso=tso(y = timeSeries3, types = c("AO", "LS", "TC","SLS"),
                         tsmethod = "stsm", args.tsmodel = list(model = "BSM")
)
end_time <- Sys.time()
print(end_time - start_time)

start_time <- Sys.time()
indicator2.WOM=outlier(indicator2)
end_time <- Sys.time()
print(end_time - start_time)

## Breakdowns #age: 06-30 days
lts3.11.zoo=log(ts3.11.zoo)
lts3.11.zoo[lts3.11.zoo==-Inf]=0 #replace -inf by 0
indicator2.11 = ts(lts3.11.zoo,frequency =7) 
start_time <- Sys.time()
indicator2.11.res=isat_my(indicator2.11, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                             vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                             optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
end_time <- Sys.time()
print(end_time - start_time)

plot(indicator2.11.res)

#level
indicator2.11 = ts(ts3.11.zoo,frequency =7) 
indicator2.11.res.level=isat_my(indicator2.11, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                          vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                          optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
plot(indicator2.11.res.level)

#version without autocorrelated residuals
out.bsc.eps2.11=isat(indicator2.11.res$pure$std.eps,sis=FALSE,iis=TRUE,ar=0)
plot.my2(out.bsc.eps2.11)
out.bsc.eta2.11=isat(diff(indicator2.11.res$pure$std.eta),sis=FALSE,iis=TRUE)#for the level
plot(out.bsc.eta2.11)
out.bsc.zeta2.11=isat(diff(indicator2.11.res$pure$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
plot(out.bsc.zeta2.11)
out.bsc.omega2.11=isat(indicator2.11.res$pure$std.omega,sis=FALSE,iis=TRUE,ar=5)#for seasonality
plot(out.bsc.omega2.11)

timeSeries3.11=indicator2.11
timeSeries3.11[is.na(timeSeries3.11)]=0
#TS=decompose(timeSeries)
TS2.11=decompose(timeSeries3.11,type="additive")
y=TS2.11$x-TS2.11$seasonal-TS2.11$trend

#testing for automatic order of ARIMA
auto.arima(timeSeries3.11)
auto.arima(y)

###IS with gets package
indicator2.11.IS.level=isat(timeSeries3.11, ar=5,t.pval=1/T, turbo=TRUE,iis=T,sis=T,tis=T,
                   vcov.type = "ordinary",info.method = 'aic')

start_time <- Sys.time()
indicator2.11.IS=isat(timeSeries3.11, ar=5,t.pval=1/T, turbo=TRUE,iis=T,sis=T,tis=T,
                            vcov.type = "ordinary",info.method = 'aic')
end_time <- Sys.time()
print(end_time - start_time)

indicator2.11.tso.level=tso(y = timeSeries3.11, types = c("AO", "LS", "TC","SLS"),
                   tsmethod = "stsm", args.tsmodel = list(model = "BSM")
)

start_time <- Sys.time()
indicator2.11.tso=tso(y = timeSeries3.11, types = c("AO", "LS", "TC","SLS"),
                            tsmethod = "stsm", args.tsmodel = list(model = "BSM")
)
end_time <- Sys.time()
print(end_time - start_time)

indicator2.11.WOM.level=outlier(indicator2.11)
start_time <- Sys.time()
indicator2.11.WOM=outlier(indicator2.11)
end_time <- Sys.time()
print(end_time - start_time)
## Breakdowns #age: #age: 01-05 days
lts3.12.zoo=log(ts3.12.zoo)
lts3.12.zoo[lts3.12.zoo==-Inf]=0 #replace -inf by 0
indicator2.12 = ts(lts3.12.zoo,frequency =7) 
indicator2.12.res=isat_my(indicator2.12, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                          vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                          optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
plot(indicator2.12.res)

#level
indicator2.12 = ts(ts3.12.zoo,frequency =7) 
indicator2.12.res.level=isat_my(indicator2.12, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                          vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                          optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
plot.my(indicator2.12.res.level)

#version without autocorrelated residuals
out.bsc.eps2.12=isat(indicator2.12.res$pure$std.eps,sis=FALSE,iis=TRUE)
plot(out.bsc.eps2.12)
out.bsc.eta2.12=isat(indicator2.12.res$pure$std.eta,sis=FALSE,iis=TRUE)#for the level
plot(out.bsc.eta2.12)
out.bsc.zeta2.12=isat(indicator2.12.res$pure$std.zeta,sis=FALSE,iis=TRUE)#for the slope
plot(out.bsc.zeta2.12)
out.bsc.omega2.12=isat(indicator2.12.res$pure$std.omega,sis=FALSE,iis=TRUE,ar=5)#for seasonality
plot(out.bsc.omega2.12)

timeSeries3.12=indicator2.12
timeSeries3.12[is.na(timeSeries3.12)]=0
#TS=decompose(timeSeries)
TS2.12=decompose(timeSeries3.12,type="additive")
y=TS2.12$x-TS2.12$seasonal-TS2.12$trend

#testing for automatic order of ARIMA
auto.arima(timeSeries3.12)
auto.arima(y)

###IS with gets package
indicator2.12.IS.level=isat(timeSeries3.12, ar=1,t.pval=1/T, turbo=TRUE,iis=T,sis=T,tis=T,
                      vcov.type = "ordinary",info.method = 'aic')

start_time <- Sys.time()
indicator2.12.IS=isat(timeSeries3.12, ar=1,t.pval=1/T, turbo=TRUE,iis=T,sis=T,tis=T,
                            vcov.type = "ordinary",info.method = 'aic')
end_time <- Sys.time()
print(end_time - start_time)

indicator2.12.tso.level=tso(y = timeSeries3.12, types = c("AO", "LS", "TC","SLS"),
                      tsmethod = "stsm", args.tsmodel = list(model = "BSM")
)
start_time <- Sys.time()
indicator2.12.tso=tso(y = timeSeries3.12, types = c("AO", "LS", "TC","SLS"),
                            tsmethod = "stsm", args.tsmodel = list(model = "BSM")
)

end_time <- Sys.time()
print(end_time - start_time)

indicator2.12.WOM.level=outlier(indicator2.12)
indicator2.12.WOM=outlier(indicator2.12)

start_time <- Sys.time()
indicator2.12.WOM=outlier(indicator2.12)
end_time <- Sys.time()
print(end_time - start_time)

## Breakdowns  #age: 30-90 days
lts3.13.zoo=log(ts3.13.zoo)
lts3.13.zoo[lts3.13.zoo==-Inf]=0 #replace -inf by 0
indicator2.13 = ts(lts3.13.zoo,frequency =7) 
indicator2.13.res=isat_my(indicator2.13, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                          vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                          optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
plot.my2(indicator2.13.res)

#level
indicator2.13 = ts(ts3.13.zoo,frequency =7) 
indicator2.13.res.level=isat_my(indicator2.13, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                          vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                          optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
plot(indicator2.13.res.level)


#version without autocorrelated residuals
out.bsc.eps2.13=isat(indicator2.13.res$pure$std.eps,sis=FALSE,iis=TRUE)
plot(out.bsc.eps2.13)
out.bsc.eta2.13=isat(diff(indicator2.13.res$pure$std.eta),sis=FALSE,iis=TRUE)#for the level
plot(out.bsc.eta2.13)
out.bsc.zeta2.13=isat(diff(indicator2.13.res$pure$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
plot(out.bsc.zeta2.13)
out.bsc.omega2.13=isat(indicator2.13.res$pure$std.omega,sis=FALSE,iis=TRUE)#for seasonality
plot(out.bsc.omega2.13)

timeSeries3.13=indicator2.13
timeSeries3.13[is.na(timeSeries3.13)]=0
#TS=decompose(timeSeries)
TS2.13=decompose(timeSeries3.13,type="additive")
y=TS2.13$x-TS2.13$seasonal-TS2.13$trend

#testing for automatic order of ARIMA
auto.arima(timeSeries3.13)
auto.arima(y)

###IS with gets package
indicator2.13.IS.level=isat(timeSeries3.13, ar=5,t.pval=1/T, turbo=TRUE,iis=TRUE,
                      sis=TRUE,tis=TRUE,
                      vcov.type = "ordinary",info.method = 'aic')

start_time <- Sys.time()
indicator2.13.IS=isat(timeSeries3.13, ar=5,t.pval=1/T, turbo=TRUE,iis=TRUE,
                      sis=TRUE,tis=TRUE,
                      vcov.type = "ordinary",info.method = 'aic')
plot.my(indicator2.13.IS)
end_time <- Sys.time()
print(end_time - start_time)

plot.my(indicator2.13.IS)
indicator2.13.tso.level=tso(y = timeSeries3.13, types = c("AO", "LS", "TC","SLS"),
                      tsmethod = "stsm", args.tsmodel = list(model = "BSM")
)

start_time <- Sys.time()
indicator2.13.tso=tso(y = timeSeries3.13, types = c("AO", "LS", "TC","SLS"),
                            tsmethod = "stsm", args.tsmodel = list(model = "BSM")
)
end_time <- Sys.time()
print(end_time - start_time)

indicator2.13.WOM.level=outlier(indicator2.13)
indicator2.13.WOM=outlier(indicator2.13)

#########################
## Indicator 3
lts4.zoo=log(ts4.zoo)
lts4.zoo[lts4.zoo==-Inf]=0 #replace -inf by 0
indicator4 = ts(lts4.zoo,frequency =7) 
indicator4.res=isat_my(indicator4, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                          vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                          optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
plot.my(indicator4.res)

#level
indicator4 = ts(ts4.zoo,frequency =7) 
indicator4.res.level=isat_my(indicator4, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                       vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                       optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
#indicator4.res.level=indicator4.res
plot(indicator4.res.level)

#version without autocorrelated residuals
out.bsc.eps4=isat(indicator4.res$pure$std.eps,sis=FALSE,iis=TRUE)
plot(out.bsc.eps4)
out.bsc.eta4=isat(diff(indicator4.res$pure$std.eta),sis=FALSE,iis=TRUE)#for the level
plot(out.bsc.eta4)
out.bsc.zeta4=isat(diff(indicator4.res$pure$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
plot(out.bsc.zeta4)
out.bsc.omega4=isat(indicator4.res$pure$std.omega,sis=FALSE,iis=TRUE)#for seasonality
plot(out.bsc.omega4)

#version without autocorrelated residuals
out.bsc.eps4=isat(indicator4.res.level$pure$std.eps,sis=FALSE,iis=TRUE,ar=1)
plot(out.bsc.eps4)
out.bsc.eta4=isat(diff(indicator4.res.level$pure$std.eta),sis=FALSE,iis=TRUE)#for the level
plot(out.bsc.eta4)
out.bsc.zeta4=isat(diff(indicator4.res.level$pure$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
plot(out.bsc.zeta4)
out.bsc.omega4=isat(indicator4.res.level$pure$std.omega,sis=FALSE,iis=TRUE)#for seasonality
plot(out.bsc.omega4)

timeSeries4=indicator4
timeSeries4[is.na(timeSeries4)]=0
#require(imputeTS)
#timeSeries4=na.interpolation(timeSeries4, option = "linear")

#TS=decompose(timeSeries)
TS4=decompose(timeSeries4,type="additive")
y=TS4$x-TS4$seasonal-TS4$trend

#testing for automatic order of ARIMA
auto.arima(timeSeries4)
auto.arima(y)

###IS with gets package
indicator4.IS.level=isat(timeSeries4, ar=1,t.pval=1/T, turbo=TRUE,iis=TRUE,
                      sis=TRUE,tis=TRUE,
                      vcov.type = "ordinary",info.method = 'aic')
start_time <- Sys.time()
indicator4.IS=isat(timeSeries4, ar=5,t.pval=1/T, turbo=TRUE,iis=TRUE,
                   sis=TRUE,tis=TRUE,
                   vcov.type = "ordinary",info.method = 'aic')
end_time <- Sys.time()
print(end_time - start_time)

plot.my(indicator4.IS)
indicator4.tso.level=tso(y = timeSeries4, types = c("AO", "LS", "TC","SLS"),
                      tsmethod = "stsm", args.tsmodel = list(model = "BSM")
)

start_time <- Sys.time()
indicator4.tso=tso(y = timeSeries4, types = c("AO", "LS", "TC","SLS"),
                         tsmethod = "stsm", args.tsmodel = list(model = "BSM")
)
end_time <- Sys.time()
print(end_time - start_time)

indicator4.WOM.level=outlier(indicator4)

start_time <- Sys.time()
indicator4.WOM=outlier(indicator4)
end_time <- Sys.time()
print(end_time - start_time)

#########################
## Indicator 3
lts5.zoo=log(ts5.zoo)
lts5.zoo[lts5.zoo==-Inf]=0 #replace -inf by 0
indicator5 = ts(lts5.zoo,frequency =7) 
start_time <- Sys.time()
indicator5.res=isat_my(indicator5, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                       vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                       optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
end_time <- Sys.time()
print(end_time - start_time)

plot.my(indicator5.res)

#level
indicator5 = ts(ts5.zoo,frequency =7) 
indicator5.res.level=isat_my(indicator5, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                             vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                             optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
#indicator5.res.level=indicator5.res
plot(indicator5.res.level)

#version without autocorrelated residuals
out.bsc.eps5=isat(indicator5.res$pure$std.eps,sis=FALSE,iis=TRUE)
plot(out.bsc.eps5)
out.bsc.eta5=isat(diff(indicator5.res$pure$std.eta),sis=FALSE,iis=TRUE)#for the level
plot(out.bsc.eta5)
out.bsc.zeta5=isat(diff(indicator5.res$pure$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
plot(out.bsc.zeta5)
out.bsc.omega5=isat(indicator5.res$pure$std.omega,sis=FALSE,iis=TRUE)#for seasonality
plot(out.bsc.omega5)

#version without autocorrelated residuals
out.bsc.eps5=isat(indicator5.res.level$pure$std.eps,sis=FALSE,iis=TRUE,ar=1)
plot.my2(out.bsc.eps5)
out.bsc.eta5=isat((indicator5.res.level$pure$std.eta),sis=FALSE,iis=TRUE)#for the level
plot.my(out.bsc.eta5)
out.bsc.zeta5=isat(diff(indicator5.res.level$pure$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
plot.my(out.bsc.zeta5)
out.bsc.omega5=isat(indicator5.res.level$pure$std.omega,sis=FALSE,iis=TRUE)#for seasonality
plot.my(out.bsc.omega5)

timeSeries5=indicator5
timeSeries5[is.na(timeSeries5)]=0
#TS=decompose(timeSeries)
TS5=decompose(timeSeries5,type="additive")
y=TS5$x-TS5$seasonal-TS5$trend

#testing for automatic order of ARIMA
auto.arima(timeSeries5)
auto.arima(y)

###IS with gets package
indicator5.IS.level=isat(timeSeries5, ar=5,t.pval=1/T, turbo=TRUE,iis=TRUE,
                   sis=TRUE,tis=TRUE,
                   vcov.type = "ordinary",info.method = 'aic')
plot.my(indicator5.IS)
indicator5.tso.level=tso(y = timeSeries4, types = c("AO", "LS", "TC","SLS"),
                   tsmethod = "stsm", args.tsmodel = list(model = "BSM")
)

start_time <- Sys.time()
indicator5.IS=isat(timeSeries5, ar=5,t.pval=1/T, turbo=TRUE,iis=TRUE,
                   sis=TRUE,tis=TRUE,
                   vcov.type = "ordinary",info.method = 'aic')
end_time <- Sys.time()
print(end_time - start_time)

start_time <- Sys.time()
indicator5.tso=tso(y = timeSeries4, types = c("AO", "LS", "TC","SLS"),
                         tsmethod = "stsm", args.tsmodel = list(model = "BSM")
)
end_time <- Sys.time()
print(end_time - start_time)

indicator5.WOM.level=outlier(indicator5)

start_time <- Sys.time()
indicator5.WOM=outlier(indicator5)
end_time <- Sys.time()
print(end_time - start_time)

results=list(indicator1.res,indicator1.IS,indicator1.tso,indicator1.WOM,
             indicator1.11.res,indicator1.11.IS,indicator1.11.tso,indicator1.11.WOM,
             indicator1.21.res,indicator1.21.IS,indicator1.21.tso,indicator1.21.WOM,
             indicator1.22.res,indicator1.22.IS,indicator1.22.tso,indicator1.22.WOM,
             indicator2.res,indicator2.IS,indicator2.tso,indicator2.WOM,
             indicator2.11.res,indicator2.11.IS,indicator2.11.tso,indicator2.11.WOM,
             indicator2.12.res,indicator2.12.IS,indicator2.12.tso,indicator2.12.WOM,
             indicator2.13.res,indicator2.13.IS,indicator2.13.tso,indicator2.13.WOM,
            indicator4.IS,indicator4.tso,indicator4.WOM,
             indicator5.res,indicator5.IS,indicator5.tso,indicator5.WOM)
save(results,file="indicator_results.RData")

results.level=list(indicator1.res.level,indicator1.IS.level,indicator1.tso.level,indicator1.WOM.level,
             indicator1.11.res.level,indicator1.11.IS.level,indicator1.11.tso.level,indicator1.11.WOM.level,
             indicator1.21.res.level,indicator1.21.IS.level,indicator1.21.tso.level,indicator1.21.WOM.level,
             indicator1.22.res.level,indicator1.22.IS.level,indicator1.22.tso.level,indicator1.22.WOM.level,
             indicator2.res.level,indicator2.IS.level,indicator2.tso.level,indicator2.WOM.level,
             indicator2.11.res.level,indicator2.11.IS.level,indicator2.11.tso.level,indicator2.11.WOM.level,
             indicator2.12.res.level,indicator2.12.IS.level,indicator2.12.tso.level,indicator2.12.WOM.level,
             indicator2.13.res.level,indicator2.13.IS.level,indicator2.13.tso.level,indicator2.13.WOM.level,
             indicator4.res,indicator4.res.level,indicator4.IS.level,indicator4.tso.level,indicator4.WOM.level,
             indicator5.res.level,indicator5.IS.level,indicator5.tso.level,indicator5.WOM.level)
save(results.level,file="indicator_results_level")
######## check the stationarity of TS
y[is.na(y)]=0 #only weekly trend included
#y2[is.na(y2)]=0 #weekly + monthly trend included

acf(y)
acf(y2)
acf(timeSeries)
pacf(y)
pacf(y2)
pacf(timeSeries) # from pacf it is obvious there is weekly seasonality and first order
#autocorrelation (rising prob from daily data)

#adf test
require(tseries)
adf.test(y) #no unit root; linear trend included, but lag 5 => mighht be seasonality
adf.test(y2) #no unit root; linear trend included, but lag 5 => mighht be seasonality
adf.test(timeSeries)

#kpss test
kpss.test(y)#no unit root
kpss.test(y2)#no unit root
kpss.test(timeSeries,"Trend") # time series without seasonality adjustment
