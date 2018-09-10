start_time <- Sys.time()
indicator4.res.short=isat_my_short(indicator4, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                       vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                       optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
#plot.my(indicator4.res)

#version without autocorrelated residuals
out.bsc.eps4=isat(indicator4.res.short$std.eps,sis=FALSE,iis=TRUE)
#plot(out.bsc.eps4)
out.bsc.eta4=isat(diff(indicator4.res.short$std.eta),sis=FALSE,iis=TRUE)#for the level
#plot(out.bsc.eta4)
out.bsc.zeta4=isat(diff(indicator4.res.short$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
#plot(out.bsc.zeta4)
out.bsc.omega4=isat(indicator4.res.short$std.omega,sis=FALSE,iis=TRUE)#for seasonality
#plot(out.bsc.omega4)
end_time <- Sys.time()
print(end_time - start_time)



start_time <- Sys.time()
indicator5.res.short=isat_my_short(indicator5, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                                   vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                                   optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
#plot.my(indicator5.res)

#version without autocorrelated residuals
out.bsc.eps5=isat(indicator5.res.short$pure$std.eps,sis=FALSE,iis=TRUE)
#plot(out.bsc.eps4)
out.bsc.eta5=isat(diff(indicator5.res.short$std.eta),sis=FALSE,iis=TRUE)#for the level
#plot(out.bsc.eta4)
out.bsc.zeta5=isat(diff(indicator5.res.short$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
#plot(out.bsc.zeta4)
out.bsc.omega5=isat(indicator5.res.short$std.omega,sis=FALSE,iis=TRUE)#for seasonality
#plot(out.bsc.omega4)
end_time <- Sys.time()
print(end_time - start_time)

start_time <- Sys.time()
indicator1.res.short=isat_my_short(indicator1, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                                   vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                                   optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
#plot.my(indicator5.res)

#version without autocorrelated residuals
out.bsc.eps1=isat(indicator1.res.short$std.eps,sis=FALSE,iis=TRUE)
#plot(out.bsc.eps4)
out.bsc.eta1=isat(diff(indicator1.res.short$std.eta),sis=FALSE,iis=TRUE)#for the level
#plot(out.bsc.eta4)
out.bsc.zeta1=isat(diff(indicator1.res.short$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
#plot(out.bsc.zeta4)
out.bsc.omega1=isat(indicator1.res.short$std.omega,sis=FALSE,iis=TRUE)#for seasonality
#plot(out.bsc.omega4)
end_time <- Sys.time()
print(end_time - start_time)

start_time <- Sys.time()
indicator2.res.short=isat_my_short(indicator2, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                                   vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                                   optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
#plot.my(indicator5.res)

#version without autocorrelated residuals
out.bsc.eps2=isat(indicator2.res.short$std.eps,sis=FALSE,iis=TRUE)
#plot(out.bsc.eps4)
out.bsc.eta2=isat(diff(indicator2.res.short$std.eta),sis=FALSE,iis=TRUE)#for the level
#plot(out.bsc.eta4)
out.bsc.zeta2=isat(diff(indicator2.res.short$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
#plot(out.bsc.zeta4)
out.bsc.omega2=isat(indicator2.res.short$std.omega,sis=FALSE,iis=TRUE)#for seasonality
#plot(out.bsc.omega4)
end_time <- Sys.time()
print(end_time - start_time)

start_time <- Sys.time()
indicator2.11.res.short=isat_my_short(indicator2.11, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                                   vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                                   optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
#plot.my(indicator5.res)

#version without autocorrelated residuals
out.bsc.eps2.11=isat(indicator2.11.res.short$std.eps,sis=FALSE,iis=TRUE)
#plot(out.bsc.eps4)
out.bsc.eta2.11=isat(diff(indicator2.11.res.short$std.eta),sis=FALSE,iis=TRUE)#for the level
#plot(out.bsc.eta4)
out.bsc.zeta2.11=isat(diff(indicator2.11.res.short$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
#plot(out.bsc.zeta4)
out.bsc.omega2.11=isat(indicator2.11.res.short$std.omega,sis=FALSE,iis=TRUE)#for seasonality
#plot(out.bsc.omega4)
end_time <- Sys.time()
print(end_time - start_time)

start_time <- Sys.time()
indicator2.12.res.short=isat_my_short(indicator2.12, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                                      vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                                      optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
#plot.my(indicator5.res)

#version without autocorrelated residuals
out.bsc.eps2.12=isat(indicator2.12.res.short$std.eps,sis=FALSE,iis=TRUE)
#plot(out.bsc.eps4)
out.bsc.eta2.12=isat(diff(indicator2.12.res.short$std.eta),sis=FALSE,iis=TRUE)#for the level
#plot(out.bsc.eta4)
out.bsc.zeta2.12=isat(diff(indicator2.12.res.short$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
#plot(out.bsc.zeta4)
out.bsc.omega2.12=isat(indicator2.12.res.short$std.omega,sis=FALSE,iis=TRUE)#for seasonality
#plot(out.bsc.omega4)
end_time <- Sys.time()
print(end_time - start_time)

start_time <- Sys.time()
indicator2.13.res.short=isat_my_short(indicator2.13, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                                      vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                                      optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
#plot.my(indicator5.res)

#version without autocorrelated residuals
out.bsc.eps2.13=isat(indicator2.13.res.short$std.eps,sis=FALSE,iis=TRUE)
#plot(out.bsc.eps4)
out.bsc.eta2.13=isat(diff(indicator2.13.res.short$std.eta),sis=FALSE,iis=TRUE)#for the level
#plot(out.bsc.eta4)
out.bsc.zeta2.13=isat(diff(indicator2.13.res.short$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
#plot(out.bsc.zeta4)
out.bsc.omega2.13=isat(indicator2.13.res.short$std.omega,sis=FALSE,iis=TRUE)#for seasonality
#plot(out.bsc.omega4)
end_time <- Sys.time()
print(end_time - start_time)

start_time <- Sys.time()
indicator1.11.res.short=isat_my_short(indicator1.11, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                                      vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                                      optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
#plot.my(indicator5.res)

#version without autocorrelated residuals
out.bsc.eps1.11=isat(indicator1.11.res.short$std.eps,sis=FALSE,iis=TRUE)
#plot(out.bsc.eps4)
out.bsc.eta1.11=isat(diff(indicator1.11.res.short$std.eta),sis=FALSE,iis=TRUE)#for the level
#plot(out.bsc.eta4)
out.bsc.zeta1.11=isat(diff(indicator1.11.res.short$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
#plot(out.bsc.zeta4)
out.bsc.omega1.11=isat(indicator1.11.res.short$std.omega,sis=FALSE,iis=TRUE)#for seasonality
#plot(out.bsc.omega4)
end_time <- Sys.time()
print(end_time - start_time)

start_time <- Sys.time()
indicator1.21.res.short=isat_my_short(indicator1.21, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                                      vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                                      optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
#plot.my(indicator5.res)

#version without autocorrelated residuals
out.bsc.eps1.21=isat(indicator1.21.res.short$std.eps,sis=FALSE,iis=TRUE)
#plot(out.bsc.eps4)
out.bsc.eta1.21=isat(diff(indicator1.21.res.short$std.eta),sis=FALSE,iis=TRUE)#for the level
#plot(out.bsc.eta4)
out.bsc.zeta1.21isat(diff(indicator1.21.res.short$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
#plot(out.bsc.zeta4)
out.bsc.omega1.21=isat(indicator1.21.res.short$std.omega,sis=FALSE,iis=TRUE)#for seasonality
#plot(out.bsc.omega4)
end_time <- Sys.time()
print(end_time - start_time)


start_time <- Sys.time()
indicator1.22.res.short=isat_my_short(indicator1.22, t.pval=1/T, model.ssmod = "llm+seas",turbo=TRUE,
                                      vcov.type = "ordinary",info.method = 'aic', gr.kfs ="numerical", 
                                      optim.kfs = list(lower = 0, upper = Inf, hessian = TRUE))
#plot.my(indicator5.res)

#version without autocorrelated residuals
out.bsc.eps1.22=isat(indicator1.22.res.short$std.eps,sis=FALSE,iis=TRUE)
#plot(out.bsc.eps4)
out.bsc.eta1.22=isat(diff(indicator1.22.res.short$std.eta),sis=FALSE,iis=TRUE)#for the level
#plot(out.bsc.eta4)
out.bsc.zeta1.22=isat(diff(indicator1.22.res.short$std.zeta,2),sis=FALSE,iis=TRUE)#for the slope
#plot(out.bsc.zeta4)
out.bsc.omega1.22=isat(indicator1.22.res.short$std.omega,sis=FALSE,iis=TRUE)#for seasonality
#plot(out.bsc.omega4)
end_time <- Sys.time()
print(end_time - start_time)