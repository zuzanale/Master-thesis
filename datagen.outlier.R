datagen.outlier=function (n, model = list(), SigmaEV, labels, n0 = 20, freq = 1, 
          old.version = FALSE, AOplus=NULL,AOplusmag=NULL, AOminus=NULL,AOminusmag=NULL,LSplus=NULL,LSplusmag=NULL
          ,LSminus=NULL, LSminusmag=NULL, 
          TCplus=NULL, TCminus=NULL, SAOminusmag=NULL,TCplusmag=NULL,TCminusmag=NULL,
          SAOplus=NULL, SAOminus=NULL,SAOplusmag=NULL) 
  
{ #########################################
  # A0 is a vector of indexes for additive outliers, 
  #LS is a vector of indexes for level shifts, 
  #TS iis a vector of indexes for trend changes
  #SAO iis a vector of indexes for seasonality additive outliers.
  #########################################
  
  n <- n + n0
  y <- rep(NA, n)
  a <- matrix(nrow = n, ncol = ncol(model$Z))
  if (is.null(model$H)) {
    eps1 <- rep(0, n)
  }
  else {eps1 <- rnorm(n, sd = sqrt(model$H[1]))
  }
  if (!is.null(AOplus)){
    for (j in AOplus){
      j=j+n0
      if (!is.null(model$H)){
        # 4 is a benchmark magnitude of an outlier but we should make different senarious ())
        eps1[j]= eps1[j]+ sqrt(model$H)*AOplusmag
      }else{
        eps1[j]= eps1[j]+ max(eps1[j])*AOplusmag
      }
    }
  }
  #same for negative outliers
  if (!is.null(AOminus)){
    for (j in AOminus){
      j=j+n0
      if (!is.null(model$H)){
        # 4 is a benchmark magnitude of an outlier but we should make different senarious ())
        eps1[j]= eps1[j]- sqrt(model$H)*AOminusmag
      }else{
        eps1[j]= eps1[j]- min(eps1[j])*AOminusmag
      }
    }
  }
  if (old.version) {
    Meps <- mvtnorm::rmvnorm(n, mean = rep(0, nrow(model$Q)), 
                             sigma = model$Q, method = "eigen", pre0.9_9994 = TRUE)
  }
  else {
    if (missing(SigmaEV)) 
      SigmaEV <- eigen(model$Q)
    eps <- rnorm(n * ncol(SigmaEV$vectors))
    Meps <- matrix(eps, ncol = n, byrow = TRUE)
    Meps <- SigmaEV$vectors %*% diag(sqrt(SigmaEV$values)) %*% 
      Meps
    Meps <- t(Meps)
  }
  if (!is.null(model$a0)) {
    a0 <- model$a0
  }
  else a0 <- rep(0, ncol(model$Z))
    #here to create outliers for states
    #LS
    if (!is.null(LSplus)){
      for (j in LSplus){
        j=j+n0
        if (!is.null(diag(model$Q)[1])){
          # 5 is a benchmark magnitude of an outlier but we should make different senarious ())
        Meps[j,1]=Meps[j,1] + sqrt(diag(model$Q)[1])*LSplusmag
        }else{
          Meps[j,1]=Meps[j,1] + max(Meps[j,1])*LSplusmag
        }
      }
    }
    #same for negative outliers
    if (!is.null(LSminus)){
      for (j in LSminus){
        j=j+n0
        if (!is.null(diag(model$Q)[1])){
          # 5 is a benchmark magnitude of an outlier but we should make different senarious ())
          Meps[j,1]=Meps[j,1] - sqrt(diag(model$Q)[1])*LSminusmag
        }else{
          Meps[j,1]=Meps[j,1] - min(Meps[j,1])*LSminusmag
        }
      }
    }
    
    if (!is.null(TCplus)){
      for (j in TCplus){
        j=j+n0
        if (!is.null(diag(model$Q)[2])){
          # 5 is a benchmark magnitude of an outlier but we should make different senarious ())
          Meps[j,2]=Meps[j,2] + sqrt(diag(model$Q)[2])*TCplusmag
        }else{
          Meps[j,2]=Meps[j,2] + max(Meps[j,2])*TCplusmag
        }
      }
    }
    #same for negative outliers
    if (!is.null(TCminus)){
      for (j in TCminus){
        j=j+n0
        if (!is.null(diag(model$Q)[2])){
          # 5 is a benchmark magnitude of an outlier but we should make different senarious ())
          Meps[j,2]=Meps[j,2] - sqrt(diag(model$Q)[2])*TCminusmag
        }else{
          Meps[j,2]=Meps[j,2] - min(Meps[j,2])*TCminusmag
        }
      }
    }
     
    if (!is.null(SAOplus)){
      for (j in SAOplus){
        j=j+n0
        if (!is.null(diag(model$Q)[3])){
          # 5 is a benchmark magnitude of an outlier but we should make different senarious ())
          Meps[j,3]=Meps[j,3] + sqrt(diag(model$Q)[3])*SAOplusmag
        }else{
          Meps[j,3]=Meps[j,3] + max(Meps[j,3])*SAOplusmag
        }
      }
    }
    #same for negative outliers
    if (!is.null(SAOminus)){
      for (j in SAOminus){
        j=j+n0
        if (!is.null(diag(model$Q)[3])){
          # 5 is a benchmark magnitude of an outlier but we should make different senarious ())
          Meps[j,3]=Meps[j,3] - sqrt(diag(model$Q)[3])*SAOminusmag
        }else{
          Meps[j,3]=Meps[j,3] - min(Meps[j,3])*SAOminusmag
        }
      }
    }
  for (i in seq(n)) {
    a0 <- a[i, ] <- model$T %*% a0 + Meps[i, ]
    y[i] <- model$Z %*% a[i, ] + eps1[i]
  }
  if (n0 > 0) {
    y <- y[-seq(n0)]
    a <- a[-seq(n0), ]
  }
  y <- ts(y, frequency = freq)
  a <- ts(a, frequency = freq)
  if (is.matrix(a)) 
    a <- a[, diag(model$Q) != 0]
  if (!missing(labels)) {
    if (length(labels) == ncol(a)) {
      colnames(a) <- labels
    }
    else warning("The length of 'labels' is not equal to the number of", 
                 "  components with non-zero variance.")
  }
  list(data = y, components = a)
}