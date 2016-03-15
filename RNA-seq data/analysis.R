  rm(list=ls());
  dd = "2015-11-06";
  file.out = sprintf("%s.Coaker.Rdata", dd);
  print(load(file.out));

  obs = full$beta[top, 1:5];
  errary = v.beta[,,top];
  betasd = sd.beta[top,];
  
  betasd
  
  # Check order of magnitude of off-diagonal elements
  dim(errary)
  offdiag = numeric()
  for(i in 1:1000){
    offdiag = c(offdiag,gdata::lowerTriangle(errary[,,i],diag=F))
  }
  offdiag
  order = floor(log10(abs(offdiag)))
  hist(order)
  # Cannot say that the covariance matrices are diagonal (have order -5 & -6)
  
  
  library(mclust)
  res = Mclust(obs,modelNames="VVV")
  summary(res)
  plot(res,what="classification")
  
  cov(obs)
  cor(obs)
  