

### Suppose "param" is a list containing all estimated parameters from M-step.
###
### Now we attempt to write the E-step for VVV structure, with estimation error

EstepVVV.err = function(param, data, err){
  n = nrow(data)
  p = ncol(data)
  G = ncol(param[[1]])
  
  muhat = param[[1]]
  tauhat = param[[2]]
  sigmahat = param[[3]] # here sigmahat is an array consisting of G p*p cov matrices
  inv.sum = param[[4]]
  
  # set initial z matrix:
  zhat = matrix(rep(0,n*G),nrow=n)
    
  denom = rep(0,n)
  for (i in 1:n){
    for (k in 1:G){
      L = chol(inv.sum[,,k,i])
      temp = as.numeric(data[i,]-muhat[,k])
      denom[i] = denom[i] + tauhat[k]*(det(sigmahat[,,k]+err[,,i]))^(-1/2)*exp(-1/2*crossprod(L%*%temp))
    }
  }
  
  for(i in 1:n){
    for(k in 1:G){
      L = chol(inv.sum[,,k,i])
      temp = as.numeric(data[i,]-muhat[,k])
      num = tauhat[k]*(det(sigmahat[,,k]+err[,,i]))^(-1/2)*exp(-1/2*crossprod(L%*%temp))
      zhat[i,k] = num/denom[i]
    }
  }
  
  parameters = list(muhat=muhat, tauhat=tauhat, sigmahat=sigmahat, inv.sum=inv.sum)
  
  loglikelihood = log.lik.vvv.err(parameters, data, err) # observed log likelihood
  
  out = list(zhat, parameters, loglikelihood)
  
  return(out) # output estimated z's, parameter estimates and estimated log-likelihood
  
}
