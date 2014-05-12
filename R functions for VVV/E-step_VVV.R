

### Suppose "param" is a list containing all estimated parameters from M-step.
###
### Now we attempt to write the E-step for VVV structure:

EstepVVV = function(param, data){
  n = nrow(data)
  p = ncol(data)
  G = ncol(param[[1]])
  
  muhat = param[[1]]
  tauhat = param[[2]]
  sigmahat = param[[3]] # here sigmahat is an array consisting of G p*p cov matrices
  
  
  # set initial z matrix:
  zhat = matrix(rep(0,n*G),nrow=n)
  
  piconst = (2*pi)^(-p/2)
  
  denom = rep(0,n)
  for (i in 1:n){
    for (k in 1:G){
      denom[i] = denom[i] + tauhat[k]*piconst*(det(sigmahat[,,k]))^(-1/2)*exp(-1/2*t(data[i,]-muhat[,k])%*%solve(sigmahat[,,k])%*%(data[i,]-muhat[,k]))
    }
  }
  
  for(i in 1:n){
    for(k in 1:G){
      num = tauhat[k]*piconst*(det(sigmahat[,,k]))^(-1/2)*exp(-1/2*t(data[i,]-muhat[,k])%*%solve(sigmahat[,,k])%*%(data[i,]-muhat[,k]))
      zhat[i,k] = num/denom[i]
    }
  }
  
  parameters = list(muhat=muhat, tauhat=tauhat, sigmahat=sigmahat)
  
  loglikelihood = log.lik.vvv(param, data)
  
  out = list(zhat, parameters, loglikelihood)
  
  return(out) # output estimated z's, parameter estimates and estimated log-likelihood
  
}

