

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
  parsigma = param[[4]] # vector of lower triangular elements of cholesky factors of cov matrices
  
  # set initial z matrix:
  zhat = matrix(rep(0,n*G),nrow=n)
  
  piconst = (2*pi)^(-p/2)
  
  inv.sum = array(0, dim=c(p,p,G,n))
  denom = rep(0,n)
  for (i in 1:n){
    for (k in 1:G){
      inv.sum[,,k,i] = solve(sigmahat[,,k]+err[,,i])
      L = chol(inv.sum[,,k,i])
      temp = data[i,]-muhat[,k]
      denom[i] = denom[i] + tauhat[k]*piconst*(det(sigmahat[,,k]+err[,,i]))^(-1/2)*exp(-1/2*crossprod(L%*%temp))
    }
  }
  
  for(i in 1:n){
    for(k in 1:G){
      L = chol(inv.sum[,,k,i])
      temp = data[i,]-muhat[,k]
      num = tauhat[k]*piconst*(det(sigmahat[,,k]+err[,,i]))^(-1/2)*exp(-1/2*crossprod(L%*%temp))
      zhat[i,k] = num/denom[i]
    }
  }
  
  parameters = list(muhat=muhat, tauhat=tauhat, sigmahat=sigmahat)
  
  loglikelihood = log.lik.vvv.err(parameters, data, err) # observed log likelihood
  
  out = list(zhat, parameters, loglikelihood, parsigma)
  
  return(out) # output estimated z's, parameter estimates and estimated log-likelihood
  
}
