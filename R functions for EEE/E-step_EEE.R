
### Suppose "param" is a list containing all estimated parameters from M-step.
###
### Now we attempt to write the E-step for EEE structure:

EstepEEE = function(param, data){
  n = nrow(data)
  p = ncol(data)
  G = ncol(param[[1]])
  
  muhat = param[[1]]
  tauhat = param[[2]]
  sigmahat = param[[3]]
  
  
  # set initial z matrix:
  zhat = matrix(rep(0,n*G),nrow=n)
  
  piconst = (2*pi)^(-p/2)
  
  denom = rep(0,n)
  for (i in 1:n){
    for (j in 1:G){
      denom[i] = denom[i] + tauhat[j]*piconst*(det(sigmahat))^(-1/2)*exp(-1/2*t(data[i,]-muhat[,j])%*%solve(sigmahat)%*%(data[i,]-muhat[,j]))
    }
  }
  
  for(i in 1:n){
    for(k in 1:G){
      num = tauhat[k]*piconst*(det(sigmahat))^(-1/2)*exp(-1/2*t(data[i,]-muhat[,k])%*%solve(sigmahat)%*%(data[i,]-muhat[,k]))
      zhat[i,k] = num/denom[i]
    }
  }
  
  parameters = list(muhat=muhat, tauhat=tauhat, sigmahat=sigmahat)
  
  loglikelihood = log.lik.eee(param, data)
  
  out = list(zhat, parameters, loglikelihood)
  
  return(out) # output estimated z's for M-step
  
}




#----------------------------------------------------------------------------------------------#



## debugging
EstepEEE.test.fun = function(){
  # source("Estep_EEE.R")
  n = 100
  p = 3
  G = 4
  
  set.seed(999)
  data = matrix(rnorm(n*p),n,p) # data
  
  set.seed(998)
  mean = matrix(rnorm(p*G),p,G) # estimated means
  
  prop = matrix(0, n, G) # estimated mixing proportions
  prop[,1] = 1
  
  sigma = diag(1:p) # estimated cov matrix
  
  param = list(mean,prop,sigma) # put everything in a list
  
  zhat = EstepEEE(param,data) # evaluate
  
  zhat
  
  rowSums(zhat)
  
}
