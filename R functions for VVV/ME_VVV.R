


### manual EM algorithm for VVV, starting with an M-step

ME.VVV = function(data, z){
  
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)

  k = 2 # keeps track of number of iteration

  loglikelihood = NA
  parameters = list()

  llike = rep(0,10000)
  llike[2] = 0.01
  
  tol = 1e-6

  # while loop for iteration:
  while(abs(llike[k]-llike[k-1])>tol){
    
    print(paste("iteration = ",k-1))
  
    thetahat = MstepVVV(z, data) # M-step
    
    temp = EstepVVV(thetahat, data) # E-step
  
    zhat = temp[[1]] # membership estimates
    
    parameters = temp[[2]] # parameter estimates
    
    loglikelihood = temp[[3]] # records log likelihood
    
    print(loglikelihood)
      
    z = zhat # update membership matrix
    
    llike[k+1] = loglikelihood
      
    k = k+1 # increment k
  }
  
  # edit output so it's consistent with meEEE() from MCLUST:
  out = list(modelname="VVV", n=n, d=p, G=G, z=zhat,
    parameters=parameters, loglik=loglikelihood, iteration=k-1)
  
  return(out)
  

}