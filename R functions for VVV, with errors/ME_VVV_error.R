


### manual EM algorithm for VVV with estimation error

ME.VVV.err = function(data, z, err){
  
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)

  k = 2 # keeps track of number of iteration

  piconst = n*p*log(2*pi)/2
  
  loglikelihood = NA
  parameters = list()
  zhat = matrix(0,n,G)

  llike = rep(0, 10000) # set convergence criterion
  llike[2] = 1e-5
 
  tol = 1e-7
  # while loop for iteration:
  while(abs(llike[k]-llike[k-1]) > tol){
    
    print(paste("iteration =",k-1)) # prints number of iterations
  
    thetahat = MstepVVV.err(z, data, err) # M-step
    
    temp = EstepVVV.err(thetahat, data, err) # E-step
  
    zhat = temp[[1]] # membership estimates
    
    parameters = temp[[2]] # parameter estimates
    
    loglikelihood = temp[[3]] # records log likelihood
    
    parsigma = temp[[4]]
      
    z = zhat # update membership matrix
    
    llike[k+1] = loglikelihood # update log likelihood of complete data
    
    print(loglikelihood) # prints current evaulated complete data log-likelihood
      
    k = k+1 # increment k
  }
  
  # edit output so it's basically consistent with meVVV() from MCLUST:
  out = list(modelname="VVV with est error", n=n, d=p, G=G, z=zhat,
    parameters=parameters, loglik=loglikelihood, iteration=k-1)
  
  return(out)
  

}