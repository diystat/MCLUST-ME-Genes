



### manual EM algorithm for VVV with estimation error

ME.VVV.err = function(data, z, err){
  
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)

  k = 2 # keeps track of number of iteration

  
  loglikelihood = NA
  parameters = list()

  llike = rep(0, 10000)
  llike[2] = 0.0011
  

  # while loop for iteration:
  while(abs(llike[k]-llike[k-1])>0.001){
    
    print(paste("iteration =",(k-1)))
  
    thetahat = MstepVVV.err(z, data, err) # M-step
    
    temp = EstepVVV.err(thetahat, data, err) # E-step
  
    zhat = temp[[1]] # membership estimates
    
    parameters = temp[[2]] # parameter estimates
    
    loglikelihood = temp[[3]] # records log likelihood
      
    z = zhat # update membership matrix
    
    llike[k+1] = loglikelihood # update log likelihood
      
    k = k+1 # increment k
  }
  
  # edit output so it's consistent with meEEE() from MCLUST:
  out = list(modelname="VVV with est error", n=n, d=p, G=G, z=zhat,
    parameters=parameters, loglik=loglikelihood, iteration=k-1)
  
  return(out)
  

}