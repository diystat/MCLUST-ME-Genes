
## Step-by-step iteration for MCLUST VVV algorithm:
stepwiseME = function(data, z, itmax=Inf){
  ## data---data matrix
  ## z---initial membership matrix
  
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)
  
  member = z # matrix for storing membership estimates
  center = matrix(0,p,G) # matrix for storing mean estimates

  k = 2 # keeps track of number of iteration

  loglikelihood = NA
  parameters = list()
  zhat = matrix(0,n,G)

  llike = rep(0, 1000) # set convergence criterion
  FLMAX = 1.7976931348623157e308
  llike[1] = FLMAX/2
  llike[2] = FLMAX

  tol = 1e-5
  
  # loop for iteration: 
  repeat{
    
    #print(paste("iteration =",k-1)) # prints number of iterations
    
    thetahat = mstepVVV(data, z)$parameters # M-step
    
    temp = estepVVV(data, thetahat) # E-step
  
    zhat = temp$z # membership estimates
    
    muhat = temp$parameters$mean # mean estimates
    
    parameters = temp$parameters
    
    loglikelihood = temp$loglik # records log likelihood
          
    z = zhat # update membership matrix
    
    member = rbind(member,zhat) # store membership estimates
    
    center = rbind(center,muhat) # store mean estimates
    
    llike[k+1] = loglikelihood # update log likelihood of observed data
    
    #print(loglikelihood) # prints current evaulated observed data log-likelihood
      
    k = k+1 # increment k
        
    delta = abs(llike[k-1]-llike[k])/(1+abs(llike[k]))
    
    it = k-2
    
    if(delta<tol || it>=itmax) break;
  }
  
  
  # edit output so it's basically consistent with meVVV() from MCLUST:
  out = list(modelname="VVV with est error", n=n, d=p, G=G, z=zhat, member=member, center=center,
    parameters=parameters, loglik=loglikelihood, iteration=k-2, likvec=llike[1:k][-(1:2)])
  
  return(out)
  
}