


### manual EM algorithm for VVV, starting with an M-step

ME.VVV = function(data, z){
  
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)

  k = 2 # keeps track of number of iteration

  dim = c(n,G,10000)
  arr = array(dim=dim)
  loglikelihood = NA
  parameters = list()

  arr[,,1] = matrix(rep(0,n*G),nrow=n)
  arr[,,2] = z

  # num = .Machine$double.xmin

  # while loop for iteration:
  while(sum(abs(arr[,,k]-arr[,,(k-1)]))>0.001){
  
    thetahat = MstepVVV(z, data) # M-step
    
    temp = EstepVVV(thetahat, data) # E-step
  
    zhat = temp[[1]] # membership estimates
    
    parameters = temp[[2]] # parameter estimates
    
    loglikelihood = temp[[3]] # records log likelihood
      
    z = arr[,,(k+1)] = zhat # update membership matrix
      
    k = k+1 # increment k
  }
  
  # edit output so it's consistent with meEEE() from MCLUST:
  out = list(modelname="VVV", n=n, d=p, G=G, z=zhat,
    parameters=parameters, loglik=loglikelihood, iteration=k-1)
  
  return(out)
  

}