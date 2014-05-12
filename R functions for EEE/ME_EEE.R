

### manual EM algorithm for EEE, starting with an M-step

ME.EEE = function(data, z){
  
  # source("M-step_EEE.R")
  # source("E-step_EEE.R")
  # source("log likelihood.R")
  # source("W matrix.R")
  # source("objective function_EEE.R")
  
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)

  k = 2 # keeps track of number of iteration

  dim = c(n,G,1000)
  arr = array(dim=dim)
  onevec = rep(1,G)
  loglikelihood = NA
  parameters = list()

  arr[,,1] = matrix(rep(0,n*G),nrow=n)
  arr[,,2] = z

  # num = .Machine$double.xmin

  # while loop for iteration:
  while(sum(abs(arr[,,k]-arr[,,(k-1)]))>0.001){
  
    thetahat = MstepEEE(z, data) # M-step
    
    temp = EstepEEE(thetahat, data) # E-step
  
    zhat = temp[[1]] # membership estimates
    
    parameters = temp[[2]] # parameter estimates
    
    loglikelihood = temp[[3]] # records log likelihood
      
    z = arr[,,(k+1)] = zhat # update membership matrix
      
    k = k+1 # increment k
  }
  
  # edit output so it's consistent with meEEE() from MCLUST:
  out = list(modelname="EEE", n=n, d=p, G=G, z=zhat,
    parameters=parameters, loglik=loglikelihood, iteration=k-1)
  
  return(out)
  

}
  
