

meVVV.fc = function(data, z,itmax=Inf, mu1, mu2){
  # Argument "errstr" is user-specified error structure.
  # If errstr="identical", all errors are the same.
  # If errstr="cluster", errors are the same within each cluster.
  # If errstr="none", no constraints on error structure.
  # When errstr="none", "d" denotes the diagonal element, default = 1
  
  # Argument "lb" sets lower bound for diagonal elements of decomposed cov
  # matrices. Default set at 1e-3=0.001. Can be set to a larger number to
  # avoid singularity issues.
  
  source("E-step_VVV_error.R")
  source("M-step_VVV_error.R")
  source("log likelihood_VVV_error.R")
  source("objective function_VVV_error.R")
  source("W_k matrix.R")
  source("inipar.R")
  
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)
  
  member = z # matrix for storing membership estimates
  center = matrix(0,p,G) # matrix for storing mean estimates

  k = 2 # keeps track of number of iteration

  piconst = n*p*log(2*pi)/2
  
  loglikelihood = NA
  parameters = list()
  zhat = matrix(0,n,G)

  llike = rep(0, 1000) # set convergence criterion
  FLMAX = 1.7976931348623157e308
  llike[1] = FLMAX/2
  llike[2] = FLMAX
 
  tol = 1e-5
  
  mu.fix = cbind(mu1,mu2)
  
  
  # while loop for iteration:
  repeat{
    
    print(paste("iteration =",k-1)) # prints number of iterations
      
    thetahat = mstepVVV(data, z) # M-step
    thetahat$parameters$mean = mu.fix
    
    temp = estepVVV(data, thetahat$parameters) # E-step
  
    zhat = temp$z # membership estimates
    
    parameters = temp$parameters # parameter estimates
    
    loglikelihood = temp$loglik # records log likelihood
          
    z = zhat # update membership matrix
          
    k = k+1 # increment k
    
    delta = abs(llike[k-1]-llike[k])/(1+abs(llike[k]))
    
    it = k-2
    
    if(delta<tol || it>=itmax) break;
  }
  
  error = llike[k]-llike[k-1]
  
  # edit output so it's basically consistent with meVVV() from MCLUST:
  out = list(modelname="VVV", n=n, d=p, G=G, z=zhat, parameters=parameters,
    loglik=loglikelihood, iteration=k-2, error=error)
  
  return(out)
  

}