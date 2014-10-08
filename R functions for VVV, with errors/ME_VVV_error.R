

setwd("Users/wzhang/MS final project/R functions for VVV, with errors")
### manual EM algorithm for VVV with estimation error

ME.VVV.err = function(data, z, err, errstr){
  # Argument "errstr" is user-specified error structure.
  # If errstr="identical", all errors are the same.
  # If errstr="cluster", errors are the same within each cluster.
  # If errstr="none", no constraints on error structure.
  
  source("E-step_VVV_error.R")
  source("M-step_VVV_error.R")
  source("log likelihood_VVV_error.R")
  source("objective function_VVV_error.R")
  source("W_k matrix.R")
  source("inipar.R")
  source("misclassification rate.R")

  
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)

  k = 2 # keeps track of number of iteration

  piconst = n*p*log(2*pi)/2
  
  loglikelihood = NA
  parameters = list()
  zhat = matrix(0,n,G)

  llike = rep(0, 1000) # set convergence criterion
  llike[1] = -10000
  llike[2] = llike[1] + 1e-4
 
  tol = 1e-5
  # while loop for iteration:
  while(llike[k]-llike[k-1] > tol){
    
    print(paste("iteration =",k-1)) # prints number of iterations
    
     # initial values for M-step
    if(errstr=="identical"){
      ini.par = ini.par.iderr(data, z, err)
    } else if(errstr=="cluster"){
      ini.par = ini.par.clust(data, z, err)
    } else if(errstr=="none"){
      ini.par = ini.par.no(data, z)
    }
      
    thetahat = MstepVVV.err(z, data, err, ini.par) # M-step
    
    temp = EstepVVV.err(thetahat, data, err) # E-step
  
    zhat = temp[[1]] # membership estimates
    
    parameters = temp[[2]] # parameter estimates
    
    loglikelihood = temp[[3]] # records log likelihood
          
    z = zhat # update membership matrix
    
    llike[k+1] = loglikelihood # update log likelihood of observed data
    
    print(loglikelihood) # prints current evaulated observed data log-likelihood
      
    k = k+1 # increment k
  }
  
  error = llike[k]-llike[k-1]
  
  # edit output so it's basically consistent with meVVV() from MCLUST:
  out = list(modelname="VVV with est error", n=n, d=p, G=G, z=zhat,
    parameters=parameters, loglik=loglikelihood, iteration=k-2, error=error, likvec=llike[1:k][-(1:2)])
  
  return(out)
  

}