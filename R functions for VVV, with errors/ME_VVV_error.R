


### manual EM algorithm for VVV with estimation error

ME.VVV.err = function(data, z, err, errstr){
  # Argument "errstr" is user-specified error structure.
  # If errstr="identical", all errors are the same.
  # If errstr="cluster", errors are the same within each cluster.
  # If errstr="no", no constraints on error structure.
  
  #source("~/R functions for VVV, with errors/E-step_VVV_error.R")
  #source("~/R functions for VVV, with errors/M-step_VVV_error.R")
  #source("~/R functions for VVV, with errors/log likelihood_VVV_error.R")
  #source("~/R functions for VVV, with errors/objective function_VVV_error.R")
  #source("~/R functions for VVV, with errors/W_k matrix.R")
  #source("~/R functions for VVV, with errors/inipar.R")
  #source("~/R functions for VVV, with errors/misclassification rate.R")

  
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)

  k = 2 # keeps track of number of iteration

  piconst = n*p*log(2*pi)/2
  
  loglikelihood = NA
  parameters = list()
  zhat = matrix(0,n,G)

  llike = rep(0, 10000) # set convergence criterion
  llike[2] = 1e-4
 
  tol = 1e-5
  # while loop for iteration:
  while(abs(llike[k]-llike[k-1]) > tol){
    
    print(paste("iteration =",k-1)) # prints number of iterations
    
     # initial values for M-step
    if(errstr=="identical"){
      ini.par = ini.par.iderr(data, z, err)
    } else if(errstr=="cluster"){
      ini.par = ini.par.clust(data, z, err)
    } else if(errstr=="no"){
      ini.par = ini.par.no(data, z)
    }
      
    thetahat = MstepVVV.err(z, data, err, ini.par) # M-step
    
    temp = EstepVVV.err(thetahat, data, err) # E-step
  
    zhat = temp[[1]] # membership estimates
    
    parameters = temp[[2]] # parameter estimates
    
    loglikelihood = temp[[3]] # records log likelihood
          
    z = zhat # update membership matrix
    
    llike[k+1] = loglikelihood # update log likelihood of complete data
    
    print(loglikelihood) # prints current evaulated complete data log-likelihood
      
    k = k+1 # increment k
  }
  
  # edit output so it's basically consistent with meVVV() from MCLUST:
  out = list(modelname="VVV with est error", n=n, d=p, G=G, z=zhat,
    parameters=parameters, loglik=loglikelihood, iteration=k-2)
  
  return(out)
  

}