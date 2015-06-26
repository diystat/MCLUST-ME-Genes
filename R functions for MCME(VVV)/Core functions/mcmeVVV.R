
### manual EM algorithm for VVV with estimation error

mcmeVVV = function(data, z, err, d=1, itmax=Inf, lb=1e-3){
  # "d" denotes the diagonal element for initial value of input
  # decomposed matrix, default = 1
  
  # Argument "lb" sets lower bound for diagonal elements of decomposed cov
  # matrices. Default set at 1e-3=0.001. Can be set to a larger number to
  # avoid singularity issues.
  
  #source("E-step_VVV_error.R")
  #source("M-step_VVV_error.R")
  #source("log likelihood_VVV_error.R")
  #source("objective function_VVV_error.R")
  #source("W_k matrix.R")
  #source("inipar.R")
  #source("inv_sum.R")

  ## debug(obj.fun.VVV.err);
  
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
  
  # while loop for iteration:
  repeat{
    
    #print(paste("iteration =",k-1)) # prints number of iterations
    
    # initial values for M-step
    ini.par = ini.par.no(data, z, d)
    
    thetahat = MstepVVV.err(z, data, err, ini.par, lb) # M-step
    
    temp = EstepVVV.err(thetahat, data, err) # E-step
  
    zhat = temp[[1]] # membership estimates
    
    parameters = temp[[2]] # parameter estimates
    
    loglikelihood = temp[[3]] # records log likelihood
          
    z = zhat # update membership matrix
    
    member = rbind(member,zhat) # store membership estimates
    
    center = rbind(center,parameters$muhat) # store mean estimates
    
    llike[k+1] = loglikelihood # update log likelihood of observed data
              
    k = k+1 # increment k
    
    delta = abs(llike[k-1]-llike[k])/(1+abs(llike[k]))
    
    it = k-2
    
    if(delta<tol || it>=itmax) break;
  }
  
  error = llike[k]-llike[k-1]
  
  uncertainty = numeric() # records classification uncertainty of each obs.
  for(i in 1:n){
    rowmax = max(zhat[i,])
    uncertainty[i] = 1-rowmax
  }
  
  # edit output so it's basically consistent with meVVV() from MCLUST:
  out = list(modelname="VVV with est error", n=n, d=p, G=G, z=zhat, parameters=parameters, uncertainty=uncertainty,
    loglik=loglikelihood, iteration=k-2, error=error, member=member, center=center, likvec=llike[1:k][-(1:2)])
  
  return(out)
  

}