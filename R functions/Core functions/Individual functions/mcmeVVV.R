
### manual EM algorithm for VVV with estimation error

mcmeVVV = function(data, z, err, d=1, itmax=Inf, lb=1e-3){
  # "d" denotes the diagonal element for initial value of input
  # decomposed matrix, default = 1
  
  # Argument "lb" sets lower bound for diagonal elements of decomposed cov
  # matrices. Default set at 1e-3=0.001. Can be set to a larger number to
  # avoid singularity issues.
  
  # Argument "itmax" sets an upper bound of the number of iterations
  # of the EM algorithm.
  
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

    ini.par = ini.par.no(data, z, d) # initial values for M-step
    
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
    print(paste("iteration=",it,sep=""))
    
    if(delta<tol || it>=itmax) break;
  }
  
  error = llike[k]-llike[k-1] # log-lik difference of last two iterations
  
  uncertainty = numeric() # records classification uncertainty of each obs.
  for(i in 1:n){
    rowmax = max(zhat[i,])
    uncertainty[i] = 1-rowmax
  }
  
  nu = G*p*(p+1)/2 + G*p + (G-1)
  bic.me = 2*loglikelihood - nu*log(n) # Calculates BIC for the model
  
  # edit output so it's basically consistent with meVVV() from MCLUST:
  out = list(modelname="VVV with est error", BIC=bic.me, n=n, d=p, G=G, z=zhat, parameters=parameters, uncertainty=uncertainty,
    loglik=loglikelihood, iteration=k-2, error=error, member=member, center=center, likvec=llike[1:k][-(1:2)])
  
  return(out)
  

}