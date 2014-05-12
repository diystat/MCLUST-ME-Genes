

### M-step function for VVV with estimation error

MstepVVV.iderr = function(z, data, err){
  library(gdata)
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)
  
  
  ## obtain mixing proportion and estimates:
  clustcount = colSums(z) + .Machine$double.xmin # n_k
  phat = clustcount/n # n_k/n
  muhat = matrix(rep(0,p*G),nrow=p) # muhat
  for(k in 1:G){
    muhat[,k] = colSums(data*z[,k])/clustcount[k]
  }
  
    
  ## acquire covariance matrix esimates, using analytic result from Celeux and Govaert (1995):
  wk = wkmat(z, data)
  dim = c(p, p, G)
  sigmahat = array(0,dim=dim)
    for(k in 1:G){
      sigmahat[,,k] = wk[,,k]/clustcount[k] - err
    }
  
  
  
  ## put estimates together into a list:
  parameters = list(muhat=muhat, phat=phat, sigma=sigmahat)
  
  return(parameters)
}




