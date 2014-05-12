

### M-step function:

MstepVVV = function(z, data){
  library(gdata)
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)
  
  
  
  clustcount = colSums(z) + .Machine$double.xmin # n_k
  phat = clustcount/n # n_k/n
  muhat = matrix(rep(0,p*G),nrow=p) # muhat
  for(k in 1:G){
    muhat[,k] = colSums(data*z[,k])/clustcount[k]
  }
  
  
  ## acquire covariance matrix esimates, using analytic result from Celeux and Govaert (1995):
  wk = wkmat(z, data)
  dim = c(p, p, G)
  sigma = array(0,dim=dim)
    for(k in 1:G){
      sigma[,,k] = wk[,,k]/clustcount[k]
    }
  
  
  
  ini.mat = sigma # use analytic solution as initial value
  # ini.mat = matrix(0, p, p)
  # diag(ini.mat) = 1
  
  
  #ini.par = numeric()
  #  for(k in 1:G){
  #    ini.par = c(ini.par, lowerTriangle(t(chol(ini.mat[,,k])), diag=TRUE)) # extract initial lower triangular elements
  #  }
  
  
  
  ini.mat = matrix(0, p, p)
  diag(ini.mat) = 1
  ini.par = rep(lowerTriangle(ini.mat, diag=TRUE), G)
  
  
  
  library(gdata)
  lower.bound = array(-Inf, dim=c(p, p, G))
  lower = numeric()
    for(k in 1:G){
      diag(lower.bound[,,k]) = 0.001
      lower = c(lower, lowerTriangle(lower.bound[,,k], diag=TRUE)) # set lower bound for lower triangle elements
    }
  
  
  
  ## find arg max for objective functin w.r.t. the covariance matrix
  param.est = optim(par=ini.par, obj.fun.VVV, z=z, data=data, 
    lower=lower, method="L-BFGS-B")$par
  
  
  ## transform back to covariance matrices:
  m = p*(p+1)/2
  sigmahat = array(0, dim=c(p, p, G))
    for(k in 1:G){     
      lowerTriangle(sigmahat[,,k], diag=TRUE) = param.est[(m*k-m+1):(m*k)]
      sigmahat[,,k] = sigmahat[,,k] %*% t(sigmahat[,,k])
    }
  
  
  ## put estimates together into a list:
  parameters = list(muhat=muhat, phat=phat, sigma=sigmahat)
  
  return(parameters)
}

