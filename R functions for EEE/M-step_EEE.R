


### M-step function:

MstepEEE = function(z, data){
  library(gdata)
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)
  
  ini.mat = wmat(z, data)/n # use analytic solution as initial value
  # ini.mat = matrix(0, p, p)
  # diag(ini.mat) = 1
  ini.par = lowerTriangle(t(chol(ini.mat)), diag=TRUE)
  
  
  lower.bound = matrix(-Inf, p, p)
  library(gdata)
  diag(lower.bound) = 0.001
  lower = lowerTriangle(lower.bound, diag=TRUE)
  
  
  ## find arg max for objective functin w.r.t. the covariance matrix
  param.est = optim(par=ini.par, obj.fun, z=z, data=data, 
    lower=lower, method="L-BFGS-B")$par
   
  
  
  clustcount = colSums(z) + .Machine$double.xmin # n_k
  phat = clustcount/n # n_k/n
  muhat = matrix(rep(0,p*G),nrow=p) # muhat
  for(k in 1:G){
    muhat[,k] = colSums(data*z[,k])/clustcount[k]
  }
  
  ## convert parameter estimates back to cov matrix:
  L = matrix(0, p, p)
  lowerTriangle(L, diag=TRUE) = param.est
  sigmahat = L %*% t(L)
  
  ## put estimates together into a list:
  parameters = list(muhat=muhat, phat=phat, sigmahat=sigmahat)
  
  return(parameters)
}




#------------------------------------------------------------------------------------------#




### debugging
MstepEEE.test.fun = function(){
  
    library(gdata)
    n = 100;
    p = 3;
    G = 3;

    set.seed(999);
    x = rnorm(n*p);
    dim(x) = c(n, p);
    data = x

    z = matrix(0, n, G);
    z[,1] = 1;
  
  ini.mat = wmat(z, data)/n
  ini.par = lowerTriangle(t(chol(ini.mat)), diag=TRUE)
  
  
  lower.bound = matrix(-Inf, p, p)
  diag(lower.bound) = 0.001
  lower = lowerTriangle(lower.bound, diag=TRUE)
  
  
  ## find arg max for objective functin w.r.t. the covariance matrix
  param.est = optim(par=ini.par, obj.fun, z=z, data=data, 
    lower=lower, method="L-BFGS-B")$par
   
}



















