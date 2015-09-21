
### M-step function for VVV with estimation error
MstepVVV.err = function(z, data, err, ini.par, lb){
  # argument "ini.par" defines a user-specified initial value for L-BFGS-B algorithm
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)
  
  
  ## obtain mixing proportion estimate:
  clustcount = colSums(z) + .Machine$double.xmin # n_k
  phat = clustcount/n # n_k/n
  
  
  # set lower bounds for parameters:
  #library(gdata)
  lower.bound = array(-Inf, dim=c(p, p, G))
  lower = numeric()
    for(k in 1:G){
      diag(lower.bound[,,k]) = lb
      lower = c(lower, gdata::lowerTriangle(lower.bound[,,k], diag=TRUE)) # set lower bound for lower triangle elements
    }
  
    
  ## find arg max for objective function w.r.t. the covariance matrix
  est.res = optim(par=ini.par, obj.fun.VVV.err, gr=NULL, z=z, data=data, err=err,
     lower=lower, method="L-BFGS-B")
  
   param.est = est.res$par
  
  
  
  ## transform back to covariance matrices:
  m = p*(p+1)/2
  sigmahat = array(0, dim=c(p, p, G))
    for(k in 1:G){     
      gdata::lowerTriangle(sigmahat[,,k], diag=TRUE) = param.est[(m*k-m+1):(m*k)]
      sigmahat[,,k] = tcrossprod(sigmahat[,,k])
    }
  
  
   
  ## obtain mean estimate:
  inv.sum = inv_sum(z,data,err,sigmahat)
  #print(inv.sum)
  muhat = matrix(0, p, G)
  temp1 = temp2 = 0
      for(k in 1:G){
        for(i in 1:n){         
          tem = inv.sum[,,k,i]
          temp1 = temp1 + z[i,k]*tem
          temp2 = temp2 + z[i,k]*tem%*%data[i,]
        }
        muhat[,k] = chol2inv(chol(temp1))%*%temp2
        temp1 = temp2 = 0
      }
  
    
  ## put estimates together into a list:
  parameters = list(muhat=muhat, phat=phat, sigma=sigmahat, inv.sum=inv.sum)
  
  return(parameters)
}



#------------------------------------------------------------------------------------#



mstep.test = function(){
  
  
  ### Case 1: set error to zero. see if results are same as MstepVVV
  library(MASS) 
  library(gdata)
  set.seed(0)  
  mu1 = c(9,9)
  mu2 = c(4,-9)
  mu3 = c(-9,4)
  
  nvec = c(30,30,40)
  n = sum(nvec)
  p = 2
  G = 3
  
  sigma1 = matrix(c(2,1,1,2),nrow=2)
  sigma2 = matrix(c(4,-1,-1,3),nrow=2)
  sigma3 = matrix(c(5,3,3,5),nrow=2)
  
  sigma = array(0,dim=c(p,p,G))
  sigma[,,1] = sigma1
  sigma[,,2] = sigma2
  sigma[,,3] = sigma3
  
  s1 = mvrnorm(nvec[1], mu1, sigma1)
  s2 = mvrnorm(nvec[2], mu2, sigma2)
  s3 = mvrnorm(nvec[3], mu3, sigma3)
  
  # membership matrix:
  temp = c(rep(c(1,0,0),30),rep(c(0,1,0),30),rep(c(0,0,1),40))
  z = matrix(temp, nrow=n, byrow=TRUE)
  
  data = rbind(s1,s2,s3)
  
  err = array(0, dim=c(p,p,n))
  
  inipar = ini.par.no(data,z)
  MstepVVV.err(z, data, err, ini.par=inipar)
  mstepVVV(data,z)
  # results are the same
  
  
  
  ### Case 2: when errors are identical, cov matrix estimates should be equal to W_k/n_k - error
  iderr = array(0, dim=c(p,p,n))
  for(i in 1:n){
    iderr[,,i] = matrix(0,p,p)
    diag(iderr[,,i]) = 2
  }
  
  sigma1.e = sigma1 + iderr[,,1]
  sigma2.e = sigma2 + iderr[,,1]
  sigma3.e = sigma3 + iderr[,,1]
  
  s1.e = mvrnorm(nvec[1], mu1, sigma1.e)
  s2.e = mvrnorm(nvec[2], mu2, sigma2.e)
  s3.e = mvrnorm(nvec[3], mu3, sigma3.e)
  
  data.e = rbind(s1.e,s2.e,s3.e)
  
  nk = colSums(z)
  
  ww = wkmat(z, data.e)
  ww[,,1] = ww[,,1]/nk[1] - iderr[,,1]
  ww[,,2] = ww[,,2]/nk[2] - iderr[,,1]
  ww[,,3] = ww[,,3]/nk[3] - iderr[,,1]
  
  
  MstepVVV.err(z, data.e, iderr, ini.par=inipar)
  ww
  # results are the same
  
  
  
  
  ### Case 3: see if initial values affect cov estimates
  ## obtain mixing proportion estimate:
  clustcount = colSums(z.ini) + .Machine$double.xmin # n_k
  #clustcount = colSums(z)
  phat = clustcount/n # n_k/n
  
    
  ## acquire covariance matrix esimates, using analytic result from Celeux and Govaert (1995):
   wk = wkmat(z.ini, samp)
   dim = c(p, p, G)
   sigma = array(0,dim=dim)
     for(k in 1:G){
       sigma[,,k] = wk[,,k]/clustcount[k] - err[,,1]
     }
  
  
   ini.mat = sigma # use analytic solution as initial value
  
   ini.par = numeric()
     for(k in 1:G){
       ini.par = c(ini.par, lowerTriangle(t(chol(ini.mat[,,k])), diag=TRUE)) # extract initial lower triangular elements
     }

  # set lower bounds for parameters:
  #library(gdata)
  lower.bound = array(-Inf, dim=c(p, p, G))
  lower = numeric()
    for(k in 1:G){
      diag(lower.bound[,,k]) = 1e-8
      lower = c(lower, lowerTriangle(lower.bound[,,k], diag=TRUE)) # set lower bound for lower triangle elements
    }
    
  ## find arg max for objective functin w.r.t. the covariance matrix
   param.est = optim(par=ini.par, obj.fun.VVV.err, z=z.ini, data=samp, err=err,
     lower=lower, method="L-BFGS-B")$par
   param.est
  ## transform back to covariance matrices:
  m = p*(p+1)/2
  sigmahat = array(0, dim=c(p, p, G))
    for(k in 1:G){     
      lowerTriangle(sigmahat[,,k], diag=TRUE) = param.est[(m*k-m+1):(m*k)]
      sigmahat[,,k] = tcrossprod(sigmahat[,,k])
    }
  sigmahat
  
  
}




