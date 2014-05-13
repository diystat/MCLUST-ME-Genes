

cc = function(param){
  
  nvec = c(3,3,4) * 5
  n = sum(nvec)
  p = 2
  
  mu1 = c(9,9)
  mu2 = c(25,-25)
  mu3 = c(-49,49)
  
  err = array(0, dim=c(p,p,n))  
  for(i in 1:n){
    err[,,i] = matrix(0,p,p)
    diag(err[,,i]) = 0.01
  }
  
  sigma1 = matrix(c(4,1,1,9),nrow=2)
  sigma2 = matrix(c(16,2,2,25),nrow=2)
  sigma3 = matrix(c(36,3,3,49),nrow=2)
  
  s1 = mvrnorm(nvec[1], mu1, (sigma1+err[,,1]))
  s2 = mvrnorm(nvec[2], mu2, (sigma2+err[,,1]))
  s3 = mvrnorm(nvec[3], mu3, (sigma3+err[,,1]))
  
  # membership matrix:
  temp = c(rep(c(1,0,0),nvec[1]),rep(c(0,1,0),nvec[2]),rep(c(0,0,1),nvec[3]))
  z.true = matrix(temp, nrow=n, byrow=TRUE)
  
  data = rbind(s1,s2,s3)
  
  z = z.true
  
  G=ncol(z)

  piconst = n*p*log(2*pi)/2
  m = p*(p+1)/2
  
  L = array(0, dim=c(p, p, G))
  for(k in 1:G){
    lowerTriangle(L[,,k],diag=TRUE) = param[(m*k-m+1):(m*k)]
  }
   
  cov.mat = array(0, dim=c(p, p, G))
  for(k in 1:G){
    cov.mat[,,k] = L[,,k] %*% t(L[,,k])
  }
    
  muhat = matrix(0, p, G)
  inv.sum = array(0, dim=c(p,p,G,n))
  temp1 = 0
  temp2 = rep(0, p)
  for(k in 1:G){
    for(i in 1:n){
      inv.sum[,,k,i] = solve(cov.mat[,,k]+err[,,i])
      temp1 = temp1 + z[i,k]*inv.sum[,,k,i]
      temp2 = temp2 + z[i,k]*inv.sum[,,k,i]%*%data[i,]
    }
    muhat[,k] = solve(temp1) %*% temp2
  }
    
  maxfun = 0
  for(i in 1:n){
    for(k in 1:G){
      temp3 = data[i,]-muhat[,k]
      L = chol(inv.sum[,,k,i])
      maxfun = maxfun + z[i,k] * log(det(cov.mat[,,k]+err[,,i])) + z[i,k] * crossprod(L%*%temp3)
    }
  }
  
  clustcount = colSums(z) + .Machine$double.xmin # n_k
  phat = clustcount/n # n_k/n
  
  out = (1/2) * maxfun + piconst - sum(colSums(z)*log(phat))
  
  
  return(out)
}






  grid = fmin.gridsearch(fun=cc, x0=ini.par, xmin=lower, xmax=upper, npts=30)











































