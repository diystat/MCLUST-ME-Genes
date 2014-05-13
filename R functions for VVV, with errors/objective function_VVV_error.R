

### objective function for VVV case, with estimation error

obj.fun.VVV.err =  function(param,z,data,err){
  ### Input:
  ### param --- vector of elements in lower triangular matrix of length Gp(p+1)/2 
  ### z --- matrix of membership probabilities. dimension = n*G. z[i,k]=P(obs. i is in cluster k)
  ### data --- matrix of data. dimension = n*d
  ### err --- array of estimation error matrices, of dimension (p, p, n), one for each observation
  
    n=nrow(data)
    p=ncol(data)
    G=ncol(z)
    length(param) = G*p*(p+1)/2
    piconst = n*p*log(2*pi)/2
  
  ### Construct covariance matrix:
  ### First convert param into lower triangular matrices, using lowerTriangle() from package GDATA
    # library(gdata)
  
    m = p*(p+1)/2
  
    L = array(0, dim=c(p, p, G))
      for(k in 1:G){
        lowerTriangle(L[,,k],diag=TRUE) = param[(m*k-m+1):(m*k)]
      }
    
  
  ### Then obtain the cov matrices:
    cov.mat = array(0, dim=c(p, p, G))
      for(k in 1:G){
        cov.mat[,,k] = L[,,k] %*% t(L[,,k])
      }
  
    
  ### After rewritting the log-likelihood for complete data, the objective function
  ### becomes:
  ###    F(Sigma) = sum_i sum_k z_ik*t(x_i-muhat_k)*(Sigma_k+Sigma_i)^(-1)*(x_i-muhat_k) + sum_i sum_k z_ik*log(det(Sigma_k+Sigma_i))
  ###
  ### Notice that in this case, MLE for mean is no longer trivial. Instead, the expression for muhat_k is:
  ###    muhat_k = (sum_i z_ik(Sigma_k+Sigma_i)^(-1))^(-1) * sum_i z_ik(Sigma_k+Sigma_i)^(-1)x_i
    
    
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
  
  
  
  ### From the last step, the objective function is thus:
    
    maxfun = 0
      for(i in 1:n){
        for(k in 1:G){
          temp3 = data[i,]-muhat[,k]
          L = chol(inv.sum[,,k,i])
          maxfun = maxfun + z[i,k] * log(det(cov.mat[,,k]+err[,,i])) + z[i,k] * crossprod(L%*%temp3)
        }
      }
  
  
  ## obtain mixing proportion estimate:
    clustcount = colSums(z) + .Machine$double.xmin # n_k
    phat = clustcount/n # n_k/n
  
  
  ## maximizing the likelihood is minimizing its opposite:
    out = (1/2) * maxfun + piconst - sum(colSums(z)*log(phat))
  
    return(out)
    
}



#------------------------------------------------------------------#



objfun.test = function(){
  
  # when error is zero, expect value to be the same as loglik in error-free case.
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

  param = numeric()
    for(k in 1:G){
      param = c(param, lowerTriangle(t(chol(sigma[,,k])), diag=TRUE))
    }
  
  obj.fun.VVV(param,z,data)
  obj.fun.VVV.err(param,z,data,err)
  
  
}

















































