

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
  
  ### Construct covariance matrix:
  ### First convert param into lower triangular matrices, using lowerTriangle() from package GDATA
    library(gdata)
  
    m = p*(p+1)/2
  
    L = array(0, dim=c(p, p, G))
      for(k in 1:G){
        lowerTriangle(L[,,k],diag=TRUE) = param[(m*k-m+1):(m*k)]
      }
    
    # print(L)
  
  ### Then obtain the cov matrices:
    cov.mat = array(0, dim=c(p, p, G))
      for(k in 1:G){
        cov.mat[,,k] = L[,,k] %*% t(L[,,k])
      }
  
    # print(cov.mat)
  
  
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
  
  
    return(maxfun)
    
}



#--------------------------------------------------------------------------#




## debugging
obj.fun.VVV.err.test = function(){
  
  data = samp
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)
  
# error array:
  err = array(dim=c(p,p,n))  
    for(i in 1:n){
      err[,,i] = matrix(0,p,p)
      diag(err[,,i]) = 1
    }
  
  
   ini.mat = matrix(0, p, p)
   diag(ini.mat) = 1
  
  
  ini.par = numeric()
    for(k in 1:G){
      ini.par = c(ini.par, lowerTriangle(t(chol(ini.mat[,,k])), diag=TRUE)) # extract initial lower triangular elements
    }
  
  
  
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
  system.time(optim(par=ini.par, obj.fun.VVV.err, z=z, data=data, err=err, 
    lower=lower, method="L-BFGS-B")) # 33s, converged
}






