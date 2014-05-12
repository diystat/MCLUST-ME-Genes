

obj.fun.VVV = function(param,z,data){
  ### Input:
  ### param --- vector of elements in lower triangular matrix of length Gp(p+1)/2 
  ### z --- matrix of membership probabilities. dimension = n*K. z[i,k]=P(obs. i is in cluster k)
  ### data --- matrix of data. dimension = n*d.
  
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
  ###    F(Sigma) = sum_k trace(W_k %*% Sigma^-1) + sum_k n_k*log(det(Sigma_k))
  ### where W = sum_i sum_k z(i,k)(x_i-mu_k)(x_i-mu_k)^T
    
    W = wkmat(z, data) # see function wmat in "W matrix.R"
  
    clustcount = colSums(z) + .Machine$double.xmin # n_k
  
  
  ### From the last step, the objective function is thus:
  
    maxfun = 0
      for(k in 1:G){
        maxfun = maxfun + sum(diag(W[,,k]%*%solve(cov.mat[,,k]))) + clustcount[k]*log(det(cov.mat[,,k]))
      }
  
  
    return(maxfun)
    
}





#----------------------------------------------------------------------------------------#




## debugging
obj.fun.VVV.test = function(){
  
  p = 2
  n = 100
  G = 3
  
  ini.mat = matrix(0, p, p)
  diag(ini.mat) = 1
  ini.par = rep(lowerTriangle(ini.mat, diag=TRUE), G)
  
  obj.fun.VVV(ini.par, z, samp) # error: system singular
  
  
  ## check if L is correct
  m = p*(p+1)/2
  
    L = array(0, dim=c(p, p, G))
      for(k in 1:G){
        lowerTriangle(L[,,k],diag=TRUE) = ini.par[m*k-m+1:m*k]
      }
  
  k = 2
  ini.par[4:6]
  ini.par[(m*k-m+1):(m*k)] # should include parentheses, otherwise incorrect
  
}



