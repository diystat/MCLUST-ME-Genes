

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
    library(gdata)
  
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





