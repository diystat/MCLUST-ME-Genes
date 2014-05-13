

## within cluster scatter matrix:
wkmat = function(z, data){
 
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)
  
  clustcount = colSums(z) + .Machine$double.xmin # prevent a/0 situation
  
  muhat = matrix(rep(0,p*G),nrow=p)
    for(k in 1:G){
      muhat[,k] = colSums(data*z[,k])/clustcount[k]
    }
  
  
  dim = c(p,p,G)
  W = array(0, dim=dim)
    
  for(k in 1:G){
    for(i in 1:n){
      temp = data[i,] - muhat[,k]
      W[,,k] = W[,,k] + z[i,k]*tcrossprod(temp)
    }
  }
  
  return(W)
}