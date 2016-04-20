
### Function that outputs inverses of sum of covariances for each observation
inv_sum = function(z, data, err, sigmahat){
  
  p = ncol(data)
  n = nrow(data)
  G = ncol(z)
  
  ue = unique(err,MARGIN=3)
  n_unq = dim(ue)[3]
  unique_inv = array(dim=c(p,p,G,n_unq))
  for(k in 1:G){
    for(i in 1:n_unq){
      unique_inv[,,k,i] = chol2inv(chol(sigmahat[,,k]+ue[,,i]))
    }
  }
  
  dim(ue) = c(p*p,n_unq)
  dim(err) = c(p*p,n)
  vec = numeric(n)
  for(j in 1:n_unq){
    vec[colSums(err==ue[,j])==p*p]=j
  }
  
  out = array(dim=c(p,p,G,n))
  for(k in 1:G){
    for(i in 1:n){
      out[,,k,i] = unique_inv[,,k,(vec[i])]
    }
  }
  
  
  return(out)  
}

