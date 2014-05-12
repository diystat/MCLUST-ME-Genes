

## log likelihood of mvn mixture, EEE case:
log.lik.eee = function(param, data){
    
  muhat = param[[1]]
  tauhat = param[[2]]
  sigmahat = param[[3]]
  
  n = nrow(data)
  p = ncol(data)
  G = ncol(param[[1]])
  
  piconst = -n*p*log(2*pi)/2
  
  temp = rep(0,n)
  lik = 0
  for(i in 1:n){
    for(k in 1:G){
      temp[i] = temp[i] + tauhat[k] * exp((-1/2) * 
          t(data[i,]-muhat[,k]) %*% solve(sigmahat) %*% (data[i,]-muhat[,k]))      
    }
    lik = lik + log(temp[i])
  }
  
  out = lik - (n/2)*log(det(sigmahat)) + piconst
  
  return(out)
}


