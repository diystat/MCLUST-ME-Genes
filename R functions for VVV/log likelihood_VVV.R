


## log likelihood of mvn mixture, VVV case, with estimation error:
log.lik.vvv = function(param, data){
    
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
      temp1 = sigmahat[,k]
      L = chol(temp)
      temp2 = data[i,]-muhat[,k]
      temp[i] = temp[i] + tauhat[k] * (det(sigmahat[,,k]))^(-1/2) * exp((-1/2) * 
          crossprod(L%*%temp2))      
    }
    lik = lik + log(temp[i])
  }
  
  out = lik + piconst
  
  return(out)
}