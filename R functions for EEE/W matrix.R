
wmat = function(z, data){
 
  n = nrow(data)
  p = ncol(data)
  G = ncol(z)
  
  clustcount = colSums(z) + .Machine$double.xmin # prevent a/0 situation
  
  muhat = matrix(rep(0,p*G),nrow=p)
    for(k in 1:G){
      muhat[,k] = colSums(data*z[,k])/clustcount[k]
    }
  
  W = matrix(0, p, p)
    
  for(i in 1:n){
    for(k in 1:G){
      W = W + z[i,k]*(data[i,]-muhat[,k])%*%t(data[i,]-muhat[,k])
    }
  }
  
  return(W)
}



#-----------------------------------------------------------------------------------------#


## debugging
wmat.test = function(){
    n = 100;
    p = 3;
    G = 3;

    set.seed(999);
    x = rnorm(n*p);
    dim(x) = c(n, p);

    z = matrix(0, n, G);
    z[,1] = 1;

    ## Test:
    W = wmat(z,x);W
    det(W)
   
   ## The original function returns NaN for muhat, we found out that
   ## it results from zero counts in some clusters. So we added a small
   ## number to each cluster count:
   clustcount = colSums(z) + .Machine$double.xmin
  
   muhat = matrix(rep(0,p*G),nrow=p)
     for(k in 1:G){
       muhat[,k] = colSums(x*z[,k])/clustcount[k]
     }
  
   W1 = matrix(0, p, p)    
    for(i in 1:n){
      for(k in 1:G){
        W1 = W1 + z[i,k]*(x[i,]-muhat[,k])%*%t(x[i,]-muhat[,k])
      }
    }
  
   W1 # the NaN problem is solved
  
  W == W1
  
  ## and we made the same modifications to the same part for other
  ## functions
}

