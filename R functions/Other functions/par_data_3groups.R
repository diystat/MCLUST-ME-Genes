
## Simulate 3-component data
## Specifying model parameters used to generate data
sim.par3 = function(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p){
  out = list(N=N,tau1=tau1,tau2=tau2,mu1=mu1,mu2=mu2,mu3=mu3,
    sig1=sig1,sig2=sig2,sig3=sig3,k=k,p=p)
  return(out)
}


## Generate data using given parameters
sim.data3 = function(simpar){  
  N = simpar$N; tau1 = simpar$tau1; tau2 = simpar$tau2
  mu1 = simpar$mu1; mu2 = simpar$mu2; mu3 = simpar$mu3
  sig1 = simpar$sig1; sig2 = simpar$sig2; sig3 = simpar$sig3
  k = simpar$k; p = simpar$p  
  E = k*matrix(c(1,0,0,1),nrow=2)
    
  ## Randomly assign constant measurement error to 100p% of observations:
  index = rbinom(N,1,p)  
  errmat = array(0,dim=c(2,2,N))
  for(i in 1:N){
    errmat[,,i] = index[i]*E
  }
 
  ## Generate sample from mixture distribution:
  U = runif(N)
  z.ini = matrix(0,N,2)
  rand.samples = matrix(0,N,2)
  for(i in 1:N){
    if(U[i]<tau1){
      rand.samples[i,] = MASS::mvrnorm(1,mu1,(sig1+errmat[,,i]))
      z.ini[i,] = c(1,0)
    } else if(U[i]<(tau1+tau2)){
      rand.samples[i,] = MASS::mvrnorm(1,mu2,(sig2+errmat[,,i]))
      z.ini[i,] = c(0,1)
    } else{
      rand.samples[i,] = MASS::mvrnorm(1,mu3,(sig3+errmat[,,i]))
    }
  }

  out = list(data=rand.samples, z.ini=z.ini, err=errmat, index=index, k=k)
  return(out)
  
}






