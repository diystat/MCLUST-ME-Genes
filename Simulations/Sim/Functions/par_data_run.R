
## Specifying model parameters used to generate data
sim.par = function(N,tau,mu1,mu2,sig1,sig2,k,p){
  out = list(N=N,tau=tau,mu1=mu1,mu2=mu2,sig1=sig1,sig2=sig2,k=k,p=p)
  return(out)
}


## Generate data using given parameters
sim.data = function(simpar){  
  N = simpar$N; tau = simpar$tau
  mu1 = simpar$mu1; mu2 = simpar$mu2
  sig1 = simpar$sig1; sig2 = simpar$sig2
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
    if(U[i]<tau){
      rand.samples[i,] = MASS::mvrnorm(1,mu1,(sig1+errmat[,,i]))
      z.ini[i,] = c(1,0)
    } else{
      rand.samples[i,] = MASS::mvrnorm(1,mu2,(sig2+errmat[,,i]))
      z.ini[i,] = c(0,1)
    }
  }

  out = list(data=rand.samples, z.ini=z.ini, err=errmat, index=index, k=k)
  return(out)
  
}


## Obtain results from MCME and meVVV
sim.run = function(simdata){
  rand.samples = simdata$data
  z.ini = simdata$z.ini
  errmat = simdata$err
  index = simdata$index
  k = simdata$k
  
  ## Run MCME:
  res.mcme = mcmeVVV(rand.samples, z.ini, errmat)
  
  ## Run mevvv:
  res.mevvv = mclust::meVVV(rand.samples,z.ini)
  
  out = list(res.mcme=res.mcme, res.mevvv=res.mevvv, z.ini=z.ini, rand.samples=rand.samples, errmat=errmat, index=index, k=k)
  return(out)
}

