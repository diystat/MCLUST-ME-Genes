N = 100
p = 0.5
tau = 0.5
  mu1 = c(0,0)
  mu2 = c(0,0)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)

  ## Randomly assign constant measurement error to half of observations:
  index = rbinom(N,1,p)
  errmat = array(0,dim=c(2,2,N))
  E = matrix(c(16,0,0,16),nrow=2)
    for(i in 1:N){
      errmat[,,i] = index[i]*E
    }

  ## Generate sample from mixture distribution:
  U = runif(N)
  data = matrix(0,N,2)
  for(i in 1:N){
    if(U[i]<tau){
      data[i,] = mvrnorm(1,mu1,(sig1+errmat[,,i]))
    } else{
      data[i,] = mvrnorm(1,mu2,(sig2+errmat[,,i]))
    }
  }

## True membership:
  z.ini = matrix(0,N,2)
  for(i in 1:N){
    if(U[i]<tau){
      z.ini[i,] = c(1,0)
    } else{
      z.ini[i,] = c(0,1)
    }
  }

source("oneit.R")
l <- lineprof(oneit(data,z.ini,errmat))
l
shine(l)

focus(l,f="MstepVVV.err")

