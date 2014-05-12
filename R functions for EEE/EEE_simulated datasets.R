
library(MASS)
### simulate from normal mixture:

  mu1 = c(2,2)
  mu2 = c(1,-2)
  mu3 = c(-2,1)

  sigma = matrix(c(2,1,1,4),nrow=2) # same cov matrix for all components

  p = c(0.1,0.3,0.6)

# The number of samples from the mixture distribution
  N = 100               

# Sample N random uniforms U
  U = runif(N)

#Variable to store the samples from the mixture distribution                                             
  rand.samples = matrix(rep(NA,2*N),nrow=N)

# Sampling from the mixture
  for(i in 1:N){
      if(U[i]<p[1]){
          rand.samples[i,] = mvrnorm(1,mu1,sigma)
      }else if(U[i]<(1-p[1])){
          rand.samples[i,] = mvrnorm(1,mu2,sigma)
      }else{
          rand.samples[i,] = mvrnorm(1,mu3,sigma)
      }
  }


# initial value for membership matrix:
  temp = c(rep(0,200),rep(1,100))
  z0 = matrix(temp, nrow=100, byrow=T)



#---------------------------------------------------------------------------------#




## simulate three groups of multivariate normal sample:

  library(MASS)

  mu1 = c(9,9)
  mu2 = c(4,-9)
  mu3 = c(-9,4)

  nvec = c(30,30,40)
  n = sum(nvec)

  sigma = matrix(c(2,1,1,2),nrow=2) # same cov matrix for all groups

  s1 = mvrnorm(nvec[1], mu1, sigma)
  s2 = mvrnorm(nvec[2], mu2, sigma)
  s3 = mvrnorm(nvec[3], mu3, sigma)

  # membership matrix:
  temp = c(rep(c(1,0,0),30),rep(c(0,1,0),30),rep(c(0,0,1),40))
  z = z0 = matrix(temp, nrow=n, byrow=TRUE)

  samp = rbind(s1,s2,s3)
  plot(samp) # very few overlaps between groups



#-------------------------------------------------------------------------------------#


## run EM on the second dataset:
  ME.EEE(samp, z)

## compare with meEEE() from MCLUST
  library(mclust)
  meEEE(samp, z) # same results here



