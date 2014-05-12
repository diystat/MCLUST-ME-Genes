

## simulate three groups of multivariate normal sample, each with different cov matrix:

  library(MASS)

  mu1 = c(9,9)
  mu2 = c(4,-9)
  mu3 = c(-9,4)

  nvec = c(30,30,40)
  n = sum(nvec)

  sigma1 = matrix(c(2,1,1,2),nrow=2)
  sigma2 = matrix(c(4,-1,-1,3),nrow=2)
  sigma3 = matrix(c(5,3,3,5),nrow=2)

  s1 = mvrnorm(nvec[1], mu1, sigma1)
  s2 = mvrnorm(nvec[2], mu2, sigma2)
  s3 = mvrnorm(nvec[3], mu3, sigma3)

  # membership matrix:
  temp = c(rep(c(1,0,0),30),rep(c(0,1,0),30),rep(c(0,0,1),40))
  z = matrix(temp, nrow=n, byrow=TRUE)

  samp = rbind(s1,s2,s3)
  plot(samp) # very few overlaps between groups



#-------------------------------------------------------------------------------------#


## run EM on the second dataset:
  ME.VVV(samp, z)

## compare with meEEE() from MCLUST
  library(mclust)
  meVVV(samp, z) # same results here
