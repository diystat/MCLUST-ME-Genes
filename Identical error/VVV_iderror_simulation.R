

library(MASS)
library(mclust)

## simulate well-separated clusters, with correct membership as initial membership matrix:  
sim.goodini = function(s){
  nvec = c(3,3,4) * s
  n = sum(nvec)
  p = 2
  
  mu1 = c(9,9)
  mu2 = c(25,-25)
  mu3 = c(-49,49)
  
  err = matrix(0,p,p) 
  diag(err) = 1
  
  sigma1 = matrix(c(4,0,0,9),nrow=2)
  sigma2 = matrix(c(16,0,0,25),nrow=2)
  sigma3 = matrix(c(36,0,0,49),nrow=2)
  
  s1 = mvrnorm(nvec[1], mu1, (sigma1+err))
  s2 = mvrnorm(nvec[2], mu2, (sigma2+err))
  s3 = mvrnorm(nvec[3], mu3, (sigma3+err))
  
  # membership matrix:
  temp = c(rep(c(1,0,0),nvec[1]),rep(c(0,1,0),nvec[2]),rep(c(0,0,1),nvec[3]))
  z.true = matrix(temp, nrow=n, byrow=TRUE)
  
  samp = rbind(s1,s2,s3)
  
  z.ini = z.true
  
  ptm1 = proc.time()  
  my.result = ME.VVV.iderr(samp, z.ini, err)
  my.class = my.result$z
  my.time = proc.time() - ptm1
  
  my.mcr = MCR(my.class, z.true)
  
  ptm2 = proc.time()
  mc.result = meVVV(samp, z.ini)
  mc.class = mc.result$z
  mc.time = proc.time() - ptm2
  
  mc.mcr = MCR(mc.class, z.true)  
  
  out = list(my.time=my.time, my.mcr=my.mcr, my.result=my.result, mc.time=mc.time, mc.mcr=mc.mcr, mc.result=mc.result)
  
  return(out)
}
  
  
  

## simulate the same clusters as above, but with incorrect initial membership matrix:  
sim.badini = function(s){
  
  nvec = c(3,3,4) * s
  n = sum(nvec)
  p = 2
  
  mu1 = c(9,9)
  mu2 = c(25,-25)
  mu3 = c(-49,49)
  
  err = matrix(0,p,p) 
  diag(err) = 1
  
  sigma1 = matrix(c(4,0,0,9),nrow=2)
  sigma2 = matrix(c(16,0,0,25),nrow=2)
  sigma3 = matrix(c(36,0,0,49),nrow=2)
  
  s1 = mvrnorm(nvec[1], mu1, (sigma1+err))
  s2 = mvrnorm(nvec[2], mu2, (sigma2+err))
  s3 = mvrnorm(nvec[3], mu3, (sigma3+err))
  
  # membership matrix:
  temp = c(rep(c(1,0,0),nvec[1]),rep(c(0,1,0),nvec[2]),rep(c(0,0,1),nvec[3]))
  z.true = matrix(temp, nrow=n, byrow=TRUE)
  
  samp = rbind(s1,s2,s3)
  
  temp2 = c(rep(c(0,1,0),nvec[1]),rep(c(0,0,1),nvec[2]),rep(c(1,0,0),nvec[3]))
  z.ini = matrix(temp2, nrow=n, byrow=TRUE)
  
  ptm1 = proc.time()
  my.result = ME.VVV.iderr(samp, z.ini, err)
  my.class = my.result$z
  my.time = proc.time() - ptm1
  
  my.mcr = MCR(my.class, z.true)
  
  ptm2 = proc.time()
  mc.result = meVVV(samp, z.ini)
  mc.class = mc.result$z
  mc.time = proc.time() - ptm2
  
  mc.mcr = MCR(mc.class, z.true)
  
  out = list(my.time=my.time, my.mcr=my.mcr, my.result=my.result, mc.time=mc.time, mc.mcr=mc.mcr, mc.result=mc.result)
  
  return(out)
}
  
  
 
  
  