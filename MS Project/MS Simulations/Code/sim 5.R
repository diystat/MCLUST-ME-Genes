
### Simulation 5: Effect of sample size on running times
library(MASS)
library(gdata)
library(mclust)


p = 2

mu1 = c(1,1)
mu2 = c(5,-5)
mu3 = c(-7,7)

sigma1 = matrix(c(5,-2,-2,5),nrow=2)
sigma2 = matrix(c(7,3,3,7),nrow=2)
sigma3 = matrix(c(9,-4,-4,9),nrow=2)

s = c(1,2,5,10,15,20)

result = list()

for(i in 1:length(s)){
  nvec = c(3,3,4) * s[i]
  n = sum(nvec)

  # set errors:
  err = array(0, dim=c(p,p,n))  
  for(i in 1:n){
    err[,,i] = matrix(0,p,p)
    diag(err[,,i]) = 10
  }
  
  set.seed(0)
  s1 = mvrnorm(nvec[1], mu1, (sigma1+err[,,1]))
  s2 = mvrnorm(nvec[2], mu2, (sigma2+err[,,1]))
  s3 = mvrnorm(nvec[3], mu3, (sigma3+err[,,1]))

  # true membership matrix:
  temp = c(rep(c(1,0,0),nvec[1]),rep(c(0,1,0),nvec[2]),rep(c(0,0,1),nvec[3]))
  z.ini = z.true = matrix(temp, nrow=n, byrow=TRUE)

  samp = rbind(s1,s2,s3)

  # record running times:
  tmp = proc.time()
  my.result.run1 = ME.VVV.err(samp, z.ini, err)
  my.time.run1 = proc.time() - tmp

  tmp = proc.time()
  my.result.run2 = ME.VVV.err(samp, z.ini, err)
  my.time.run2 = proc.time() - tmp

  tmp = proc.time()
  my.result.run3 = ME.VVV.err(samp, z.ini, err)
  my.time.run3 = proc.time() - tmp

  my.time = list(my.time.run1,my.time.run2,my.time.run3)


  tmp = proc.time()
  mc.result1 = meVVV(samp,z.ini)
  mc.time1 = proc.time() - tmp

  tmp = proc.time()
  mc.result2 = meVVV(samp,z.ini)
  mc.time2 = proc.time() - tmp

  tmp = proc.time()
  mc.result3 = meVVV(samp,z.ini)
  mc.time3 = proc.time() - tmp

  mc.time = list(mc.time1,mc.time2,mc.time3)


  result[[i]] = list(my.time,mc.time)
  
}

save(result,file="sim5_result.RData")
















