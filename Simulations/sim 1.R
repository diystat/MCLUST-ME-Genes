
### Simulation 1: zero errors
library(gdata)
library(mclust)
library(MASS)

# set sample size
nvec = c(3,3,4) * 10
n = sum(nvec)
p = 2

mu1 = c(1,1)
mu2 = c(10,-10)
mu3 = c(-25,25)

# set errors to be zero
err = array(0, dim=c(p,p,n))  
for(i in 1:n){
  err[,,i] = matrix(0,p,p)
  diag(err[,,i]) = 0
}

sigma1 = matrix(c(3,0,0,3),nrow=2)
sigma2 = matrix(c(5,0,0,5),nrow=2)
sigma3 = matrix(c(7,0,0,7),nrow=2)

set.seed(0)
s1 = mvrnorm(nvec[1], mu1, (sigma1+err[,,1]))
s2 = mvrnorm(nvec[2], mu2, (sigma2+err[,,1]))
s3 = mvrnorm(nvec[3], mu3, (sigma3+err[,,1]))

# true membership matrix:
temp = c(rep(c(1,0,0),nvec[1]),rep(c(0,1,0),nvec[2]),rep(c(0,0,1),nvec[3]))
z.true = matrix(temp, nrow=n, byrow=TRUE)

samp = rbind(s1,s2,s3)

x1 = min(samp[,1])-2
x2 = max(samp[,1])+2
y1 = min(samp[,2])-2
y2 = max(samp[,2])+2

plot(s1,xlim=c(x1,x2),ylim=c(y1,y2))
points(s2, col="blue")
points(s3, col="red")

z.ini = z.true

library(gdata)
library(mclust)

my.result = ME.VVV.err(samp, z.ini, err)

mc.result = meVVV(samp,z.ini)

save(my.result,mc.result,file="sim1_result.RData")
