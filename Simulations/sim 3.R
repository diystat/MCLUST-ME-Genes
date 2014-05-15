
### Simulation 3: varying error structures

library(MASS)
library(gdata)
library(mclust)

## Structure 1: identical errors across observations

# set sample size
nvec = c(3,3,4) * 10
n = sum(nvec)
p = 2

set.seed(0)
mu1 = c(1,1)
mu2 = c(10,-10)
mu3 = c(-25,25)

# set errors to be zero
err = array(0, dim=c(p,p,n))  
for(i in 1:n){
  err[,,i] = matrix(0,p,p)
  diag(err[,,i]) = 20
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
z.ini = z.true = matrix(temp, nrow=n, byrow=TRUE)

samp = rbind(s1,s2,s3)

x1 = min(samp[,1])-2
x2 = max(samp[,1])+2
y1 = min(samp[,2])-2
y2 = max(samp[,2])+2

plot(s1,xlim=c(x1,x2),ylim=c(y1,y2),xlab="",ylab="")
points(s2, col="blue")
points(s3, col="red")

my.result = ME.VVV.err(samp, z.ini, err)
my.class = my.result$z
my.mcr = MCR(my.class,z.true)

mc.result = meVVV(samp,z.ini)
mc.class = mc.result$z
mc.mcr = MCR(mc.class,z.true)


#--------------------------------------------------------------------#


## Structure 2: errors the same within clusters

err = array(0, dim=c(p,p,n))
err1 = matrix(c(9,1,1,9),nrow=2)
err2 = matrix(c(1,0,0,1),nrow=2)
err3 = matrix(c(20,4,4,20),nrow=2)

for(i in 1:nvec[1]){
  err[,,i] = err1
}

for(i in (nvec[1]+1):(n-nvec[3])){
  err[,,i] = err2
}

for(i in (n-nvec[3]+1):n){
  err[,,i] = err3
}

set.seed(0)
s1 = mvrnorm(nvec[1], mu1, (sigma1+err1))
s2 = mvrnorm(nvec[2], mu2, (sigma2+err2))
s3 = mvrnorm(nvec[3], mu3, (sigma3+err3))

# true membership matrix:
temp = c(rep(c(1,0,0),nvec[1]),rep(c(0,1,0),nvec[2]),rep(c(0,0,1),nvec[3]))
z.ini = z.true = matrix(temp, nrow=n, byrow=TRUE)

samp = rbind(s1,s2,s3)

x1 = min(samp[,1])-2
x2 = max(samp[,1])+2
y1 = min(samp[,2])-2
y2 = max(samp[,2])+2

plot(s1,xlim=c(x1,x2),ylim=c(y1,y2),xlab="",ylab="")
points(s2, col="blue")
points(s3, col="red")


my.result = ME.VVV.err(samp, z.ini, err)
my.class = my.result$z
my.mcr = MCR(my.class,z.true)

mc.result = meVVV(samp,z.ini)
mc.class = mc.result$z
mc.mcr = MCR(mc.class,z.true)



#----------------------------------------------------------------------#



### Structure 3: different errors for each observation
set.seed(999)
cc = abs(rnorm(300,4,3))
err = array(0, dim=c(p,p,n))

for(i in 1:n){
  L = lowertriangle(cc[(3*n-2):(3*n)],diag=TRUE)
  err[,,i] = tcrossprod(L)
}

samp = matrix(n,p)
for(i in 1:nvec[1]){
  samp[i] = mvrnorm(1, mu1, (sigma1+err[,,i]))
}

for(i in (nvec[1]+1):(n-nvec[3])){
  samp[i] = mvrnorm(1, mu2, (sigma2+err[,,i]))
}

for(i in (n-nvec[3]):n){
  samp[i] = mvrnorm(1, mu3, (sigma3+err[,,i]))
}


my.result = ME.VVV.err(samp, z.ini, err)
my.class = my.result$z
my.mcr = MCR(my.class,z.true)

mc.result = meVVV(samp,z.ini)
mc.class = mc.result$z
mc.mcr = MCR(mc.class,z.true)










