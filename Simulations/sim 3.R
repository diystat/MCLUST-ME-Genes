
### Simulation 3: varying error structures

library(MASS)
library(gdata)
library(mclust)

## Structure 1: identical errors across observations

# set sample size
nvec = c(3,3,4) * 5
n = sum(nvec)
p = 2
G = 3

tau = c(0.3,0.3,0.4)

set.seed(0)
mu1 = c(1,1)
mu2 = c(10,-10)
mu3 = c(-25,25)
mu = cbind(mu1,mu2,mu3)

# set errors to be zero
err = array(0, dim=c(p,p,n))  
for(i in 1:n){
  err[,,i] = matrix(0,p,p)
  diag(err[,,i]) = 0.1
}

sigma1 = matrix(c(15,-2,-2,15),nrow=2)
sigma2 = matrix(c(23,3,3,23),nrow=2)
sigma3 = matrix(c(31,-4,-4,31),nrow=2)

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



tmp = proc.time()
my.result1 = ME.VVV.err(samp, my.class, err, errstr="identical")
my.class1 = my.result1$z
my.mcr1 = MCR(my.class,z.true)
my.time1 = proc.time() - tmp

tmp = proc.time()
mc.result1 = meVVV(samp,mc.class)
mc.class1 = mc.result1$z
mc.mcr1 = MCR(mc.class,z.true)
mc.time1 = proc.time() - tmp



#--------------------------------------------------------------------#


## Structure 2: errors the same within clusters

err = array(0, dim=c(p,p,n))
err1 = matrix(c(0.1,0,0,0.1),nrow=2)
err2 = matrix(c(0.3,0,0,0.3),nrow=2)
err3 = matrix(c(0.5,0,0,0.5),nrow=2)

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

tmp = proc.time()
my.result2 = ME.VVV.err(samp, z.ini, err, "cluster")
my.class2 = my.result$z
my.mcr2 = MCR(my.class,z.true)
my.time2 = proc.time() - tmp

tmp = proc.time()
mc.result2 = meVVV(samp,z.ini)
mc.class2 = mc.result$z
mc.mcr2 = MCR(mc.class,z.true)
mc.time2 = proc.time() - tmp




#----------------------------------------------------------------------#



### Structure 3: different errors for each observation
set.seed(999)
m = n*p*(p+1)/2
q = p*(p+1)/2
cc = abs(rnorm(m,1,1))
err = array(0, dim=c(p,p,n))

for(i in 1:n){
  L = matrix(0,p,p)
  lowerTriangle(L,diag=T) = cc[(i*q-q+1):(i*q)]
  err[,,i] = tcrossprod(L)
}

samp = matrix(0,n,p)
for(i in 1:nvec[1]){
  samp[i,] = mvrnorm(1, mu1, (sigma1+err[,,i]))
}

for(i in (nvec[1]+1):(n-nvec[3])){
  samp[i,] = mvrnorm(1, mu2, (sigma2+err[,,i]))
}

for(i in (n-nvec[3]):n){
  samp[i,] = mvrnorm(1, mu3, (sigma3+err[,,i]))
}


plot(s1,xlim=c(x1,x2),ylim=c(y1,y2),xlab="",ylab="")
points(s2, col="blue")
points(s3, col="red")


tmp = proc.time()
my.result3 = ME.VVV.err(samp, z.ini, err, errstr="no")
my.class3 = my.result$z
my.mcr3 = MCR(my.class,z.true)
my.time3 = proc.time() - tmp

tmp = proc.time()
mc.result3 = meVVV(samp,z.ini)
mc.class3 = mc.result$z
mc.mcr3 = MCR(mc.class,z.true)
mc.time3 = proc.time() - tmp

save(my.result1,my.class1,my.mcr1,my.time1,
  mc.result1,mc.class1,mc.mcr1,mc.time1,
  my.result2,my.class2,my.mcr2,my.time2,
  mc.result2,mc.class2,mc.mcr2,mc.time2,
  my.result3,my.class3,my.mcr3,my.time3,
  mc.result3,mc.class3,mc.mcr3,mc.time3, file="sim3_results.RData")






