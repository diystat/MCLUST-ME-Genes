
### Simulation 2: Increasing, but identical and diagonal errors
library(MASS)
library(gdata)
library(mclust)

# set sample size
nvec = c(3,3,4) * 10
n = sum(nvec)
p = 2

mu1 = c(1,1)
mu2 = c(5,-5)
mu3 = c(-25,25)

# define cov matrices
  sigma1 = matrix(c(15,-2,-2,15),nrow=2)
  sigma2 = matrix(c(23,3,3,23),nrow=2)
  sigma3 = matrix(c(31,-4,-4,31),nrow=2)


# set error value here:
eps = c(0, 5, 10, 15, 20, 30)

# set error matrix:
err = array(0, dim=c(p, p, length(eps)))
for(i in 1:length(eps)){
  diag(err[,,i]) = eps[i]
}

err.array = array(0, dim=c(p,p,n,length(eps)))
for(i in 1:length(eps)){  
  for(j in 1:n){
    err.array[,,j,i] = err[,,i]
  } 
}


# true membership matrix:
temp = c(rep(c(1,0,0),nvec[1]),rep(c(0,1,0),nvec[2]),rep(c(0,0,1),nvec[3]))
z.ini = z.true = matrix(temp, nrow=n, byrow=TRUE)


my.mcr = mc.mcr = rep(0,length(eps))
my.time = mc.time = list()



par(mfrow=c(2,3))
for(i in 1:length(eps)){  
  # sample from dist:
  set.seed(0)
  s1 = mvrnorm(nvec[1], mu1, (sigma1+err[,,i]))
  s2 = mvrnorm(nvec[2], mu2, (sigma2+err[,,i]))
  s3 = mvrnorm(nvec[3], mu3, (sigma3+err[,,i]))
  
  samp = rbind(s1,s2,s3)
  
  x1 = min(samp[,1])-2
  x2 = max(samp[,1])+2
  y1 = min(samp[,2])-2
  y2 = max(samp[,2])+2
  
  plot(s1,xlim=c(x1,x2),ylim=c(y1,y2),main = paste(expression(epsilon), " = ", eps[i]),xlab="",ylab="")
  points(s2, col="blue")
  points(s3, col="red")  
  
  tmp = proc.time()
  my.class = ME.VVV.err(samp, z.ini, err.array[,,,i])$z
  my.mcr[i] = MCR(my.class, z.true)
  my.time[[i]] = proc.time() - tmp
  
  tmp = proc.time()
  mc.class = meVVV(samp,z.ini)$z
  mc.mcr[i] = MCR(mc.class, z.true)
  mc.time[[i]] = proc.time() - tmp
  
}
par(mfrow=c(1,1))

save(my.mcr, mc.mcr, my.time, mc.time, file="sim2_result.RData")


