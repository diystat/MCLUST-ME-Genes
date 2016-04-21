

library(MASS)
library(mclust)
library(gdata)
library(phyclust)

N = 200
tau = 0.5
mu1 = c(0,0)
mu2 = c(8,0)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
k = 9
p = 0.1
nseed = 100

#check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)

setwd("Z://Research Project/simulations/Sim/example 3/p=0.1")

ptm <- proc.time()
sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
t01 = proc.time() - ptm






N = 200
tau = 0.5
mu1 = c(0,0)
mu2 = c(8,0)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
k = 9
p = 0.3
nseed = 100

#check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)
  
setwd("Z://Research Project/simulations/Sim/example 3/p=0.3")

ptm <- proc.time()
sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
t03 = proc.time() - ptm






N = 200
tau = 0.5
mu1 = c(0,0)
mu2 = c(8,0)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
k = 9
p = 0.7
nseed = 100

#check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)

setwd("Z://Research Project/simulations/Sim/example 3/p=0.7")

ptm <- proc.time()
sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
t07 = proc.time() - ptm






N = 200
tau = 0.5
mu1 = c(0,0)
mu2 = c(8,0)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
k = 9
p = 0.9
nseed = 100

#check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)

setwd("Z://Research Project/simulations/Sim/example 3/p=0.9")

ptm <- proc.time()
sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
t09 = proc.time() - ptm






## Examine change in clustering uncertainties
plot.uncr = function(seed){
  setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 3/p=0.5/Results")
  load("res.RData")

  # Extract clustering result for seed="seed"
  ivec = rep((seed-1)*7,7)+(1:7)

  res.mcme = out$sim.result[[ivec[1]]]
  res.mevvv = out$sim.result[[ivec[2]]]
  z.ini = out$sim.result[[ivec[3]]]
  rand.samples = out$sim.result[[ivec[4]]]
  errmat = out$sim.result[[ivec[5]]]
  index = out$sim.result[[ivec[6]]]
  k = out$sim.result[[ivec[7]]]
  simrun = list(res.mcme=res.mcme,res.mevvv=res.mevvv,errmat=errmat,z.ini=z.ini,rand.samples=rand.samples,index=index,k=k)

  mcme7 = res.mcme
  mclust7 = res.mevvv
  unc.mcme = mcme7$uncertainty
  unc.mclust = numeric(200)
  for(i in 1:200){
    unc.mclust[i] = 1-max(mclust7$z[i,])
  }

  # Ratio of uncertainties
  uncr = log((unc.mcme+0.001)/(unc.mclust+0.001))

  # Compare ratio with measurement error
  # fr.new = data.frame(ratio=uncr,index=simdata$index)
  # boxplot(ratio~index,data=fr.new)

  #point.size = numeric(200)
  #for(i in 1:200){
  #  if(uncr[i]>log(2)){
  #    point.size[i] = 3
  #  }else if(uncr[i]<log(1/2)){
  #    point.size[i] = 0.7
  #  }else{
  #    point.size[i] = 1.3
  #  }
  #}
  
  point.size = 0.3+unc.mcme*4
  point.size.mclust = 0.3+unc.mclust*4

  plot.boundary.new(simrun,point.size,point.size.mclust)
}






# Read simulation results
setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 3/p=0.1/Results")
load("res.RData")
rr = out$rand.raw
r1.mcme = rr[,1]
r1.mevvv = rr[,4]
fr1.mcme = rr[,7]
fr1.mevvv = rr[,8]

setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 3/p=0.3/Results")
load("res.RData")
rr = out$rand.raw
r3.mcme = rr[,1]
r3.mevvv = rr[,4]
fr3.mcme = rr[,7]
fr3.mevvv = rr[,8]

setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 3/p=0.5/Results")
load("res.RData")
rr = out$rand.raw
r5.mcme = rr[,1]
r5.mevvv = rr[,4]
fr5.mcme = rr[,7]
fr5.mevvv = rr[,8]

setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 3/p=0.7/Results")
load("res.RData")
rr = out$rand.raw
r7.mcme = rr[,1]
r7.mevvv = rr[,4]
fr7.mcme = rr[,7]
fr7.mevvv = rr[,8]

setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 3/p=0.9/Results")
load("res.RData")
rr = out$rand.raw
r9.mcme = rr[,1]
r9.mevvv = rr[,4]
fr9.mcme = rr[,7]
fr9.mevvv = rr[,8]



rand.mcme = c(r1.mcme,r3.mcme,r5.mcme,r7.mcme,r9.mcme)
rand.mevvv = c(r1.mevvv,r3.mevvv,r5.mevvv,r7.mevvv,r9.mevvv)
frand.mcme = c(fr1.mcme,fr3.mcme,fr5.mcme,fr7.mcme,fr9.mcme)
frand.mevvv = c(fr1.mevvv,fr3.mevvv,fr5.mevvv,fr7.mevvv,fr9.mevvv)

group = c(rep(0.1,100),rep(0.3,100),rep(0.5,100),rep(0.7,100),rep(0.9,100))

# Individual Rand/Fuzzy Rand indices
par(mfrow=c(2,2))
boxplot(rand.mcme~group,main="Rand, MCLUST-ME",xlab="proportion of obs with errors")
boxplot(rand.mevvv~group,main="Rand, MCLUST",xlab="proportion of obs with errors")
boxplot(frand.mcme~group,main="Fuzzy Rand, MCLUST-ME",xlab="proportion of obs with errors")
boxplot(frand.mevvv~group,main="Fuzzy Rand, MCLUST",xlab="proportion of obs with errors")
par(mfrow=c(1,1))


# Pairwise difference for all p
rand.diff = rand.mcme - rand.mevvv
frand.diff = frand.mcme - frand.mevvv

par(mfrow=c(2,2))
boxplot(rand.diff~group,main="Rand index",xlab="proportion of obs with errors")
abline(0,0,lty="dashed")
boxplot(frand.diff~group,main="Fuzzy Rand index",xlab="proportion of obs with errors")
abline(0,0,lty="dashed")

# Pairwise comparison for p=0.5
plot(r5.mevvv,r5.mcme-r5.mevvv,xlab="MCLUST",ylab="MCLUST-ME--MCLUST",main="Rand index, p = 0.5")
abline(0,0,lty="dashed")

plot(fr5.mevvv,fr5.mcme-fr5.mevvv,xlab="MCLUST",ylab="MCLUST-ME--MCLUST",main="Fuzzy Rand index, p = 0.5")
abline(0,0,lty="dashed")

par(mfrow=c(1,1))



