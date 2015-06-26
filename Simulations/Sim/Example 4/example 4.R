N = 200
tau = 0.5
mu1 = c(0,0)
mu2 = c(8,0)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
k = 36
p = 0.3
nseed = 100

check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)

setwd("~/Research Project/R functions for MCME(VVV)/Core functions")

library(gdata)
library(phyclust)

ptm <- proc.time()
sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
proc.time() - ptm





setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 4/p=0.1/Results")
load("res.RData")
rr = out$rand.raw
r1.mcme = rr[,1]
r1.mevvv = rr[,4]
fr1.mcme = rr[,7]
fr1.mevvv = rr[,8]

setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 4/p=0.3/Results")
load("res.RData")
rr = out$rand.raw
r3.mcme = rr[,1]
r3.mevvv = rr[,4]
fr3.mcme = rr[,7]
fr3.mevvv = rr[,8]

setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 4/p=0.5/Results")
load("res.RData")
rr = out$rand.raw
r5.mcme = rr[,1]
r5.mevvv = rr[,4]
fr5.mcme = rr[,7]
fr5.mevvv = rr[,8]

setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 4/p=0.7/Results")
load("res.RData")
rr = out$rand.raw
r7.mcme = rr[,1]
r7.mevvv = rr[,4]
fr7.mcme = rr[,7]
fr7.mevvv = rr[,8]

setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 4/p=0.9/Results")
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

group = c(rep(1,100),rep(3,100),rep(5,100),rep(7,100),rep(9,100))

boxplot(rand.mcme~group)
boxplot(rand.mevvv~group)
boxplot(frand.mcme~group)
boxplot(frand.mevvv~group)

rand.diff = rand.mcme - rand.mevvv
frand.diff = frand.mcme - frand.mevvv

boxplot(rand.diff~group)
abline(0,0,lty="dashed")

boxplot(frand.diff~group)
abline(0,0,lty="dashed")












