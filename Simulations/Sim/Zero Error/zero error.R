N = 200
tau = 0.5
mu1 = c(0,0)
mu2 = c(8,0)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
k = 36
p = 0
nseed = 100

check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)

setwd("~/Research Project/R functions for MCME(VVV)/Core functions")

library(gdata)
library(phyclust)

ptm <- proc.time()
sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
proc.time() - ptm



setwd("/Users/wzhang/Research Project/Simulations/Sim/Zero Error/Results")
load("res.RData")
rr = out$rand.raw
r.mcme = rr[,1]
r.mevvv = rr[,4]
fr.mcme = rr[,7]
fr.mevvv = rr[,8]

rand.diff = r.mcme - r.mevvv
frand.diff = fr.mcme - fr.mevvv

par(mfrow=c(1,2))
boxplot(rand.diff,main="Pairwise diff of Rand indices\n (MCME-meVVV)")
abline(0,0,lty="dashed")

boxplot(frand.diff,main="Pairwise diff of Fuzzy Rand indices\n (MCME-meVVV)")
abline(0,0,lty="dashed")
par(mfrow=c(1,1))

