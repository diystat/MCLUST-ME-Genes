N = 200
tau = 0.5
mu1 = c(0,0)
mu2 = c(8,0)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
k = 36
p = 1
nseed = 100

check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)

setwd("~/Research Project/R functions for MCME(VVV)/Core functions")

library(gdata)
library(phyclust)

ptm <- proc.time()
sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
proc.time() - ptm


setwd("/Users/wzhang/Research Project/Simulations/Sim/Constant Error/Results")
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

hist(rand.diff)
hist(frand.diff)

sum(rand.diff==0)/100
sum(frand.diff==0)/100

summary(rand.diff)
summary(frand.diff)

t.test(rand.diff,mu=0)
t.test(frand.diff,mu=0)




out$seed
names(out)
names(out$sim.result)
rmcme = out$sim.result[505]$res.mcme
rmevvv = out$sim.result[506]$res.mevvv

names(rmcme)
names(rmevvv)

pmcme = rmcme$par
names(pmcme)
pmcme$sigmahat

pmevvv = rmevvv$par
names(pmevvv)
pmevvv$var













