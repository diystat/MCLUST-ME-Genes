## specify sim parameters:
N = 100
tau = 0.5
mu1 = c(0,0)
mu2 = c(0,0)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
k = 9
p = 0.5
nseed = 100

library(gdata)
library(phyclust)

ptm <- proc.time()
sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
proc.time() - ptm



setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 4/p=0.1/Results")
load("res.RData")
rr = out$rand.raw
r.mcme = rr[,1]
r.mevvv = rr[,4]
fr.mcme = rr[,7]
fr.mevvv = rr[,8]

c(min(r.mcme),median(r.mcme),mean(r.mcme),max(r.mcme))
c(min(r.mevvv),median(r.mevvv),mean(r.mevvv),max(r.mevvv))

c(min(fr.mcme),median(fr.mcme),mean(fr.mcme),max(fr.mcme))
c(min(fr.mevvv),median(fr.mevvv),mean(fr.mevvv),max(fr.mevvv))


group = c(rep(1,100),rep(2,100))
y = c(r.mcme,r.mevvv)
fy = c(fr.mcme,fr.mevvv)
par(mfrow=c(1,2))
boxplot(y~group,names=c("MCME","meVVV"),main="Rand index, 10%")
boxplot(fy~group,names=c("MCME","meVVV"),main="Fuzzy Rand index, 10%")
par(mfrow=c(1,1))

