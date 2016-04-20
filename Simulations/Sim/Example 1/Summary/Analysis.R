
## Analysis of results from Example 4:

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

group = c(rep(0.1,100),rep(0.3,100),rep(0.5,100),rep(0.7,100),rep(0.9,100))

## Boxplots for individual rand indices:
par(mfrow=c(2,2))
boxplot(rand.mcme~group,xlab="proportion of obs with errors",main="Rand index, MCME")
boxplot(rand.mevvv~group,xlab="proportion of obs with errors",main="Rand index, meVVV")
boxplot(frand.mcme~group,xlab="proportion of obs with errors",main="Fuzzy Rand index, MCME")
boxplot(frand.mevvv~group,xlab="proportion of obs with errors",main="Fuzzy Rand index, meVVV")
par(mfrow=c(1,1))


## Boxplots for pairwise differences of rand indices:
rand.diff = rand.mcme - rand.mevvv
frand.diff = frand.mcme - frand.mevvv

par(mfrow=c(1,2))
boxplot(rand.diff~group,xlab="proportion of obs with errors",main="Pairwise diff of Rand indices\n (MCME-meVVV)")
abline(0,0,lty="dashed")

boxplot(frand.diff~group,xlab="proportion of obs with errors",main="Pairwise diff of Fuzzy Rand indices\n (MCME-meVVV)")
abline(0,0,lty="dashed")
par(mfrow=c(1,1))

