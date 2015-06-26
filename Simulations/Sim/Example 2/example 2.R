N = 200
tau = 0.5
mu1 = c(0,0)
mu2 = c(8,0)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(16,0,0,16),nrow=2)
k = 9
p = 0.3
nseed = 100

library(gdata)
library(phyclust)

ptm <- proc.time()
sim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
proc.time() - ptm


setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 2/p=0.1/Results")
load("res.RData")
rr = out$rand.raw
r.mcme = rr[,1]
r.mevvv = rr[,4]
fr.mcme = rr[,7]
fr.mevvv = rr[,8]



xlow = min(r.mevvv)-0.05
xup = max(r.mevvv)+0.05
ylow = min(r.mcme)-0.05
yup = max(r.mcme)+0.05

prop = sum(r.mcme>r.mevvv)/nseed*100

plot(r.mevvv,r.mcme,xlim=c(xlow,xup),ylim=c(ylow,yup),cex=0.7,
  main=paste("Above y=x: ",prop,"%",sep=""),xlab="meVVV Rand index",ylab="MCME Rand index")
abline(0,1,lty="dashed")
















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



setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 2/p=0.3/Results")
load("res.RData")
rr = out$rand.raw
r.mcme = rr[,1]
r.mevvv = rr[,4]
fr.mcme = rr[,7]
fr.mevvv = rr[,8]



xlow = min(r.mevvv)-0.05
xup = max(r.mevvv)+0.05
ylow = min(r.mcme)-0.05
yup = max(r.mcme)+0.05

prop = sum(r.mcme>r.mevvv)/nseed*100

plot(r.mevvv,r.mcme,xlim=c(xlow,xup),ylim=c(ylow,yup),cex=0.7,
  main=paste("Above y=x: ",prop,"%",sep=""),xlab="meVVV Rand index",ylab="MCME Rand index")
abline(0,1,lty="dashed")























c(min(r.mcme),median(r.mcme),mean(r.mcme),max(r.mcme))
c(min(r.mevvv),median(r.mevvv),mean(r.mevvv),max(r.mevvv))

c(min(fr.mcme),median(fr.mcme),mean(fr.mcme),max(fr.mcme))
c(min(fr.mevvv),median(fr.mevvv),mean(fr.mevvv),max(fr.mevvv))


group = c(rep(1,100),rep(2,100))
y = c(r.mcme,r.mevvv)
fy = c(fr.mcme,fr.mevvv)
par(mfrow=c(1,2))
boxplot(y~group,names=c("MCME","meVVV"),main="Rand index, 30%")
boxplot(fy~group,names=c("MCME","meVVV"),main="Fuzzy Rand index, 30%")
par(mfrow=c(1,1))



setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 2/p=0.5/Results")
load("res.RData")
rr = out$rand.raw
r.mcme = rr[,1]
r.mevvv = rr[,4]
fr.mcme = rr[,7]
fr.mevvv = rr[,8]


xlow = min(r.mevvv)-0.05
xup = max(r.mevvv)+0.05
ylow = min(r.mcme)-0.05
yup = max(r.mcme)+0.05

prop = sum(r.mcme>r.mevvv)/nseed*100

plot(r.mevvv,r.mcme,xlim=c(xlow,xup),ylim=c(ylow,yup),cex=0.7,
  main=paste("Above y=x: ",prop,"%",sep=""),xlab="meVVV Rand index",ylab="MCME Rand index")
abline(0,1,lty="dashed")









c(min(r.mcme),median(r.mcme),mean(r.mcme),max(r.mcme))
c(min(r.mevvv),median(r.mevvv),mean(r.mevvv),max(r.mevvv))

c(min(fr.mcme),median(fr.mcme),mean(fr.mcme),max(fr.mcme))
c(min(fr.mevvv),median(fr.mevvv),mean(fr.mevvv),max(fr.mevvv))


group = c(rep(1,100),rep(2,100))
y = c(r.mcme,r.mevvv)
fy = c(fr.mcme,fr.mevvv)
par(mfrow=c(1,2))
boxplot(y~group,names=c("MCME","meVVV"),main="Rand index, 50%")
boxplot(fy~group,names=c("MCME","meVVV"),main="Fuzzy Rand index, 50%")
par(mfrow=c(1,1))



setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 2/p=0.7/Results")
load("res.RData")
rr = out$rand.raw
r.mcme = rr[,1]
r.mevvv = rr[,4]
fr.mcme = rr[,7]
fr.mevvv = rr[,8]



xlow = min(r.mevvv)-0.05
xup = max(r.mevvv)+0.05
ylow = min(r.mcme)-0.05
yup = max(r.mcme)+0.05

prop = sum(r.mcme>r.mevvv)/nseed*100

plot(r.mevvv,r.mcme,xlim=c(xlow,xup),ylim=c(ylow,yup),cex=0.7,
  main=paste("Above y=x: ",prop,"%",sep=""),xlab="meVVV Rand index",ylab="MCME Rand index")
abline(0,1,lty="dashed")














c(min(r.mcme),median(r.mcme),mean(r.mcme),max(r.mcme))
c(min(r.mevvv),median(r.mevvv),mean(r.mevvv),max(r.mevvv))

c(min(fr.mcme),median(fr.mcme),mean(fr.mcme),max(fr.mcme))
c(min(fr.mevvv),median(fr.mevvv),mean(fr.mevvv),max(fr.mevvv))


group = c(rep(1,100),rep(2,100))
y = c(r.mcme,r.mevvv)
fy = c(fr.mcme,fr.mevvv)
par(mfrow=c(1,2))
boxplot(y~group,names=c("MCME","meVVV"),main="Rand index, 70%")
boxplot(fy~group,names=c("MCME","meVVV"),main="Fuzzy Rand index, 70%")
par(mfrow=c(1,1))



setwd("/Users/wzhang/Research Project/Simulations/Sim/Example 2/p=0.9/Results")
load("res.RData")
rr = out$rand.raw
r.mcme = rr[,1]
r.mevvv = rr[,4]
fr.mcme = rr[,7]
fr.mevvv = rr[,8]



xlow = min(r.mevvv)-0.05
xup = max(r.mevvv)+0.05
ylow = min(r.mcme)-0.05
yup = max(r.mcme)+0.05

prop = sum(r.mcme>r.mevvv)/nseed*100

plot(r.mevvv,r.mcme,xlim=c(xlow,xup),ylim=c(ylow,yup),cex=0.7,
  main=paste("Above y=x: ",prop,"%",sep=""),xlab="meVVV Rand index",ylab="MCME Rand index")
abline(0,1,lty="dashed")



















c(min(r.mcme),median(r.mcme),mean(r.mcme),max(r.mcme))
c(min(r.mevvv),median(r.mevvv),mean(r.mevvv),max(r.mevvv))

c(min(fr.mcme),median(fr.mcme),mean(fr.mcme),max(fr.mcme))
c(min(fr.mevvv),median(fr.mevvv),mean(fr.mevvv),max(fr.mevvv))


group = c(rep(1,100),rep(2,100))
y = c(r.mcme,r.mevvv)
fy = c(fr.mcme,fr.mevvv)
par(mfrow=c(1,2))
boxplot(y~group,names=c("MCME","meVVV"),main="Rand index, 90%")
boxplot(fy~group,names=c("MCME","meVVV"),main="Fuzzy Rand index, 90%")
par(mfrow=c(1,1))


























data.par = sim.par(N,tau,mu1,mu2,sig1,sig2,k,p)
simdata = sim.data(data.par)

rand.samples = simdata$data
group = simdata$z.ini[,1]

col.blue = rgb(30,144,255,max=255)  
xlow = min(rand.samples[,1])-0.5
xup = max(rand.samples[,1])+0.5
ylow = min(rand.samples[,2])-0.2
yup = max(rand.samples[,2])+0.2

plot(rand.samples[,1],rand.samples[,2],col=ifelse(group==0,"red",col.blue))

sum(group)


