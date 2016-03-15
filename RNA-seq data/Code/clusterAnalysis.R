rm(list=ls())
dd = "2015-11-06"
file.out = sprintf("%s.Coaker.Rdata", dd)
print(load(file.out))

obs = full$beta[top, 1:5]
errary = v.beta[,,top]
betasd = sd.beta[top,]

save(obs,errary,file="RNASeq.RData")

setwd("../R functions for MCLUST-ME(VVV)/Core functions")

### Use entire dataset ###

# Group data with mclust first:
res.mclust = Mclust(obs,modelNames="VVV")
summary(res.mclust)

# Generate initial membership matrix with hierarchical agglomeration
N = nrow(obs)
G = res.mclust$G
hcTree = mclust::hc(obs)
cl = mclust::hclass(hcTree, G)
z.ini = matrix(0,N,G)
  for(i in 1:N){
    for(j in 1:G){
      z.ini[i,j] = ifelse(cl[i]==j,1,0)
    }
  }

# Group data with MCME:
ptm <- proc.time()
res.mcme.full = mcmeVVV(obs,z.ini,errary)
proc.time() - ptm
#      user    system   elapsed 
# 75970.298   987.514 77643.007 

# Write result to file
save(res.mcme.full,file="realDataClusterResult.RData")



###------ Analysis -----###

load("realDataClusterResult.RData")

# Predicted class labels:
zhat = res.mcme.full$z
predClass = numeric(1000)
for(i in 1:1000){
  temp = zhat[i,]
  predClass[i] = which(temp==max(temp))
}

zhat.mclust = res.mclust$z
predClass.mclust = numeric(1000)
for(i in 1:1000){
  temp = zhat.mclust[i,]
  predClass.mclust[i] = which(temp==max(temp))
}


# Confusion matrix:
library(caret)
confusionMatrix(predClass,predClass.mclust,dnn=c("MCME","MCLUST"))

# Adjusted rand index:
adjustedRandIndex(predClass,predClass.mclust)

# Scatterplot matrix:
pairs(obs,col=predClass,main="MCME")
pairs(obs,col=predClass.mclust,main="MCLUST")

# Line plots for each predicted cluster:
lineplot = function(group,method="MCLUST"){
  x = 1:5
  groupind.mcme = which(predClass==group)
  groupind.mclust = which(predClass.mclust==group)
  
  if(method=="MCLUST"){
    data = obs[predClass.mclust==group,]
    label = "MCLUST"
    redind = setdiff(groupind.mclust,groupind.mcme)
  }else{
    data = obs[predClass==group,]
    label = "MCLUST-ME"
    redind = setdiff(groupind.mcme,groupind.mclust)
  }
  
  ymin = min(data)
  ymax = max(data)
  plot(x,data[1,],xaxt="n",type="l",ylim=c(ymin,ymax),col="lightgreen",xlab="",ylab="Log2 fold change",
    main=paste(label," Group ",group,sep=""),lwd=0.5)
  axis(1,at=1:5,labels=c("10min","1h","3h","6h","12h"))
  for(i in 2:nrow(data)){
    lines(x,data[i,],col="lightgreen",lwd=0.5)
  }
  # Use red lines to show differently clustered observations
  for(j in 1:length(redind)){
    lines(x,obs[redind[j],],col="red",lwd=0.6)
  }
}

par(mfrow=c(2,3))
lineplot(1,"")
lineplot(2,"")
lineplot(3,"")
lineplot(1)
lineplot(2)
lineplot(3)
par(mfrow=c(1,1))



# Plot mean and variance estimates:
param = res.mcme.full$parameters
names(param)
mu.mcme = param$muhat
sigma.mcme = param$sigmahat

param = res.mclust$parameters
mu.mclust = param$mean
sigma.mclust = param$variance$sigma

meanLinePlot = function(group,legend=FALSE){
  x = 1:5
  mu1.mcme = mu.mcme[,group]
  mu1.mclust = mu.mclust[,group]
  ymin = min(c(mu1.mcme,mu1.mclust))
  ymax = max(c(mu1.mcme,mu1.mclust))
  plot(x,mu1.mcme,xaxt='n',type="l",main=paste("Group ",group," mean",sep=""),
    col="blue",lwd=2,xlab="",ylim=c(ymin,ymax),ylab="Log2 fold change")
  axis(1,at=1:5,labels=c("10min","1h","3h","6h","12h"))
  lines(x,mu1.mclust,lty="dashed",col="red",lwd=2)
  if(legend==TRUE){
    legend("bottomright",legend=c("MCLUST-ME","MCLUST"),col=c("blue","red"),
      lty=c("solid","dashed"),lwd=c(2,2))
  }
}

par(mfrow=c(1,3))
meanLinePlot(1)
meanLinePlot(2)
meanLinePlot(3,TRUE)
par(mfrow=c(1,1))








