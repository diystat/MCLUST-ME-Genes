setwd("~/Research Project/RNA-seq data")
rm(list=ls())
dd = "2015-11-06"
file.out = sprintf("%s.Coaker.Rdata", dd)
print(load(file.out))

obs = full$beta[top, 1:5]
errary = v.beta[,,top]
betasd = sd.beta[top,]

save(obs,errary,file="RNASeq.RData")

setwd("../R functions for MCLUST-ME(VVV)/Core functions")
source(file="all functions.R")

### Use entire dataset ###

# Group data with mclust first:
library(mclust)
res.mclust = Mclust(obs,modelNames="VVV")
summary(res.mclust)

# Define a wrapper function for MCLUST-ME:
mclust.run = function(G){
  # Generate initial membership matrix with hierarchical agglomeration
  N = nrow(obs)
  hcTree = hc(obs)
  cl = hclass(hcTree, G)
  z.ini = matrix(0,N,G)
  for(i in 1:N){
    for(j in 1:G){
      z.ini[i,j] = ifelse(cl[i]==j,1,0)
    }
  }
  
  # Group data with MCLUST-ME:
  res.mcme.full = mcmeVVV(obs,z.ini,errary)
  return(res.mcme.full)  
}

# Run MCLUST-ME assuming 1~8 clusters:
run1 = mclust.run(1)
run2 = mclust.run(2)
run3 = mclust.run(3)
run4 = mclust.run(4)
run5 = mclust.run(5)
run6 = mclust.run(6)
run7 = mclust.run(7)
run8 = mclust.run(8)

# Plot BIC against number of components:
b = numeric(8)
b[1] = run1$BIC
b[2] = run2$BIC
b[3] = run3$BIC
b[4] = run4$BIC
b[5] = run5$BIC
b[6] = run6$BIC
b[7] = run7$BIC
b[8] = run8$BIC
setwd("~/Research Project/RNA-seq data/Results/BIC")
save(b,file="rnaseq_bic.RData")

setwd("~/Research Project/RNA-seq data/Results/BIC")
load("rnaseq_bic.RData")
a = res.mclust$BIC[1:8]
ylow = min(c(a,b))
yup = max(c(a,b))
plot(1:8,b,type="b",lwd=2,xlab="Number of components",ylab="BIC",pch=0,
  ylim=c(ylow,yup))
lines(1:8,a,type="b",lwd=2,pch=1,lty=2)
legend("bottomright",legend=c("MCLUST-ME","MCLUST"),
    pch=c(0,2),lty=c("solid","dashed"),lwd=2)

# BIC is greatest with 3 components, so keep this result:
res.mcme.full = run3

# Write result to file
save(res.mcme.full,file="realDataClusterResult.RData")



###------ Analysis -----###

setwd("~/Research Project/RNA-seq data/Results")
load("realDataClusterResult.RData")

## Obtain mean and covariance estimates
param = res.mcme.full$parameters
names(param)
mu.mcme = param$muhat
sigma.mcme = param$sigmahat

param = res.mclust$parameters
mu.mclust = param$mean
sigma.mclust = param$variance$sigma


### Predicted class labels:
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


### Which genes are clustered differently?
index_diff = which(predClass!=predClass.mclust)
mem.table = round(cbind(zhat[index_diff,],zhat.mclust[index_diff,]),4)
library(xtable)
xtable(mem.table,digits=4)
unc = res.mcme.full$uncertainty
unc_mclust = res.mclust$uncertainty
new = data.frame(unc=unc[index_diff],unc_mclust=unc_mclust[index_diff],
  mclust_me=predClass[index_diff],mclust=predClass.mclust[index_diff])

library(dplyr)
new = new %>% 
  arrange(mclust_me,mclust) %>%
  mutate(o1=10^(floor(log10(unc))),o2=10^(floor(log10(unc_mclust))))

filter(new,o1!=o2) # genes with uncertainties of different order of magnitude





### Confusion matrix:
library(caret)
confusionMatrix(predClass,predClass.mclust,dnn=c("MCME","MCLUST"))

### Rand index:
phyclust::RRand(predClass,predClass.mclust)

### Adjusted rand index:
adjustedRandIndex(predClass,predClass.mclust)

### Raw data scatterplot matrix:
pairs(obs,labels=c("10min","1h","3h","6h","12h"),pch=16,cex=0.6)

### Partitioned scatterplot matrix (colorized)
pairs(obs,main="MCLUST-ME",col=c("magenta3","cyan3","orange")[predClass],
  labels=c("10min","1h","3h","6h","12h"),pch=16,cex=0.6)
pairs(obs,main="MCLUST",col=c("magenta3","cyan3","orange")[predClass.mclust],
  labels=c("10min","1h","3h","6h","12h"),pch=16,cex=0.6)

### Partitioned scatterplot matrix (black and white)
pairs(obs,labels=c("10min","1h","3h","6h","12h"),pch=c(1,2,3)[predClass],
  cex=0.6,col="#00000033") # MCLUST-ME
pairs(obs,labels=c("10min","1h","3h","6h","12h"),pch=c(1,2,3)[predClass.mclust],
  cex=0.6) # MCLUST

### Partitioned scatterplot matrix (black and white,with uncertainty)
unc = res.mcme.full$uncertainty
pointsize = 0.3+unc*4 
pairs(obs,labels=c("10min","1h","3h","6h","12h"),pch=c(1,2,3)[predClass],
  cex=pointsize) # MCLUST-ME
pairs(obs,labels=c("10min","1h","3h","6h","12h"),pch=c(1,2,3)[predClass.mclust],
  cex=pointsize) # MCLUST

### Partitioned data,with 1h and 3h only
par(mfrow=c(2,2))
plot(1:8,b,main="BIC of MCLUST-ME",type="b",pch=0,lwd=2,xlab="Number of components",ylab="BIC")
plot(1:8,a,main="BIC of MCLUST",type="b",pch=0,lwd=2,xlab="Number of components",ylab="BIC")
plot(obs[,2:3],main="MCLUST-ME",pch=c(1,2,3)[predClass],cex=0.7,xlab="1h",ylab="3h")
#abline(v=0,lty="dashed")
#abline(h=0,lty="dashed")
plot(obs[,2:3],main="MCLUST",pch=c(1,2,3)[predClass.mclust],cex=0.7,xlab="1h",ylab="3h")
#abline(v=0,lty="dashed")
#abline(h=0,lty="dashed")
par(mfrow=c(1,1))


### Partitioned data on 1h and 3h with confidence outlines
library(ellipse)
# Generate coordinates for confidence outlines(MCLUST-ME)
g1.coord = ellipse(sigma.mcme[2:3,2:3,1],centre=mu.mcme[2:3,1],level=0.9)
g2.coord = ellipse(sigma.mcme[2:3,2:3,2],centre=mu.mcme[2:3,2],level=0.9)
g3.coord = ellipse(sigma.mcme[2:3,2:3,3],centre=mu.mcme[2:3,3],level=0.9)
# Generate coordinates for confidence outlines(MCLUST)
g1.coord.mc = ellipse(sigma.mclust[2:3,2:3,1],centre=mu.mclust[2:3,1],level=0.9)
g2.coord.mc = ellipse(sigma.mclust[2:3,2:3,2],centre=mu.mclust[2:3,2],level=0.9)
g3.coord.mc = ellipse(sigma.mclust[2:3,2:3,3],centre=mu.mclust[2:3,3],level=0.9)

par(mfrow=c(2,2))
plot(1:8,b,main="BIC of MCLUST-ME",type="b",pch=0,lwd=2,xlab="Number of components",ylab="BIC")
#plot(1:8,a,main="BIC of MCLUST",type="b",pch=0,lwd=2,xlab="Number of components",ylab="BIC")
lines(1:8,a,lty="dashed",type="b",pch=0,lwd=2)

plot(obs[,2],obs[,3],main="MCLUST-ME",pch=c(1,2,3)[predClass],cex=0.7,
  xlab="1h",ylab="3h",col="#00000070")
# Plot the confidence outlines
lines(g1.coord,type="l")
lines(g2.coord,type="l")
lines(g3.coord,type="l")

plot(obs[,2],obs[,3],main="MCLUST",pch=c(1,2,3)[predClass.mclust],cex=0.7,
  xlab="1h",ylab="3h",col="#00000070")
# Plot the confidence outlines
lines(g1.coord.mc,type="l")
lines(g2.coord.mc,type="l")
lines(g3.coord.mc,type="l")
lines(g1.coord,type="l",lty="dashed")
lines(g2.coord,type="l",lty="dashed")
lines(g3.coord,type="l",lty="dashed")
par(mfrow=c(1,1))







### Plot mean and variance estimates:
meanLinePlot = function(group,legend=FALSE){
  x = 1:5
  mu1.mcme = mu.mcme[,group]
  mu1.mclust = mu.mclust[,group]
  ymin = min(c(mu.mcme,mu.mclust))
  ymax = max(c(mu.mcme,mu.mclust))
  plot(x,mu1.mcme,xaxt='n',type="l",main=paste("Group ",group," mean",sep=""),
    lwd=2,xlab="",ylim=c(ymin,ymax),ylab="log 2 fold change")
  axis(1,at=1:5,labels=c("10min","1h","3h","6h","12h"))
  lines(x,mu1.mclust,lty="dashed",lwd=2)
  abline(0,0,lty="dotted")
  if(legend==TRUE){
    legend("topright",legend=c("MCLUST-ME","MCLUST"),
      lty=c("solid","dashed"),lwd=c(2,2))
  }
}

par(mfrow=c(1,3))
meanLinePlot(1)
meanLinePlot(2)
meanLinePlot(3,TRUE)
par(mfrow=c(1,1))



### Line plot of grouped genes
g1 = obs[predClass==1,]
g2 = obs[predClass==2,]
g3 = obs[predClass==3,]
x = 1:5

ymin1 = min(g1)
ymax1 = max(g1)
ymin2 = min(g2)
ymax2 = max(g2)
ymin3 = min(g3)
ymax3 = max(g3)

mcme.g1 = which(predClass==1)
mcme.g2 = which(predClass==2)
mcme.g3 = which(predClass==3)

mclust.g1 = which(predClass.mclust==1)
mclust.g2 = which(predClass.mclust==2)
mclust.g3 = which(predClass.mclust==3)

diff1 = setdiff(mcme.g1,mclust.g1) # magenta3
diff2 = setdiff(mcme.g2,mclust.g2) # cyan3
diff3 = setdiff(mcme.g3,mclust.g3) # orange

g1m = obs[predClass.mclust==1,]
g2m = obs[predClass.mclust==2,]
g3m = obs[predClass.mclust==3,]

ymin1m = min(g1m)
ymax1m = max(g1m)
ymin2m = min(g2m)
ymax2m = max(g2m)
ymin3m = min(g3m)
ymax3m = max(g3m)

ymin = min(c(ymin1,ymin2,ymin3,ymin1m,ymin2m,ymin3m))
ymax = max(c(ymax1,ymax2,ymax3,ymax1m,ymax2m,ymax3m))


# Function for drawing line plots
bwlineplot = function(y,groupnum,method="MCLUST-ME"){
  x = 1:5
  d = rnorm(1000,0,0.01)
  eps1 = 0.01
  eps2 = 0.02
  title = paste(method," Group ",groupnum,sep="")

  # Individual lines with error bars
  plot(x,y[1,],xaxt="n",type="l",ylim=c(ymin,ymax),col="grey78",xlab="",
    ylab="Log2 fold change",main=title,lwd=0.5)
  sd1 = sqrt(diag(errary[,,1]))
  segments(x+d[1],y[1,]-sd1,x+d[1],y[1,]+sd1,lwd=0.5,col="grey78") # add 1 sd error bar
  segments(x+d[1]-eps1,y[1,]-sd1,x+d[1]+eps1,y[1,]-sd1,lwd=0.5,col="grey78")
  segments(x+d[1]-eps1,y[1,]+sd1,x+d[1]+eps1,y[1,]+sd1,lwd=0.5,col="grey78")
  axis(1,at=1:5,labels=c("10min","1h","3h","6h","12h"))
  for(i in 2:nrow(y)){
    lines(x,y[i,],col="grey78",lwd=0.5)
    sd1 = sqrt(diag(errary[,,i]))
    segments(x+d[i],y[i,]-sd1,x+d[i],y[i,]+sd1,lwd=0.5,col="grey78") # add 1 sd error bar
    segments(x+d[i]-eps1,y[i,]-sd1,x+d[i]+eps1,y[i,]-sd1,lwd=0.5,col="grey78")
    segments(x+d[i]-eps1,y[i,]+sd1,x+d[i]+eps1,y[i,]+sd1,lwd=0.5,col="grey78")
  }

  # group mean with error bar
  lines(x,mu.mcme[,groupnum],lwd=2) 
  sd1.mcme = sqrt(diag(sigma.mcme[,,groupnum]))
  mu1.mcme = mu.mcme[,groupnum]
  segments(x,mu1.mcme-sd1.mcme,x,mu1.mcme+sd1.mcme,lwd=2) # add 1 sd error bar
  segments(x-eps2,mu1.mcme-sd1.mcme,x+eps2,mu1.mcme-sd1.mcme,lwd=2)
  segments(x-eps2,mu1.mcme+sd1.mcme,x+eps2,mu1.mcme+sd1.mcme,lwd=2)
}



# Black and white line plots (for paper, differently grouped genes not shown)
par(mfrow=c(2,3))

bwlineplot(g1,1) # MCLUST-ME group 1
bwlineplot(g2,2) # MCLUST-ME group 2
bwlineplot(g3,3) # MCLUST-ME group 3
bwlineplot(g1m,1,"MCLUST") # MCLUST-ME group 1
bwlineplot(g2m,2,"MCLUST") # MCLUST-ME group 2
bwlineplot(g3m,3,"MCLUST") # MCLUST-ME group 3
  
par(mfrow=c(1,1))



# Colored line plots
par(mfrow=c(2,3))
plot(x,g1[1,],xaxt="n",type="l",ylim=c(ymin,ymax),col="magenta3",xlab="",
  ylab="Log2 fold change",main="MCLUST-ME Group 1",lwd=0.5)
axis(1,at=1:5,labels=c("10min","1h","3h","6h","12h"))
for(i in 2:nrow(g1)){
    lines(x,g1[i,],col="magenta3",lwd=0.5)
}

plot(x,g2[1,],xaxt="n",type="l",ylim=c(ymin,ymax),col="cyan3",xlab="",
  ylab="Log2 fold change",main="MCLUST-ME Group 2",lwd=0.5)
axis(1,at=1:5,labels=c("10min","1h","3h","6h","12h"))
for(i in 2:nrow(g2)){
    lines(x,g2[i,],col="cyan3",lwd=0.5)
}

plot(x,g3[1,],xaxt="n",type="l",ylim=c(ymin,ymax),col="orange",xlab="",
  ylab="Log2 fold change",main="MCLUST-ME Group 3",lwd=0.5)
axis(1,at=1:5,labels=c("10min","1h","3h","6h","12h"))
for(i in 2:nrow(g3)){
    lines(x,g3[i,],col="orange",lwd=0.5)
}

plot(x,g1m[1,],xaxt="n",type="l",ylim=c(ymin,ymax),col="grey88",xlab="",
  ylab="Log2 fold change",main="MCLUST Group 1",lwd=0.5)
axis(1,at=1:5,labels=c("10min","1h","3h","6h","12h"))
for(i in 2:nrow(g1m)){
    lines(x,g1m[i,],col="grey88",lwd=0.5)
}
for(i in intersect(diff2,mclust.g1)){
  lines(x,obs[i,],col="cyan3",lwd=1)
}
for(i in intersect(diff3,mclust.g1)){
  lines(x,obs[i,],col="orange",lwd=1)
}

plot(x,g2m[1,],xaxt="n",type="l",ylim=c(ymin,ymax),col="grey88",xlab="",
  ylab="Log2 fold change",main="MCLUST Group 2",lwd=0.5)
axis(1,at=1:5,labels=c("10min","1h","3h","6h","12h"))
for(i in 2:nrow(g2m)){
    lines(x,g2m[i,],col="grey88",lwd=0.5)
}
for(i in intersect(diff1,mclust.g2)){
  lines(x,obs[i,],col="magenta3",lwd=1)
}
for(i in intersect(diff3,mclust.g2)){
  lines(x,obs[i,],col="orange",lwd=1)
}

plot(x,g3m[1,],xaxt="n",type="l",ylim=c(ymin,ymax),col="grey88",xlab="",
  ylab="Log2 fold change",main="MCLUST Group 3",lwd=0.5)
axis(1,at=1:5,labels=c("10min","1h","3h","6h","12h"))
for(i in 2:nrow(g3m)){
    lines(x,g3m[i,],col="grey88",lwd=0.5)
}
for(i in intersect(diff1,mclust.g3)){
  lines(x,obs[i,],col="magenta3",lwd=1)
}
for(i in intersect(diff2,mclust.g3)){
  lines(x,obs[i,],col="cyan3",lwd=1)
}
par(mfrow=c(1,1))




### Look at differently clustered genes
diff1 = setdiff(mcme.g1,mclust.g1) 
diff2 = setdiff(mcme.g2,mclust.g2)
diff3 = setdiff(mcme.g3,mclust.g3) 
diff = c(diff1,diff2,diff3)

# Distribution of order of magnitude of error covariance elements
library(gdata)
od = numeric() # off-diagonal elements
for(i in 1:1000){
  od = c(od,floor(log10(abs(lowerTriangle(errary[,,i])))))
}
dg = numeric() # diagonal elements
for(i in 1:1000){
  dg = c(dg,floor(log10(diag(errary[,,i]))))
}
boxplot(od,dg,names=c("off-diagonal","diagonal"),main="Order of magnitude of error covariance elements")

# Use "total sample variance" from Johnson and Wichern (2007)
library(psych)
temp = numeric(length(diff))
for(i in 1:length(diff)){
  temp[i] = tr(errary[,,diff[i]])
}

temp1 = numeric(1000)
for(i in 1:1000){
  temp1[i] = tr(errary[,,i])
}

qq = numeric(length(temp))
for(i in 1:length(temp)){
  qq[i] = mean(temp1<temp[i])
}
a = data.frame(index=diff,quantile=qq)
attach(a)
b = a[order(quantile),];b

mean(quantile>0.75) # percentage of genes with large errors (among differently clustered ones)
mean(quantile<0.25) # percentage of genes with small errors (among differently clustered ones)
detach(a)

hist(temp1)
hist(temp)







### Uncertainty VS total error variance plot
unc = res.mcme.full$uncertainty
unc.mclust = res.mclust$uncertainty
tv = data.frame(id=1:1000,totalvar=temp1,unc=unc,unc.mclust=unc.mclust)
tv.sort = tv[order(tv$totalvar),]
vert = numeric(length(diff))
for(i in 1:length(diff)){
  vert[i] = which(tv.sort$id==diff[i])
}
plot(1:1000,tv.sort$totalvar,type="l",xlab="",ylab="Total error variance")
for(i in 1:length(vert)){
  lines(c(vert[i],vert[i]),c(0,tv.sort$unc[vert[i]]))
}
abline(h=0)

par(mfrow=c(1,2))
plot(tv.sort$totalvar,tv.sort$unc,cex=0.7,xlab="total error variance",
  ylab="uncertainty",main="MCLUST-ME")
plot(tv.sort$totalvar,tv.sort$unc.mclust,cex=0.7,xlab="total error variance",
  ylab="uncertainty",main="MCLUST")
par(mfrow=c(1,1))

mcme.sig = res.mcme.full$parameters$sigmahat
mean.error = apply(errary,c(1,2),mean)
mevvv.sig = res.mclust$parameters$variance$sigma

mcme.sig[,,1]+mean.error
mevvv.sig[,,1]

mcme.sig[,,2]+mean.error
mevvv.sig[,,2]

mcme.sig[,,3]+mean.error
mevvv.sig[,,3]

unc.ratio = (unc+0.005)/(unc.mclust+0.005)
plot(tv.sort$totalvar,(tv.sort$unc+0.005)/(tv.sort$unc.mclust+0.005),cex=0.5,log="y")
abline(h=1)
hist(log(unc.ratio))










