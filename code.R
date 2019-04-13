

## Set working directory
setwd("/Users/wzhang/Project 1 Code")

## Source all necessary functions
source("Code/core_functions.R")
source("Code/simulation_functions.R")

## Load required packages
library(MASS)
library(gdata)
library(caret)
library(mclust)
library(phyclust)
library(openintro)


## Import raw data
load("Data/rna_raw.RData")
obs = full$beta[top, 1:5] # data points
errary = v.beta[,,top] # measurement error covariances

## MCLUST clustering results
res.mclust = Mclust(obs,modelNames="VVV") # Running time < 5min
summary(res.mclust)

## MCLUST-ME clustering results
load("Results/rna_3group.RData")
res.mcme = run3



###---------------------- Analysis -----------------------###

### Obtain mean and covariance estimates
# MCLUST-ME
param = res.mcme$parameters
mu.mcme = param$muhat # centers
sigma.mcme = param$sigmahat # covariances

# MCLUST
param = res.mclust$parameters
mu.mclust = param$mean # centers
sigma.mclust = param$variance$sigma # covariances


### Predicted class labels:
# MCLUST-ME labels
zhat = res.mcme.full$z
predClass = numeric(1000)
for(i in 1:1000){
  temp = zhat[i,]
  predClass[i] = which(temp==max(temp))
}

# MCLUST labels
zhat.mclust = res.mclust$z
predClass.mclust = numeric(1000)
for(i in 1:1000){
  temp = zhat.mclust[i,]
  predClass.mclust[i] = which(temp==max(temp))
}


### Confusion matrix:
confusionMatrix(predClass,predClass.mclust,dnn=c("MCME","MCLUST"))
conf = cbind(predClass.mclust,predClass)
m = data.frame(conf)
contTable(m)




## Which points are clustered differently by the two methods?
diff_ind = which(predClass!=predClass.mclust)

## Clustering uncertainties for two methods
unc_diff_mcme = res.mcme$uncertainty[diff_ind]
unc_diff_mclust = res.mclust$uncertainty[diff_ind]

## Points clustered the same by two methods
same_ind = which(predClass==predClass.mclust)
same_obs = obs[same_ind,]
diff_obs = obs[diff_ind,]






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

#png("rna.png",width=14,height=11,units="in",res=300)
tiff("all_unc.tiff",width=1680,height=980,compression="lzw",res=240,pointsize=9)
par(mfrow=c(1,2),cex.lab=1.2,cex.axis=1.2,cex.main=1.2)

plot(diff_obs[,2],diff_obs[,3],main="MCLUST-ME Clusters",pch=c(1,2,3)[predClass[diff_ind]],cex=1.2,
  xlab="1h",ylab="3h",col="#00000070")
# Plot the confidence outlines
lines(g1.coord,type="l",lwd=1)
lines(g2.coord,type="l",lwd=1)
lines(g3.coord,type="l",lwd=1)

plot(diff_obs[,2],diff_obs[,3],main="MCLUST Clusters",pch=c(1,2,3)[predClass.mclust[diff_ind]],cex=1.2,
  xlab="1h",ylab="3h",col="#00000070")
# Plot the confidence outlines
lines(g1.coord.mc,type="l",lwd=1)
lines(g2.coord.mc,type="l",lwd=1)
lines(g3.coord.mc,type="l",lwd=1)

par(mfrow=c(1,1),cex.lab=1,cex.axis=1,cex.main=1)
dev.off()

plot(res.mclust, what = "uncertainty",  main = "")
plot(res.mcme, what="uncertainty", main="")


    


    
#### Plot clustering uncertainty of differently clustered points

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
 
       
colors = c("firebrick3", "dodgerblue3", "chartreuse3")
unc_mcme = res.mcme$uncertainty
unc_mclust = res.mclust$uncertainty

u_mcme = (unc_mcme - min(unc_mcme))/(max(unc_mcme) - 
            min(unc_mcme))
u_mclust = (unc_mclust - min(unc_mclust))/(max(unc_mclust) - 
            min(unc_mclust))


tiff("diff_unc.tiff",width=1680,height=980,compression="lzw",res=240,pointsize=9)
par(mfrow=c(1,2),cex.lab=1.2,cex.axis=1.2,cex.main=1.2)

plot(obs[diff_ind, 2], obs[diff_ind, 3], pch = 21, col = alpha(colors[predClass[diff_ind]], 0.3), 
    bg = alpha(colors[predClass[diff_ind]], 0.3), lwd = 1.5, cex = 2*u_mcme,
  xlab = "1h", ylab = "3h", main = "Clustering Uncertainty: MCLUST-ME")
lines(g1.coord, type="l", lwd=1, lty="dotted")
lines(g2.coord, type="l", lwd=1, lty="dotted")
lines(g3.coord, type="l", lwd=1, lty="dotted")  

plot(obs[diff_ind, 2], obs[diff_ind, 3], pch = 21, col = alpha(colors[predClass.mclust[diff_ind]], 0.3), 
    bg = alpha(colors[predClass.mclust[diff_ind]], 0.3), lwd = 1.5, cex = 2*u_mclust,
  xlab = "1h", ylab = "3h", main = "Clustering Uncertainty: MCLUST")  
lines(g1.coord.mc, type="l", lwd=1, lty="dotted")
lines(g2.coord.mc, type="l", lwd=1, lty="dotted")
lines(g3.coord.mc, type="l", lwd=1, lty="dotted")  
  
par(mfrow=c(1,1),cex.lab=1,cex.axis=1,cex.main=1)
dev.off()
  


  
tiff("all_unc.tiff",width=1680,height=980,compression="lzw",res=240,pointsize=9)
par(mfrow=c(1,2),cex.lab=1.2,cex.axis=1.2,cex.main=1.2)

plot(obs[, 2], obs[, 3], pch = 21, col = alpha(colors[predClass], 0.3), 
    bg = alpha(colors[predClass], 0.3), lwd = 1.5, cex = 2*u_mcme,
  xlab = "1h", ylab = "3h", main = "Clustering Uncertainty: MCLUST-ME")
lines(g1.coord, type="l", lwd=1, lty="dotted")
lines(g2.coord, type="l", lwd=1, lty="dotted")
lines(g3.coord, type="l", lwd=1, lty="dotted")  

plot(obs[, 2], obs[, 3], pch = 21, col = alpha(colors[predClass.mclust], 0.3), 
    bg = alpha(colors[predClass.mclust], 0.3), lwd = 1.5, cex = 2*u_mclust,
  xlab = "1h", ylab = "3h", main = "Clustering Uncertainty: MCLUST")  
lines(g1.coord.mc, type="l", lwd=1, lty="dotted")
lines(g2.coord.mc, type="l", lwd=1, lty="dotted")
lines(g3.coord.mc, type="l", lwd=1, lty="dotted")  
  
par(mfrow=c(1,1),cex.lab=1,cex.axis=1,cex.main=1)
dev.off()  
  
  



