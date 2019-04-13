#######################################################################
###------------------------ RNA-Seq Example ------------------------###
#######################################################################

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
library(psych)


###------------------- Importing Data ---------------------###

###################################################################
### Raw and processed data available in "Data" folder.          ###
###################################################################

# Import raw data
load("Data/rna_raw.RData")
obs = full$beta[top, 2:3] # data points
errary = v.beta[2:3,2:3,top] # measurement error covariances

## Run mevvv:
res.mclust = Mclust(obs, modelNames="VVV")
saveRDS(res.mclust, "Results/res_real2d_mclust.rds")
res.mclust = readRDS("Results/res_real2d_mclust.rds")
plot(res.mclust, what="classification")
#### Mclust gives 3 clusters  


fit.mcme = function(g){
  ## Generate initial group membership for mclust-me
  N = nrow(obs); G = g
  hcTree = hc(obs) # hc() is model-based hierarchical clustering
  cl = hclass(hcTree, G)
  z.ini = matrix(0,N,G)
  for(i in 1:N){
    for(j in 1:G){
      z.ini[i,j] = ifelse(cl[i]==j,1,0)
    }
  }
  
  ## Run MCME:
  res.mcme = mcmeVVV(obs, z.ini, errary)

  fname = paste("Results/res_real2d_mcme_g",g,".rds",sep="")
  saveRDS(res.mcme, fname)
}

for(i in 1:8){
  print(paste("G = ",i,sep=""))
  fit.mcme(i)
}





saveRDS(res.mcme, "Results/res_real2d_mcme.rds")  
res.mcme = readRDS("Results/res_real2d_mcme.rds")







## Extract cluster membership
group.mclust = res.mclust$classification
group.mcme = numeric(1000)
for(i in 1:1000){
  group.mcme[i] = which(res.mcme$member[i,]==1)
}

## Similarity between two groupings
RRand(group.mclust, group.mcme) # from package 'phyclust'
#   Rand adjRand  Eindex 
# 0.8070  0.5918  0.2900 

png("scatter_read2d.png", width=14, height=7, units="in", res=270, pointsize=14)
par(mfrow=c(1,2))
plot(obs, col = c("orangered","seagreen","dodgerblue")[group.mclust], xlab="1h", ylab="3h",
  main="MCLUST Clustering", pch=c(15,16,17)[group.mclust], cex=0.8)
abline(h=0,v=0,lty="dashed",col="gray")
plot(obs, col = c("orangered","seagreen","dodgerblue")[group.mcme], xlab="1h", ylab="3h",
  main="MCLUST-ME Clustering", pch=c(15,16,17)[group.mcme], cex=0.8)
abline(h=0,v=0,lty="dashed",col="gray")
par(mfrow=c(1,1))
dev.off()



res.mcme = readRDS("Results/res_real2d_mcme.rds")
res.mclust = readRDS("Results/res_real2d_mclust.rds")

## Extract and compare membership probabilities
g1.mcme = res.mcme$z[,1]
g1.mclust = res.mclust$z[,1]

g2.mcme = res.mcme$z[,2]
g2.mclust = res.mclust$z[,2]

g3.mcme = res.mcme$z[,3]
g3.mclust = res.mclust$z[,3]

par(mfrow=c(1,3))
plot(g1.mcme, g1.mclust, cex=0.7, main="Pr(Cluster 1)", xlab="MCLUST-ME", ylab="MCLUST")
abline(0,1,lty="dashed")
abline(h=0.5,v=0.5,lty="dashed")

plot(g2.mcme, g2.mclust, cex=0.7, main="Pr(Cluster 2)", xlab="MCLUST-ME", ylab="MCLUST")
abline(0,1,lty="dashed")
abline(h=0.5,v=0.5,lty="dashed")

plot(g3.mcme, g3.mclust, cex=0.7, main="Pr(Cluster 3)", xlab="MCLUST-ME", ylab="MCLUST")
abline(0,1,lty="dashed")
abline(h=0.5,v=0.5,lty="dashed")
par(mfrow=c(1,1))






## Extract clustering uncertainty
unc.mcme = res.mcme$uncertainty
unc.mclust = 1-apply(res.mclust$z, 1, max)


## Summarize error variation
s = numeric(1000)
for(i in 1:1000){
  s[i] = tr(errary[,,i])
}


plot(unc.mclust, unc.mcme)
abline(0,1,lty="dashed")



plot(log(s), unc.mcme)
plot(log(s), unc.mclust)


plot(log(s), unc.mcme-unc.mclust, col=ifelse(s>mean(s),"red","dodgerblue"), 
  main="Pairwise uncertainty difference",xlab="generalized error variance",
  ylab="MCLUST-ME--MCLUST")
abline(h=0, v=mean(s), lty="dashed")
legend("topleft", legend=c(">mean","<=mean"),col=c("red","blue"),pch=1)

