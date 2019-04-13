
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


unc.diff = function(){

## Set mean and covariance parameters
N = 20
tau = 0.5
mu1 = c(-6,0)
mu2 = c(6,0)
sig1 = matrix(c(64,0,0,64),nrow=2)
sig2 = matrix(c(64,0,0,64),nrow=2)


## Generate uniformly distributed error variance values
s = runif(N, 0, 36)
E = matrix(c(1,0,0,1),nrow=2)
errmat = array(0,dim=c(2,2,N))
for(i in 1:N){
  errmat[,,i] = s[i]*E
}

## Generate sample from mixture distribution:
U = runif(N)
z.ini = matrix(0,N,2)
rand.samples = matrix(0,N,2)
for(i in 1:N){
  if(U[i]<tau){
    rand.samples[i,] = MASS::mvrnorm(1,mu1,(sig1+errmat[,,i]))
    z.ini[i,] = c(1,0)
  } else{
    rand.samples[i,] = MASS::mvrnorm(1,mu2,(sig2+errmat[,,i]))
    z.ini[i,] = c(0,1)
  }
}

## Run MCME:
res.mcme = mcmeVVV(rand.samples, z.ini, errmat)
saveRDS(res.mcme, "Paper edit/Edit results/mcme_res_small.rds")
  
## Run mevvv:
res.mevvv = meVVV(rand.samples,z.ini)
saveRDS(res.mevvv, "Paper edit/Edit results/mevvv_res_small.rds")

## Extract cluster memberships
cl.mcme = cl.mevvv = numeric(N)
for(i in 1:N){
  cl.mcme[i] = which(res.mcme$z[i,]==max(res.mcme$z[i,]))
  cl.mevvv[i] = which(res.mevvv$z[i,]==max(res.mevvv$z[i,]))
}

## Plot classifications
plot(rand.samples, col=c("red","dodgerblue")[cl.mcme], pch=c(15,16)[cl.mcme])
plot(rand.samples, col=c("red","dodgerblue")[cl.mevvv], pch=c(15,16)[cl.mcme])

  
## Extract and plot clustering uncertainties
unc.mcme = res.mcme$uncertainty
unc.mevvv = apply(res.mevvv$z, 1, min)

plot(s, unc.mcme)
plot(s, unc.mevvv)

plot(unc.mevvv, unc.mcme-unc.mevvv)

plot(unc.mevvv, unc.mcme) # Use point size to show measurement error magnitude
abline(0,1,lty="dashed")


## Plot pairwise difference in uncertainty against error variance
plot(s, unc.mcme-unc.mevvv, col=ifelse(s>mean(s),"red","dodgerblue"), 
  main="Pairwise uncertainty difference",xlab="generalized error variance",
  ylab="MCLUST-ME--MCLUST")
abline(h=0, v=mean(s), lty="dashed")
#legend("topleft", legend=c(">mean","<=mean"),col=c("red","blue"),pch=1)



## Use arrows to show change in uncertainty between two methods
png("unc_change_2.png", width=14, height=10, units="in", res=270, pointsize=14)
par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
plot(s, unc.mcme, xlab="generalized error variance", ylab="membership uncertainty", 
  main="Uncertainty change: MCLUST-->MCLUST-ME", type="n")
points(s, unc.mevvv, type="n")
for(i in 1:N){
  arrows(s[i], unc.mevvv[i], s[i], unc.mcme[i], lwd=1.3, length=0.1,
    col=ifelse(unc.mevvv[i]>unc.mcme[i],"red","dodgerblue"))
}
legend("topright", legend=c("decrease", "increase"),
  col=c("red","dodgerblue"),lty=1, inset=c(-0.13,0), cex=0.7)
dev.off()

}



