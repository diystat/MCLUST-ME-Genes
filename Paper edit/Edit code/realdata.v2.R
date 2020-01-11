#######################################################################
###------------------------ RNA-Seq Example ------------------------###
#######################################################################

setwd("/home/yanming/ongoing/Model-based-Clustering-Research");

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

## Random sample 1000 genes, excluding rows with 0 relative frequnencies under any treatment
n.genes = nrow(nb.data$counts);
zero = apply(nb.data$rel.frequencies < 1/sum(nb.data$counts), 1, any);

m = 1000;

## seed = 999; set.seed(seed);

seed = 1001; set.seed(seed);

ss = sample((1:n.genes)[!zero], m);

if (FALSE) {
  par(mfrow=c(1,2));
  obs = full$beta[top, 2:3] # data points
  errary = v.beta[2:3,2:3,top] # measurement error covariances
  plot(obs);

}

obs = full$beta[ss, 2:3] # data points
errary = v.beta[2:3,2:3,ss] # measurement error covariances
plot(obs);

## Run mevvv:
res.mclust = Mclust(obs, modelNames="VVV")
plot(res.mclust, what="classification")
plot(res.mclust, what="BIC");

## Save the sampled data
file.mclust = paste(sprintf("Results/res_real2d_v2_m%d_seed%d_mclust.Rdata", m, seed));
save(m, seed, ss, obs, errary, res.mclust, file=file.mclust); 


if (FALSE) {
  ## It seems Mclust is deterministic
  res.mclust.2 = Mclust(obs, modelNames="VVV")
  identical(res.mclust, res.mclust.2)
}

## saveRDS(res.mclust, "Results/res_real2d_mclust.rds")
## res.mclust = readRDS("Results/res_real2d_mclust.rds")

#### Mclust gives 3 clusters (if we cluster the top 1000 DE genes)

#### 10/23/2019 Depending on the seeds used, Mclust gives 2--4
#### clusters (one cluster of non-DE, 1 or more clusters of DE genes)


fit.mcme = function(g){
  ## Generate initial group membership for mclust-me
  N = nrow(obs);
  G = g
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

  fname = paste(sprintf("Results/res_real2d_v2_m%d_seed%d_mcme_g%d.rds",m, seed, g))
  saveRDS(res.mcme, fname)
}


for(i in 1:8){
  print(paste("G = ",i,sep=""))
  fit.mcme(i)
}


if (FALSE) {
## Old codes

## saveRDS(res.mcme, "Results/res_real2d_mcme.rds")  
res.mcme = readRDS("Results/res_real2d_v2_m1000_seed999_mcme_g2.rds")

## res.mcme = readRDS("Results/res_real2d_v2_m1000_seed999_mcme_g4.rds")
                           

## Extract cluster membership
group.mclust = res.mclust$classification

group.mcme = numeric(1000)
for(i in 1:1000){
  ## group.mcme[i] = which(res.mcme$member[i,]==1)
  group.mcme[i] = order(-res.mcme$z[i,])[1];
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

## points(obs[id2,], pch = "o")

plot(obs, col = c("orangered","seagreen","dodgerblue")[group.mcme], xlab="1h", ylab="3h",
  main="MCLUST-ME Clustering", pch=c(15,16,17)[group.mcme], cex=0.8)
abline(h=0,v=0,lty="dashed",col="gray")
par(mfrow=c(1,1))

## points(obs[id2,], pch = "o")

dev.off()


## Identify the four points that were blue in mclust, but green/red in mclust-me
id = group.mcme == 3 & group.mclust != 3;

obs[id,];
errary[,,id];

res.mclust$z[id,];
res.mclust$classification[id];

res.mcme$z[id,];
res.mcme$member[1:1000,][id,];

## Identify the "neighbors" of the above four points
id2 = obs[,1] < 0.5 & obs[,1]>-0.5 & obs[,2]>0 & obs[,2]<1;
points(obs[id2,], pch = "o")
errary[,,id2];

res.mclust$z[id2,];
res.mclust$classification[id2];

res.mcme$z[id2,];
res.mcme$member[1:1000,][id2,];



res.mcme = readRDS("Results/res_real2d_mcme.rds")
res.mclust = readRDS("Results/res_real2d_mclust.rds")

## Extract and compare membership probabilities
e1.mcme = res.mcme$z[,1]
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

}
