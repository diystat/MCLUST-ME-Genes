#######################################################################
###------------------------ RNA-Seq Example ------------------------###
#######################################################################

setwd("/home/yanming/ongoing/Model-based-Clustering-Research");

library(mclust);
library(phyclust)
library(psych)

## library(MASS)
## Load required packages
## library(gdata)
## library(caret)
## library(openintro)

## Load raw data
load("Data/rna_raw.RData")

## Load sampled data and the mclust results
m = 1000;
seed = 999;
## seed = 1001;
file.mclust = paste(sprintf("Results/res_real2d_v2_m%d_seed%d_mclust.Rdata", m, seed));
list(load(file.mclust));

## res.mclust.4 = Mclust(obs, G = 4, modelNames="VVV");
## plot(res.mclust.4);

plot(res.mclust, "classification");

## Load MCME results
BIC = numeric(8);
for (g in 1:8) {
  fname = paste(sprintf("Results/res_real2d_v2_m%d_seed%d_mcme_g%d.rds",m, seed, g))
  BIC[g] = readRDS(fname)$BIC;
}
plot(res.mclust, "BIC", ylim = c(-3500, -2900));
points(BIC, xlab="Number of components", ylab="BIC", type="b", pch=0);
## plot(BIC, xlab="Number of components", ylab="BIC", type="b", pch=0);

## g = 2, 3, or 4 are typically the best mdoels;
g = 2;
fname = paste(sprintf("Results/res_real2d_v2_m%d_seed%d_mcme_g%d.rds",m, seed, g));
res.mcme = readRDS(fname);

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


file.png = sprintf("Graphs/real_data_seed%d_clustering_results.png", seed);
png(file.png, width=14, height=7, units="in", res=270, pointsize=14)

## Compare the two clustering results

## cols = c("magenta", "darkcyan", "brown", "black");
cols = mclust.options("classPlotColors");
pchs = mclust.options("classPlotSymbols");

par(mfrow=c(1,2))
plot(obs, col = cols[group.mclust], xlab="1h", ylab="3h", main="MCLUST Clustering", pch=pchs[group.mclust], cex=0.8)
abline(h=0,v=0,lty="dashed",col="gray")


plot(obs, col = cols[group.mcme], xlab="1h", ylab="3h", main="MCLUST-ME Clustering", pch=pchs[group.mcme], cex=0.8)
abline(h=0,v=0,lty="dashed",col="gray")

par(mfrow=c(1,1))

dev.off()

## See parameters of the fitted model
res.mclust$parameters

res.mcme$parameters$muhat;
res.mcme$parameters$tauhat;
res.mcme$parameters$sigmahat;


## Summarize error variation

## Use trace()
st = numeric(1000)
for(i in 1:1000){
  st[i] = tr(errary[,,i])
}

## Use det()
sd = numeric(1000)
for(i in 1:1000){
  sd[i] = det(errary[,,i])
}


## Identify points of disgreement
table(group.mclust, group.mcme);

id12 = group.mclust ==1 & group.mcme==2;
id21 = group.mclust ==2 & group.mcme==1;

id  = group.mclust != group.mcme;

## Examine the error distribution of among the points of disagreements  
quantile(s);
## quantile(s[id]);
quantile(s[id12]);
quantile(s[id21]);

file.png = sprintf("Graphs/real_data_seed%d_error_distribution.png", seed);

png(file.png, width=7, height=7, units="in", res=270, pointsize=14)

boxplot(st, st[id12], st[id21], names = 1:3, log="y", main="Distribution of the traces of the error covariance matrices"); 

dev.off();


## Examien the z-values
res.mcme$z[id12,];
res.mclust$z[id12,];

data.frame(z.mclust = res.mcme$z[id12,], z.mcme=res.mclust$z[id12,]);
data.frame(z.mclust = res.mcme$z[id21,], z.mcme=res.mclust$z[id21,]);


## Extract and compare membership probabilities
g1.mcme = res.mcme$z[,1]
g1.mclust = res.mclust$z[,1]

g2.mcme = res.mcme$z[,2]
g2.mclust = res.mclust$z[,2]

## g3.mcme = res.mcme$z[,3]
## g3.mclust = res.mclust$z[,3]

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

plot(unc.mclust, unc.mcme)

abline(0,1,lty="dashed")

plot(log(s), unc.mcme)
plot(log(s), unc.mclust)


plot(log(s), unc.mcme-unc.mclust, col=ifelse(s>mean(s),"red","dodgerblue"), 
  main="Pairwise uncertainty difference",xlab="generalized error variance",
  ylab="MCLUST-ME--MCLUST")
abline(h=0, v=mean(s), lty="dashed")
legend("topleft", legend=c(">mean","<=mean"),col=c("red","blue"),pch=1)

