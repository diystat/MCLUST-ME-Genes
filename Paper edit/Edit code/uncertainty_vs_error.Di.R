## Linux working directory
## setwd("/Users/wzhang/Project 1 Code")
## setwd("/home/yanming/ongoing/Model-based-Clustering-Research");

## Windows work directory
setwd("D:/ongoing/Model-based-Clustering-Research");

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


run.simulation = function(){


  ## Set mean and covariance parameters
  seed = 999; set.seed(seed);
  N = 200
  tau = 0.5
  mu1 = c(-6,0)
  mu2 = c(6,0)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(64,0,0,64),nrow=2)
  s = runif(200, 0, 36)


  ## Set mean and covariance parameters
  seed = 1001; set.seed(seed);

  N = 200
  tau = 0.5
  mu1 = c(-10,0)
  mu2 = c(10,0)
  sig1 = matrix(c(100,0,0,100),nrow=2)
  sig2 = matrix(c(100,0,0,100),nrow=2)
  s = runif(200, 0, 100)

  ## Set mean and covariance parameters
  seed = 1001; set.seed(seed);

  N = 1000
  tau = 0.5
  mu1 = c(-10,0)
  mu2 = c(10,0)
  sig1 = matrix(c(100,0,0,100),nrow=2)
  sig2 = matrix(c(100,0,0,100),nrow=2)
  s = runif(N, 0, 100)


  ## Generate uniformly distributed error variance values
  E = matrix(c(1,0,0,1),nrow=2)
  errmat = array(0,dim=c(2,2,N))
  for(i in 1:N){
    errmat[,,i] = s[i]*E
  }

  ## Generate sample from mixture distribution:
  U = runif(N)
  z = 1 * (U< tau) + 2 * (U >=tau);

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

  ## Save simulation parameters and simulated data
  file.sim.pars.data = sprintf("Paper edit/Edit results/sim.pars.data.seed%d.Rdata", seed);
  save(seed, N, tau, mu1, mu2, sig1, sig2, z, s, errmat, rand.samples,file=file.sim.pars.data); 
  

  ## Run MCME:
  res.mcme = mcmeVVV(rand.samples, z.ini, errmat);
  file.res.mcme = sprintf("Paper edit/Edit results/mcme_res.seed%d.rds", seed);
  saveRDS(res.mcme, file.res.mcme);

  ## saveRDS(res.mcme, "Paper edit/Edit results/mcme_res.seed.rds");
  
  ## Run mevvv:
  res.mevvv = meVVV(rand.samples, z.ini)
  file.res.mevvv = sprintf("Paper edit/Edit results/mevvv_res.seed%d.rds", seed);
  saveRDS(res.mevvv, file.res.mevvv);

  ## saveRDS(res.mevvv, "Paper edit/Edit results/mevvv_res.rds")
}

analyze.results = function() {
  library(latex2exp)

  ## mcme: mclust-me
  ## mevvv: mclust
  
  ## res.mcme = readRDS("Paper edit/Edit results/mcme_res.rds")
  ## res.mevvv = readRDS("Paper edit/Edit results/mevvv_res.rds")

  ## seed = 999;

  ## seed = 1001 results were presented in the GENES paper
  seed = 1001;
  
  file.sim.pars.data = sprintf("Paper edit/Edit results/sim.pars.data.seed%d.Rdata", seed);
  load(file.sim.pars.data);

  file.res.mcme = sprintf("Paper edit/Edit results/mcme_res.seed%d.rds", seed);
  res.mcme = readRDS(file.res.mcme);

  file.res.mevvv = sprintf("Paper edit/Edit results/mevvv_res.seed%d.rds", seed);
  res.mevvv = readRDS(file.res.mevvv);
  

  ## Extract cluster memberships
  N = res.mcme$n;

  cl.mcme = cl.mevvv = numeric(N)
  
  for(i in 1:N){
    cl.mcme[i] = which(res.mcme$z[i,]==max(res.mcme$z[i,]))
    cl.mevvv[i] = which(res.mevvv$z[i,]==max(res.mevvv$z[i,]))
  }
  
  all.equal(cl.mcme, (res.mcme$z[,1]<0.5) + 1);
  all.equal(cl.mevvv, (res.mevvv$z[,1]<0.5) + 1);

  table(cl.mcme, cl.mevvv);
  id.diff = cl.mcme!=cl.mevvv;
  sum(id.diff)

  table(cl.mcme, z);
  table(cl.mevvv, z);

  
  ## Figure 6 of GENES paper
  ## Last updated: 2020-Feb-03
  file.png = sprintf("Paper edit/Edit graphs/clustering_results_%d.png", seed);
  png(file.png, width=1200, height=600, units="px");

  ## Plot classifications
  par(mfrow=c(1,2))
  plot(rand.samples, col=c("magenta","darkcyan")[cl.mevvv], pch=c(0,1)[cl.mevvv], lwd=2,
       xlab ="Y1", ylab="Y2", main="MCLUST", cex=2, cex.axis=1.5, cex.lab=1.2)
  points(rand.samples[z != cl.mevvv,], pch=4, lwd=2);
  points(rand.samples[id.diff,], pch=1, cex=4,lwd=2, col="black")
  
  plot(rand.samples, col=c("magenta","darkcyan")[cl.mcme], pch=c(0,1)[cl.mcme], lwd=2,
       xlab ="Y1", ylab="Y2", main="MCLUST-ME", cex=2, cex.axis=1.5, cex.lab=1.2)
  points(rand.samples[z != cl.mcme,], pch=4, lwd=2);
##  points(rand.samples[z==1,], col="red", pch=2);
##  points(rand.samples[z==2,], col="blue", pch=2);
  points(rand.samples[id.diff,], pch=1, cex=4,lwd=2, col="black")
  
  par(mfrow=c(1,1))
  dev.off();


  ## Estimated model parameters
  str(res.mevvv);
  str(res.mcme);

  res.mevvv$parameters$mean
  res.mevvv$parameters$pro
  res.mevvv$parameters$variance$sigma

  res.mcme$parameters$muhat;
  res.mcme$parameters$tauhat;
  res.mcme$parameters$sigmahat;
  
  ## Extract and compare membership probabilities
  g1.mcme = res.mcme$z[,1]
  g1.mevvv = res.mevvv$z[,1]

  plot(g1.mevvv, g1.mcme)
  abline(0,1,lty="dashed")
  abline(h=0.5, lty=2);
  abline(v=0.5, lty=2);

  ## Extract and plot clustering uncertainties
  unc.mcme = res.mcme$uncertainty
  all.equal(apply(res.mcme$z, 1, min), res.mcme$uncertainty);
  
  unc.mevvv = apply(res.mevvv$z, 1, min)
  
  delta = res.mcme$z[,1] - res.mevvv$z[,1];
  
  ## Figure 7 of Genes paper
  ## Last updated: 2020-Feb-03
  file.png = sprintf("Paper edit/Edit graphs/change_prob1_%d.png", seed);
  ## png(file.png, width=1200, height=600, units="px");
  png(file.png, width=14, height=10, units="in", res=270, pointsize=14)

  ## Use arrows to show estimated membership probability to cluster 1.
  par(mar=c(5.1, 4.1, 4.1, 4.1), xpd=TRUE)
  plot(s, res.mevvv$z[,1], 
       xlab="magnitude of error variance", 
       ylab="membership probability to cluster 1", 
       pch = 20, 
       cex.lab = 1.4,
       main=TeX("Changes in membership probability to cluster 1: MCLUST $\\rightarrow$ MCLUST-ME"));
  
  ## points(s, res.mevvv$z[,1]);
  for(i in 1:N){
    if (abs(delta[i])>0.01) {
      arrows(s[i], res.mevvv$z[i,1], s[i], res.mcme$z[i,1], 
             angle = 15, lwd=1.3, length=0.1,
             col=ifelse(unc.mevvv[i]>unc.mcme[i],"red","dodgerblue"))
    }
  }
  ## legend("topright", legend=c("decrease", "increase"), col=c("red","dodgerblue"),lty=1, inset=c(-0.13,0), cex=0.7)
  segments(-3, 0.5, 103, 0.5, lty=2);
  dev.off();
  
 
  ## Figure 8. Points with most changes in estimated membership probabilities
  ## Last updated: Feb 3, 2020
  id = abs(delta) > 0.1;
  cols = ifelse(unc.mevvv > unc.mcme, "red", "dodgerblue");
  
  file.png = sprintf("Paper edit/Edit graphs/points_most_changes_in_uncertainty_%d.png", seed);
  ## png(file.png, width=1200, height=600, units="px");
  png(file.png, width=10, height=10, units="in", res=300, pointsize=16);
  
  plot(rand.samples, xlab="Y1", ylab="Y2", cex=2, 
       main="Points with most changes in estimated membership probabilities\n between MCLUST and MCLUST-ME");
  points(rand.samples[id,], col = cols[id], pch=19, cex=2)
  dev.off();
  

  format(data.frame(s=s, unc.mevvv=unc.mevvv, unc.mcme=unc.mcme), digits=2);
  table(unc.mcme > unc.mevvv);

  old.pars = par(mfrow=c(1,2));
  plot(s, unc.mcme)
  plot(s, unc.mevvv)
  par(old.pars);

  plot(unc.mevvv, unc.mcme-unc.mevvv)

  plot(unc.mevvv, unc.mcme) # Use point size to show measurement error magnitude
  abline(0,1,lty="dashed")


  ## Plot pairwise difference in uncertainty against error variance
  plot(s, unc.mcme-unc.mevvv, col=ifelse(s>mean(s),"red","dodgerblue"), 
       main="Pairwise uncertainty difference",xlab="generalized error variance",
       ylab="MCLUST-ME--MCLUST")
  abline(h=0, v=mean(s), lty="dashed")
  # legend("topleft", legend=c(">mean","<=mean"),col=c("red","blue"),pch=1)


  ## Use arrows to show change in uncertainty between two methods
  file.png = sprintf("Paper edit/Edit graphs/unc_change_%d.png", seed);
  png(file.png, width=1200, height=600, units="px");
  
  ## png("unc_change_2.png", width=14, height=10, units="in", res=270, pointsize=14)

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



