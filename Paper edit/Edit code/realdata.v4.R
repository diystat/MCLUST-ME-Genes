#######################################################################
###------------------------ RNA-Seq Example ------------------------###
#######################################################################

## 01/25/2020
## Try the mError method by Kumar and Patel on the real data set.

##  Linux working directory
setwd("/home/yanming/ongoing/Model-based-Clustering-Research");

## Windows working directory
## setwd("D:/ongoing/Model-Based-Clustering-Research");

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

library(NBPSeq);


###------------------- Importing Data ---------------------###

###################################################################
### Raw and processed data available in "Data" folder.          ###
###################################################################

# Import raw data
print(load("Data/rna_raw.RData"));

date = "2020-01-15";
m = 1000;
seed = 999;
file.results = sprintf("Results/res_real2d_v3_m%d_seed%d.%s.Rdata", m, seed, date);
print(load(file.results));


##
plot(obs);

## Should be 775 & 10 // 30 & 185
table(res.mclust$classification, res0$classification);

compute.iv = function(v) {
  m = dim(v)[3];
  iv = v;
  for (i in 1:m) {
    iv[,,i] = ginv(v[,,i]);
  }
  iv;
}

test.compute.iv = function(){
  iv = compute.iv(errary);

  for (i in 1:m) {
    print(errary[,,i] %*% iv[,,i]);
  }

}



compute.center = function(y, iv) {
  m = nrow(y);
  d = ncol(y);

  sivy = iv[,,1]%*% y[1,];
  siv = iv[,,1];


  if (m > 1) {
    for (i in 2:m) {
      sivy = sivy + iv[,,i]%*%y[i,];
      siv = siv + iv[,,i];
    };
  }

  ginv(siv) %*% sivy;

}

test.compute.center() {
  y = obs;

  ## use identity covariances
  iv = errary;
  iv[1,1,] = 1;
  iv[2,2,] = 1;
  iv[1,2,] = 0;
  iv[2,1,] = 0;
  
  compute.center(y, iv);
  colMeans(y);

}

## Compute distance matrix using the mError formula
## @param y observations, m by d.
## @param v covariances for each y, d x d x m.
## @mu the centers
compute.distance.matrix = function(y, iv, mu) {
  m = nrow(y);
  d = ncol(y);
  G = nrow(mu);

  dist = matrix(NA, m, G);

  for (k in 1:G) {
    for (i in 1:m) {
      ymu = matrix(y[i,] -  mu[k,], d, 1);
      dist[i,k] = t(ymu) %*% iv[,,i] %*% ymu;
    }
  }

  dist
}

test.compute.distance.matrix() {
  y = obs;
  iv = compute.iv(errary);
  mu = res0$parameters$muhat;
  
  dist = compute.distance.matrix(y, iv, mu);

  cl = numeric(m);
  cl[] = 2;
  cl[dist[,2]>dist[,1]] = 1;

  plot(obs, col=cl);

  table(cl, res0$classification);

  ## Use identity covariance matrices
  iv = errary;
  iv[1,1,] = 1;
  iv[2,2,] = 1;
  iv[1,2,] = 0;
  iv[2,1,] = 0;

  dist = compute.distance.matrix(y, iv, mu);

  y1 = t(t(y) - mu[1,]);
  dim(y1);
  d1 = rowSums(y1^2);
  all.equal(dist[,1], d1);
    
  y2 = t(t(y) - mu[2,]);
  dim(y2);
  d2 = rowSums(y2^2);
  all.equal(dist[,2], d2);

}

## One step of mError: find the centers of the clusters
##
## @param y observations, m by d.
## @param iv inverses of covariances for each y, d x d x m.
## @param cl current classifications
##
## output: new centers
mError.1 = function(y, iv, cl) {
  m = nrow(y);
  d = ncol(y);

  G = max(cl);
  mu = matrix(NA, G, d);

  sivy = iv[,,1]%*% y[1,];
  siv = iv[,,1];

 
  for (k in 1:G) {
    s = (cl==k);
    mu[k,] = compute.center(y[s,], iv[,,s]);
  }

  mu
}

## Assign points to the closed cluster
##
## @param y observations, m by d.
## @param mu, centers of the clusters
##
## output: cl, classification
mError.2 = function(y, iv, mu){
  dist = compute.distance.matrix(y, iv, mu);
  m = nrow(y);
  cl =  numeric(m);
  for(i in 1:m){
    cl[i] = order(dist[i,])[1];
  }
  cl
}

mError = function(y, iv, cl, niter=10) {
  for (i in 1:niter) {
    mu = mError.1(y, iv, cl);
    cl.old = cl;
    cl = mError.2(y, iv, mu);

    if (all(cl.old==cl)) {
      return(list(mu=mu, classification=cl, niter=i));
    }

  }

  return(list(mu=mu, classification=cl, niter=i));

}

examine.mError.stesp = function() {

  ##  Inspect steps in mError
  y = obs;
  iv = compute.iv(errary);
  cl = res0$classification;


  par(mfrow=c(2, 4));

  for (i in 1:8){ 
    plot(obs, col=cl);

    mu = mError.1(y, iv, cl);
    cl = mError.2(y, iv, mu);
  }

}
 
main = function(){
  y = obs;
  iv = compute.iv(errary);
  cl = res0$classification;

  res.mError = mError(y, iv, cl);
  cl1 = res.mError$classification;
  plot(obs, col=cl1);

  file.eps = sprintf("Graphs/real_data_seed%d_mError_results.eps", seed);
  postscript(file.eps, paper="special", width=6, height=6);

  cols = mclust.options("classPlotColors");
  pchs = mclust.options("classPlotSymbols");

  
  plot(obs,
       xlab="1h", ylab="3h", main="mError Clustering",
       col = cols[res.mError$classification],
       pch = pchs[res.mError$classification])

  dev.off()

 
}

