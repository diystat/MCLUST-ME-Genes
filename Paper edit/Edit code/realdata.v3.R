#######################################################################
###------------------------ RNA-Seq Example ------------------------###
#######################################################################

## 01/13/2020
##
## Robustness concerns: what if the error variance matrices are
## estimated wrong?  We will ranomly double or halve the estiamted
## variance-covariance matrices and see how the clustereing results
## change.
##
## 01/15/2020
##
## Or we can simualte another set of v.beta values by parametric
## boostrapping and fit our model with the modified error covariances.

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


###------------------- Importing Data ---------------------###

###################################################################
### Raw and processed data available in "Data" folder.          ###
###################################################################

# Import raw data
print(load("Data/rna_raw.RData"));
print(load("Results/1000.simulated.v.beta.RData"));

## Random sample 1000 genes, excluding rows with 0 relative frequnencies under any treatment
n.genes = nrow(nb.data$counts);

zero = apply(nb.data$rel.frequencies < 1/sum(nb.data$counts), 1, any);
m = 1000;

## seed=999 was the one used for the Genes paper.
seed = 999;
set.seed(seed);

ss = sample((1:n.genes)[!zero], m);

obs = full$beta[ss, 2:3] # data points
errary = v.beta[2:3,2:3,ss] # measurement error covariances
plot(obs);

## Run mclust
res.mclust = Mclust(obs, modelNames="VVV");


## Run MCME with original estimated errors and with two sets of modified errors.
## We will only run mcme for G=2.

## Generate initial group membership for mclust-me
generate.initial.group.membership = function(g){
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
  z.ini
}

z.ini = generate.initial.group.membership(2);

if (FALSE) {
  ## To confirm that the initial group memberships are deterministic, not random.
  z.ini.2 = generate.initial.group.membership(2);
  identical(z.ini, z.ini.2);
}

## 0. Run MCME with the originally estimated errors
res0 = mcmeVVV(obs, z.ini, errary);
res0$classification = 1 + (res0$z[,2]>0.5);

## The four cells should be 775, 10, 30, 185
table(res.mclust$classification, res0$classification);


## 1. Run MCME with modified errors 

## Now modify the error variance-covariance matrices
dim(errary);
err1 = errary;

## Generate a random vector of 0.5 or 2 values
seed1 = 999;
set.seed(seed1);
e1 = 2^((rbinom(m, 1, 0.5)*2)-1);

for (i in 1:m){
  err1[,,i] = errary[,,i] * e1[i];
}


res1 = mcmeVVV(obs, z.ini, err1);
res1$classification = 1 + (res1$z[,2]>0.5);

## 2. Run MCME with modified errors

## Now modify the error variance-covariance matrices
dim(errary);
err2 = errary;

## Generate a random vector of 0.5 or 2 values
seed2 = 999;
set.seed(seed2);

e2 = 2^(rbinom(m, 1, 0.5)-0.5);

for (i in 1:m){
  err2[,,i] = errary[,,i] * e2[i];
}
res2 = mcmeVVV(obs, z.ini, err2);
res2$classification = 1 + (res2$z[,2]>0.5);

## 3. Run MCME with modified errors (same as 2, but with different seed)

## Now modify the error variance-covariance matrices
dim(errary);
err3 = errary;

## Generate a random vector of 0.5 or 2 values
seed3 = 1001;
set.seed(seed3);

e3 = 2^(rbinom(m, 1, 0.5)-0.5);

for (i in 1:m){
  err3[,,i] = errary[,,i] * e3[i];
}
res3 = mcmeVVV(obs, z.ini, err3);
res3$classification = 1 + (res3$z[,2]>0.5);

## 01/15/2020, load previous results, then run additional sets
if (FALSE) {
  file.results = sprintf("Results/res_real2d_v3_m%d_seed%d.Rdata", m, seed);
  print(load(file.results));
}

## 4. Run MCME with simulated (bootstrap) erros in v.beta.sim
dim(v.beta.sim);
err4 = v.beta.sim[2:3,2:3,,1];
res4 = mcmeVVV(obs, z.ini, err4);
res4$classification = 1 + (res4$z[,2]>0.5);

## 5. Run MCME with simulated (bootstrap) erros in v.beta.sim
err5 = v.beta.sim[2:3,2:3,,2];
res5 = mcmeVVV(obs, z.ini, err5);
res5$classification = 1 + (res5$z[,2]>0.5);

## 6. Run MCME with simulated (bootstrap) erros in v.beta.sim
err6 = v.beta.sim[2:3,2:3,,3];
res6 = mcmeVVV(obs, z.ini, err6);
res6$classification = 1 + (res6$z[,2]>0.5);


## Save all results
file.results = sprintf("Results/res_real2d_v3_m%d_seed%d.%s.Rdata", m, seed, Sys.Date());
## save(m, seed, ss, obs, errary, seed1, err1, seed2, err2, seed3, err3, res.mclust, res0, res1, res2, res3, file=file.results);

save(m, seed, ss, obs, errary, seed1, err1, seed2, err2, seed3, err3, err4, err5, err6,
     res.mclust, res0, res1, res2, res3, res4, res5, res6, file=file.results);

## Quick comparison
table(res.mclust$classification, res0$classification);

table(res.mclust$classification, res1$classification);
table(res.mclust$classification, res2$classification);
table(res.mclust$classification, res3$classification);
table(res.mclust$classification, res4$classification);
table(res.mclust$classification, res5$classification);
table(res.mclust$classification, res6$classification);

table(res1$classification, res0$classification);
table(res2$classification, res0$classification);
table(res3$classification, res0$classification);
table(res4$classification, res0$classification);
table(res5$classification, res0$classification);
table(res6$classification, res0$classification);

table(res5$classification, res4$classification);
table(res6$classification, res4$classification);
table(res6$classification, res5$classification);

adjustedRandIndex(res.mclust$classification, res0$classification);
adjustedRandIndex(res.mclust$classification, res1$classification);
adjustedRandIndex(res.mclust$classification, res2$classification);
adjustedRandIndex(res.mclust$classification, res3$classification);

adjustedRandIndex(res1$classification, res0$classification);
adjustedRandIndex(res2$classification, res0$classification);
adjustedRandIndex(res3$classification, res0$classification);
adjustedRandIndex(res4$classification, res0$classification);
adjustedRandIndex(res5$classification, res0$classification);
adjustedRandIndex(res6$classification, res0$classification);
