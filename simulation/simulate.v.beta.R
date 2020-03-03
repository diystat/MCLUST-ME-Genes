rm(list());

#######################################################################
###------------------------ RNA-Seq Example ------------------------###
#######################################################################

setwd("/home/yanming/ongoing/Model-based-Clustering-Research");
## setwd("D:/ongoing/Model-based-Clustering-Research");

library(NBPSeq);
library(MASS);

###------------------- Importing Data ---------------------###

###################################################################
### Raw and processed data available in "Data" folder.          ###
###################################################################

# Import raw data
print(load("Data/rna_raw.RData"));

## Simulate 1000 rows of data according to the fitted model.
str(full);


## Random sample 1000 genes, excluding rows with 0 relative frequnencies under any treatment
n.genes = nrow(nb.data$counts);
zero = apply(nb.data$rel.frequencies < 1/sum(nb.data$counts), 1, any);

m = 1000;
seed = 999;
set.seed(seed);
ss = sample((1:n.genes)[!zero], m);

nb.data.0 = nb.data[ss,];
mu0 = full$mu[ss,];
phi0 = full$phi[ss];
v.beta.0 = v.beta[,,ss];

if (FALSE) {
  ## Sanity check
  beta = full$beta[ss,];
  phi = full$phi[ss];
  mu = t(exp(x %*% t(beta)) * nb.data.0$eff.lib.sizes);
  all.equal(mu, mu0);

  all.equal(nb.data.0$rel.frequencies, t(t(nb.data.0$counts)/nb.data.0$eff.lib.sizes));
}

## There is a function with the same name in NBPSeq! (That one needs to be removed or renamed.)

##' Simuldate a nb data set based on an existing nb data and estiamted mu and phi
##'
##' @title Simuldata a new nb data set
##' @param nb.data.0 a current nb.data set
##' @param mu0 a mxn matrix of mu values
##' @param ph0 a m vector of dispersion parameters
simulate.nb.data = function(nb.data.0, mu0, phi0) {
  n = length(mu0);
  y = rnbinom(n, size = 1/phi0, mu=mu0);
  nb.data.1 = nb.data.0;
  nb.data.1$counts[,] = y;
  nb.data.1$rel.frequencies = t(t(nb.data.1$counts)/nb.data.1$eff.lib.sizes);
  nb.data.1 
}

test.simulate.nb.data() = function() {
  ## To simulate a nb data set, we only need mu and kappa = 1/phi
  set.seed(1);
  nb.data.1 = simulate.nb.data(nb.data.0, mu0, phi0);

  cor(nb.data.1$counts, nb.data.0$counts);

  ## See whether the simualted data set look like the original.
  nb.data.0[1:5, 1:10];
  nb.data.1[1:5, 1:10];

}

compute.v.beta = function(nb.glm) {
  m = nrow(nb.glm$mu);
  v.beta = array(NA, c(5, 5, m));

  ## sd.beta = matrix(NA, m, 5);

  for (i in 1:m) {
    ## for (i in ss) {
    j = nb.glm$j[,,i];

    if (!any(is.na(j))) {
      v.beta[,,i] = ginv(j)[-1,-1][1:5, 1:5];
      ## sd.beta[i,] = sqrt(diag(v.beta[,,i]));
    }
  }

  v.beta
}

main() {
  
  ## Compare the nb.glm on the original data sets to the saved results: they should be the same
  nb.glm0 = fit.nb.glm(nb.data.0, x);
  all.equal(nb.glm0$mu, mu0);

  ## Compare the v.beta to the v.beta computed for the original data set
  v.beta.new = compute.v.beta(nb.glm0);
  all.equal(v.beta.new, v.beta.0);


  ## Simulate new nb data and compute the v.beta
  n.sim = 100;
  v.beta.sim = array(dim=c(5, 5, m, n.sim));
  
  
  for (i in 1:n.sim) {
    set.seed(i);
    nb.data.1 = simulate.nb.data(nb.data.0, mu0, phi0);
    nb.glm1 = fit.nb.glm(nb.data.1, x);
    v.beta.sim[,,,i] = compute.v.beta(nb.glm1);
  }

  dim(v.beta.sim)

  file.v.beta = "Results/1000.simulated.v.beta.Rdata";
  save(v.beta.0, v.beta.sim, file=file.v.beta);
  
}

analyze.v.beta = function() {
  file.v.beta = "Results/1000.simulated.v.beta.Rdata";
  print(load(file.v.beta));

  v.beta.0[1,1,1];
  v.beta.sim[1,1,1,];
  hist(v.beta.sim[1,1,1,]/ v.beta.0[1,1,1]);

  i = 1;
  i = 2;
  hist(v.beta.sim[2,2,i,]/ v.beta.0[2,2,i]);
  hist(v.beta.sim[3,3,i,]/ v.beta.0[3,3,i]);

  plot(v.beta.0[2,2,], v.beta.sim[2,2,,1], log="xy");

  lines(v.beta.0[2,2,], v.beta.0[2,2,]);

  i = 3
  par.old = par(mfrow=c(1,2));
  hist(v.beta.sim[2,2,,1]/v.beta.0[2,2,], nclass=30);
  hist(v.beta.sim[3,3,,1]/v.beta.0[3,3,], nclass=30);
  par(par.old):
  

  plot(v.beta.0[3,3,], v.beta.sim[3,3,,1]);
  lines(v.beta.0[3,3,], v.beta.0[3,3,]);

  lines(v.beta.0[2,2,], v.beta.0[2,2,]*0.5, log="xy");
  lines(v.beta.0[2,2,], v.beta.0[2,2,]*1.4, log="xy");

  plot(v.beta.0[2,2,], v.beta.sim[2,2,,2]/v.beta.0[2,2,], log="x");
  plot(v.beta.0[2,2,], v.beta.sim[2,2,,3]/v.beta.0[2,2,], log="x");

  abline(0,1);

}
