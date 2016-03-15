## Thu Nov  5 11:06:35 PST 2015

## Provide
## 
##   1. Normalized nb.data
##   2. treatment, hpi, and the model matrix x
##   3. Fitted full and reduced models
##   4. Estimated variance/covaraince matrices.
##   5. LR test results.

library(NBPSeq);
library(MASS);

## Input file names
file.nb.data = "../data/2014-01-29_normalized.nb.data.Rdata";
file.design = "../data/2014-01-30_design.Rdata";
file.nb.glm = "../output/2014-05-06.nb.glm.Rdata";

main = function() {
  ## 2014-05-07
  rm(list=ls());
  source("top.genes.with.CI.R");

  ## Load NB data
  print(load(file.nb.data));
  class(nb.data) = "nb.data";
  print(nb.data);

  m = nrow(nb.data$counts);
  n = ncol(nb.data$counts);
 
  ## Load design
  print(load(file.design));

  ## In this model, the first 5 columns of the regression coefficients
  ## correspond to fold changes. The dispersion was unknown.
  print(load(file.nb.glm));

  ## See the structure of fitted full model
  str(full);
  class(full)="nb.glm";

  ## Specify the model matrix
  
  ## The model matrix used to fit NB regression model: the first 5
  ## columns correspond the log fold changes at 5 time points and last
  ## 6 columns correspond to the mean expression levels at each time
  ## point
  p = 11;
  x = matrix(0, n, p);
  x[,6:11] = model.matrix(~0+ factor(hpi));
  v = rep(c(-0.5, 0.5), each=3); 
  for(i in 1:5) {
    x[(4:9) + (i-1)*6, i] = v;
  }

  ## In this model, the first 5 columns of the regression coefficients
  ## correspond to fold changes. The dispersion was unknown.
  print(load(file.nb.glm));

  ## See the structure of fitted full model
  str(full);
  full.old = full;
  reduced.old = reduced;

  full = fit.nb.glm(nb.data, x);
  str(full);

  beta0 = rep(NA, p);
  beta0[1:5] = 0;
  reduced  = fit.nb.glm(nb.data, x, beta0);
  reduced$beta0 = beta0;
  str(reduced);
  
  compare.mu = function(new, old) {
    re.mu = new$mu/old$mu - 1;
    plot(old$mu, re.mu, xlab="Fitted mu from old model", ylab="Relative error of mu between new and old models",
         log="x");
    print("Range of relative error when fitted mu (from old model) is greater than 0.01:");
    print(range(re.mu[old$mu > 0.01]));
    invisible(re.mu);
  }

  ## Verify that the new and old models give consistent results
  par(mfrow=c(1,2));
  compare.mu(full, full.old);
  compare.mu(reduced, reduced.old);

  ## Compute CI

  ## Compute the variance-covariance matrix and sd of the first five
  ## regression coefficients (corresponding to log fold changes at five
  ## time points). The first method is correct.
  v.beta = array(NA, c(5, 5, m));
  sd.beta = matrix(NA, m, 5);

  ## DEBUG: when i = 12, some beta.hat tend to -Inf.
  ## for (i in ss[1:100]) {
  for (i in 1:m) {
    ## for (i in ss) {
    j = full$j[,,i];

    if (!any(is.na(j))) {
      v.beta[,,i] = ginv(j)[-1,-1][1:5, 1:5];
      sd.beta[i,] = sqrt(diag(v.beta[,,i]));
    }
  }

  ## Perform LRT test
  t = 2 * (full$l - reduced$l);

  ## The type-I error is inflated
  p = pchisq(t, 5, lower.tail=FALSE);
  hist(p);

  ## Identify top 1000 genes
  filter = apply(full$mu>1, 1, all);
  table(filter);
  t[!filter] = NA;

  or = order(t, decreasing=TRUE);
  top = or[1:1000];

  t[top];
 
  file.out = sprintf("../output/%s.Coaker.Rdata", Sys.Date());
  save(nb.data, treatment, hpi, x, full, reduced, v.beta, sd.beta, t, or, top, file=file.out);

}


plot.top.genes = function() {
  ## Fri Nov  6 09:48:24 PST 2015
  rm(list=ls());
  dd = "2015-11-06";
  file.out = sprintf("../output/%s.Coaker.Rdata", dd);
  print(load(file.out));

  full$beta[top, 1:5];
  v.beta[,,top];
  sd.beta[top,];

  ## MA plots at five time points and highlight top 100 DE genes
  id = filter = apply(full$mu>1, 1, all);
  top100 = or[1:100];
  
  par(mfrow=c(2, 3));
  xlabs = c("10 min", "1 hr", "3 hr", "6 hr", "12 hr");
  for (i in 1:5) {
    smart.plot(full$beta[id,i+6]/log(10), full$beta[id,i]/log(2), pch=19, clip=32,
               xlab = xlabs[i],
               ylab = "log (base 2) fold change" 
               );
    smart.points(full$beta[top100,i+6]/log(10), full$beta[top100,i]/log(2), col="cyan", pch=19);
  }

  ## matplot
  x = c(1/6, 1, 3, 6, 12);
  top = or[1:10];
  matplot(x, t(full$beta[top, 1:5]), type="l", col=1);

  ##
  require(ggplot2)
  par(mfrow=c(2,5));
  top = or[2:11];

  for (i in top) {

    ## i = top[2];
    v = solve(full$j[, , i])[-1,-1][1:5,1:5];
    sd = sqrt(diag(v));

    y = full$beta[i,1:5];
    L = y - sd * 2;
    U = y + sd * 2;
    matplot(x , cbind(L, y, U), type=c("l", "b", "l"), pch=1, lty=c(2, 1, 2), col=1);

    if (FALSE) {

      matpoints(x , cbind(L, U));

      graphics.off();

      df = data.frame(x=y, y=y, L=L, U=U);
      ggplot(df, aes(x = x, y = y)) +
        geom_point(size = 4) +
          geom_errorbar(aes(ymax = U, ymin = L))

    }
  }

}


  



compute.variance = function() {
  ## 2014-05-07
  rm(list=ls());
  source("summary.R");

  ## Load NB data
  print(load(file.nb.data));
  class(nb.data) = "nb.data";
  print(nb.data);

  m = nrow(nb.data$counts);
  n = ncol(nb.data$counts);
 
  ## Load design
  print(load(file.design));

  ## In this model, the first 5 columns of the regression coefficients
  ## correspond to fold changes. The dispersion was unknown.
  print(load(file.nb.glm));

  ## See the structure of fitted full model
  str(full);
  class(full)="nb.glm";

  ## The model matrix used to fit NB regression model: the first 5
  ## columns correspond the log fold changes at 5 time points and last
  ## 6 columns correspond to the mean expression levels at each time
  ## point
  p = 11;
  x = matrix(0, n, p);
  x[,6:11] = model.matrix(~0+ factor(hpi));
  v = rep(c(-0.5, 0.5), each=3); 
  for(i in 1:5) {
    x[(4:9) + (i-1)*6, i] = v;
  }
  
  ## Compute variance from the observed information matrix
    
  ## Compare there ways to compute the CI:  inverse, generalized
  ## inverse, and inverting the diaganoal elements

  ## Filter rows with all zeros
  ## ss = (1:m)[rowSums(nb.data$counts)>0];
  ss = (1:m)[apply(full$mu>1, 1, all)];
  length(ss);

  ## Compute the variance-covariance matrix and sd of the first five
  ## regression coefficients (corresponding to log fold changes at five
  ## time points). The first method is correct.
  v1 = array(NA, c(5, 5, m));
  v2 = array(NA, c(5, 5, m));

  sd1 = matrix(NA, m, 5);
  sd2 = matrix(NA, m, 5);

  ## DEBUG: when i = 12, some beta.hat tend to -Inf.
  for (i in ss[1:100]) {
  ## for (i in ss) {
    j = full$j[,,i];
    v1[,,i] = ginv(j)[-1,-1][1:5, 1:5];
    sd1[i,] = sqrt(diag(v1[,,i]));
  }

  ## If we treat the estimated dispersion as known, the observed info is always diagonal?
  ## Yes for the log fold changes.
  for (i in ss[1:100]) {
  ## for (i in ss[1:100]) {
    j = full$j[,,i];
    v2[,,i] = solve(j[-1,-1])[1:5,1:5];
    sd2[i,] = sqrt(diag(v2[,,i]));
  }

  range(sd2-sd1, na.rm=TRUE);
  range(sd2/sd1, na.rm=TRUE);

  i = 7;
  obj = full[7,];

  y = c(1, 1, 1, 0, 0, 0);

  eps = 1e-8;

  x2 = matrix(0, 6, 2);
  x2[,1] = 1;
  x2[,2] = c(0.5, 0.5, 0.5, -0.5, -0.5, -0.5);

  ## debug(fit.nb.glm.1u);
  ## undebug(fit.nb.glm.1u);
  y = c(1, 2, 3, eps, eps, eps);

  y = c(1, 2, 3, 0, 0, 0);
  o = fit.nb.glm.1u(y, s=1e6, x2);
  o$mu;
  solve(o$j[-1,-1]);
  ginv(o$j)[-1,-1];

  ##

  "Done."
}

