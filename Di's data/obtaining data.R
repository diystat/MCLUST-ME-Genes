


library(MASS);

## File names
file.nb.data = "data/2014-01-29_normalized.nb.data.Rdata";
file.design = "data/2014-01-30_design.Rdata";
file.nb.glm = "data/2014-05-06.nb.glm.Rdata";

compute.variance = function() {
  ## 2014-05-07
  rm(list=ls());
  source("summary.R");

  ## Load NB data
  print(load(file.nb.data));
  class(nb.data) = "nb.data";
  ## print(nb.data);
  m = nrow(nb.data$counts);
  n = ncol(nb.data$counts);
  m;n
 
  ## Load design
  print(load(file.design));

  ## In this model, the first 5 columns of the regression coefficients
  ## correspond to fold changes. The dispersion was unknown.
  print(load(file.nb.glm));

  ## See the structure of fitted full model
  str(full);

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
    
  ## Compare three ways to compute the CI:  inverse, generalized
  ## inverse, and inverting the diaganoal elements

  ## Filter rows with any estimated mean < 1
  ss = (1:m)[apply(full$mu>1, 1, all)];
  length(ss);

  ## Compute the variance-covariance matrix and sd of the first five
  ## regression coefficients (corresponding to log fold changes at five
  ## time points). The first method is correct.
  v1 = array(NA, c(5, 5, m));
  sd1 = matrix(NA, m, 5);

  ## DEBUG: when i = 12, some beta.hat tend to -Inf.
  for (i in ss[1:100]) {
  ## for (i in ss) {
    j = full$j[,,i];
    v1[,,i] = ginv(j)[-1,-1][1:5, 1:5]; # get rid of dispersion parameter(1st row/col)
    sd1[i,] = sqrt(diag(v1[,,i]));
  }

}





## Find top genes using LR test
find.top.genes = function() {
  ## 2014-05-08
  rm(list=ls());
  source("summary.R");

  ## Load NB data
  print(load(file.nb.data));
  class(nb.data) = "nb.data";
  ## print(nb.data);

  ## Load design
  print(load(file.design));

  ## In this model, the first 5 columns of the regression coefficients
  ## correspond to fold changes. The dispersion was unknown.
  print(load(file.nb.glm));

  ## In this model, each column of the regression coefficient
  ## correspond to one group. The dipsersion was estimated from a NBQ
  ## model.
  ## 
  ## print(load(file.group.model));

  ## LRT statistic:
  t = 2 * (full$l - reduced$l);

  ## Assign NAs to test stat whose mean level <=1
  filter = apply(full$mu>1, 1, all);
  t[!filter] = NA;

  ## The type-I error is inflated
  p = pchisq(t, 5, lower.tail=FALSE); # obtain p-value
  hist(p);

  ## Obtain top 500 most DE genes
  or = order(-t);
  top = or[1:500];
  
 

}


  ## obtain observations:
  coef = full$beta[top,(1:5)]
  
  ## obtain estimation error matrices:
  err = v1[,,top]

  save(coef, err, file = "rnaseq_data.RData")


  print(load("rnaseq_data.RData"))


















