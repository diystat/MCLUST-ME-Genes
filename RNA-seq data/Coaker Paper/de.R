library(NBPSeq);
source('nb.glm.test.contrast.R');

file.counts = "../data/2014-01-30_counts.Rdata";
file.design = sprintf("../data/2014-01-30_design.Rdata", Sys.Date());
file.nb.data = "../data/2014-01-29_normalized.nb.data.Rdata";
file.nbq = "../output/2014-01-30_nbq.Rdata";
file.nbs = "../output/2014-01-30_nbs.Rdata";
file.group.model = "../output/2014-02-06_group.model.Rdata";
## file.test = "../output/2014-02-06.test.4.contrasts.res.Rdata";
file.test = "../output/2014-02-08.test.4.contrasts.res.Rdata";

## Results
## For Jeff
file.patterns = "../output/2014-02-08_patterns.top1000.G7.ps";
file.patterns.ps = "../output/2014-02-08_patterns.top1000.G7.ps";
file.mclust.l2fc = "../output/2014-02-08_mclust.l2fc.Rdata";
file.mclust.l2fc.ps = "../output/2014-02-08_mclust.l2fc.ps";

## For Eugene
file.l2fc = "../output/2014-02-08_matrix.txt";
file.cor = "../output/2014-02-08_correlation.txt";


##' Create file.counts, file.design
##'
##' @title  Process the raw data
##' @return  save counts to file.count, save treatment, hpi to file.design
read.data = function() {

  rm(list=ls());
  
  ## 2013-12-04
  ## 2014-01-30
  file.data = "../data/raw_gene_counts.txt";

  ## Read the counts
  raw  = read.table(file.data, header=TRUE, row.names=1);
  counts = as.matrix(raw);

  ## tags = rownames(counts);
  n = dim(counts)[2];

  ## treatment types
  header = colnames(counts);
  treatment = character(n);
  treatment[grep("Untreated", header)]="Untreated";
  treatment[grep("H20", header)]="H20";
  treatment[grep("flg22", header)]="flg22";

  treatment = factor(treatment, levels=c("Untreated", "H20", "flg22"));

  ## Hour after treatment (as factor?)
  hpi = rep(c(0, 1/6, 1, 3, 6, 12), c(3, 6, 6, 6, 6, 6));

  file.counts = sprintf("../data/%s_counts.Rdata", Sys.Date());
  save(counts, file=file.counts);

  file.design = sprintf("../data/%s_design.Rdata", Sys.Date());
  save(treatment, hpi, file=file.design);

  file.out
}

##' file.nb.data
##'
##' @title Normalize the data and run preapre.nb.data
##' @return save nb.data to file.nb.data
summary.statistics = function() {
  ## 2013-01-16
  ## 2013-01-29
  
  rm(list=ls());

  source('de.R');
  
  load(file.counts);

  ## Prepare NB data structure without normalization
  nb.data = prepare.nb.data(counts, tags=NULL);

  idp = rowSums(nb.data$counts)>0;
  eps = 1/sum(nb.data$counts);
  boxplot(log(nb.data$rel.freq[idp,]+eps));

  file.eps = sprintf("../output/%s_nb.data.unnormalized.eps",Sys.Date());
  postscript(file.eps, width=18, pointsize=1);
  plot.nb.data(nb.data);
  dev.off();

  ## Gene IDs
  genes = rownames(counts);

  ## A list of stably-expressed genes
  m = dim(nb.data$counts)[1];
  stable = as.matrix(read.table("../data/stably-expressed-genes.txt"));
  ## id = match(stable, genes);
  id = match(stable[6:24], genes);
  id = id[!is.na(id)]

  ## See their relative frequencies);
  matplot(t(nb.data$rel.freq[id,]), type="l");

  ## This data set requires normalization
  sorted = -apply(-nb.data$rel, 2, sort);
  colSums(sorted[1:10,]);

  colSums(sorted[11:99,]);
  colSums(sorted[100:999,]);

  colSums(sorted[11:20,]);
  colSums(sorted[21:30,]);
  colSums(sorted[31:40,]);
  colSums(sorted[41:50,]);
  colSums(sorted[100:200,]);
  colSums(sorted[1000:2000,]);
  colSums(sorted[2000:3000,]);

  sorted = -apply(-nb.data$rel, 2, sort);
  colSums(head(sorted, 10));

  head(sorted, 20);
  colSums(head(sorted, 20));

  ## Estimate norm factors based on all genes
  norm.factors.1 = estimate.norm.factors(counts, method="AH2010");

  ## Estimate norm factors based on the selected set of stably expressed genes
  norm.factors.2 = estimate.norm.factors(counts[id,], lib.sizes=colSums(counts), method="AH2010");

  ## The results are "comparable".
  plot(norm.factors.1, norm.factors.2);
  abline(0, 1);

  ## Prepare NB data structure with estimated normalization factors
  nb.data = prepare.nb.data(counts, norm.factors = norm.factors.1, tags=NULL);
  ## See the relative frequencies of the set of stably expressed genes.
  matplot(t(nb.data$rel.freq[id,]), type="l");
  boxplot(log(nb.data$rel.freq[idp,]+eps));
   
  nb.data = prepare.nb.data(counts, norm.factors = norm.factors.2, tags=NULL);
  matplot(t(nb.data$rel.freq[id,]), type="l");
  boxplot(log(nb.data$rel.freq[idp,]+eps));

  ## See the expression of a set of randomly selected 20 genes
  set.seed(999);
  id = sample(m, 24);
  matplot(t(nb.data$rel.freq[id,]),  type="l");

  ## Estimate normalizaiton factors using AH2010
  norm.factors.3 = estimate.norm.factors(counts[idp,]+1, method="AH2010");
  plot(norm.factors.1, norm.factors.3);

  ## If we only normalize samples 4:15, the resulting norm factors will be similar
  norm.factors.4 = estimate.norm.factors(counts[, -(1:3)], method="AH2010");
  r = norm.factors.1[-(1:3)]/norm.factors.4;
  r / mean(r) - 1;
  
  ##
  nb.data = prepare.nb.data(counts, norm.factors = norm.factors.1, tags=NULL);
  nb.data$eff.lib.sizes;

  ## plot.nb.data
  file.eps = sprintf("../output/%s_nb.data.normalized.eps",Sys.Date());
  postscript(file.eps, height=18, width=18, pointsize=1);
  plot.nb.data(nb.data);
  dev.off();

  file.nb.data = sprintf("../data/%s_normalized.nb.data.Rdata", Sys.Date());
  save(nb.data, file=file.nb.data);

  file.out
}

##' @title Fit NBQ, NBS dispersion models
##' @return save nbq, nbs to file.nbq, file.nbs
dipsersion.models = function(){
  ## 2014-01-29
  rm(list=ls());
  source("de.R");

  ## Load normalized NB data
  vars = load(file.nb.data);

  ## Take a subset
  ## nb.data = prepare.nb.data(nb.data$counts[1:1000,], nb.data$eff.lib.sizes); 

  ## Estimate dispersion
  grp.ids = factor(rep(0:10, each=3));
  x = model.matrix(~0+grp.ids);

  ## NBQ
  ## debug(estimate.dispersion);
  ## debug(disp.nbq);
  ## debug(irls.nb);
  nbq = estimate.dispersion(nb.data, x, model="NBQ", control=list(trace=3));
  file.nbq = sprintf("../output/%s_nbq.Rdata", Sys.Date());
  save(nbq, file=file.nbq);

  nbs = estimate.dispersion(nb.data, x, model="NBS", control=list(trace=3));
  ## system.time({nbq = estimate.dispersion(nb.data, x, model="NBQ", control=list(trace=6));});
  file.nbs = sprintf("../output/%s_nbs.Rdata", Sys.Date());
  save(nbs, file=file.nbs);

  nbq;
  plot(nbq);
  plot(nbs);

  quantile(nbq$estimates, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99));

  z = nbq$model$pi.pre * 1e6;
  smart.plot(z, nbq$est, log="xy", clip=64, ylim=c(1e-2, 1));

  hist(log10(nbq$model$pi.pre));

  invisible();
}

##' @title Fit group model
##' @return save res0 to file.group.model
fit.group.model =  function(){
  ## 2014-02-06
  rm(list=ls());
  source("de.R");

  ## Load NB data
  var.name = load(file.nb.data);

  m = dim(nb.data$counts)[1];
  n = dim(nb.data$counts)[2];

  ## Load design
  var.name = load(file.design);

  ## Load dispersion estimates
  var.name = load(file.nbq);

  ## Fit NB regressions to all genes
  ## res = fit.nb.regression(nb.data, disp, x);

  ## Specify the model matrix with group structure
  grp.ids = factor(rep(0:10, each=3));
  x = model.matrix(~0+grp.ids);
  p = 11;

  res0 = irls.nb(nb.data$counts, nb.data$eff.lib.sizes, x, nbq$estimates, beta0 = rep(NA, dim(x)[2]));

  file.group.model = sprintf("../output/%s_group.model.Rdata", Sys.Date());

  save(res0, file = file.group.model);

}

##' Run four separate tests to compare the log fold changes between
##' 10min and 1 hour, 1 hour and 3 hours, 
##' @title comparing log fold changes between consecutive time points 
##' @param pattern a string vector giving the alternative hypotheses,
##' each component can be "greater", "less" or "two.sided". 
##' @return a matrix of p-values
test.pattern = function(pattern) {
  ## 2014-01-30
  
  ## rm(list=ls());
  source("de.R");
  pattern = c("greater", "less", "less", "less"); 

  ## Load NB data
  file.nb.data = "../data/2014-01-29_normalized.nb.data.Rdata";
  var.name = load(file.nb.data);

  m = dim(nb.data$counts)[1];
  n = dim(nb.data$counts)[2];

  ## Load design
  var.name = load(file.design);

  ## Load dispersion estimates
  var.name = load(file.nbq);

  ## Specify the model matrix
  grp.ids = factor(rep(0:10, each=3));
  x = model.matrix(~0+grp.ids);
  p = 11;

  p.values = matrix(0, m, 4);

  for (i in 1:4) {
    ## Test (5-4) greater than (3-2)
    ## Test (7-6) less than (5-4)
    ## Test (9-8) less than (7-6)
    ## Test (11-10) less than (9-8)

    A = numeric(p);
    cols = (0:3) + 2*i
    A[cols] = c(1, -1, -1, 1);
    debug(test.contrasts);
    ## debug(complete.matrix);
    res1 = test.contrasts(nb.data, nbq, x, A, 0, tests = c("HOA", "LR", "Wald"), alternative=pattern[i]);
    p.values[,i] = res1$HOA$p.values;
  }

  p.values
}

##' Run one single composite test to compare the log fold changes between
##' 10min, 1 hour, 3 hours, 6 hours and 12 hours 
##'
##' The test is "two.sided".
##' 
##' @title comparing log fold changes between different time points 
##'
##' @return a list, output from test.contrasts (which calls
##' test.coefficient after applying a linear transformation on the
##' design matrix with group strucutre)
test.multiple.contrasts = function() {
  ## 2014-02-06
  
  ## rm(list=ls());
  ## source("de.R");

  ## Load NB data
  var.name = load(file.nb.data);

  m = dim(nb.data$counts)[1];
  n = dim(nb.data$counts)[2];

  ## Load design
  var.name = load(file.design);

  ## Load dispersion estimates
  var.name = load(file.nbq);

  ## Specify the model matrix
  grp.ids = factor(rep(0:10, each=3));
  x = model.matrix(~0+grp.ids);
  p = 11;

  p.values = numeric(m);

  ## The contrast matrix
  q = 4;
  A = matrix(0, q, p);

  for (i in 1:4) {
    ## Test (5-4) greater than (3-2)
    ## Test (7-6) less than (5-4)
    ## Test (9-8) less than (7-6)
    ## Test (11-10) less than (9-8)

    cols = (0:3) + 2*i
    A[i, cols] = c(1, -1, -1, 1);
  }

  ## debug(complete.matrix);
  res = test.contrasts(nb.data, nbq, x, A, 0, tests = c("HOA", "LR", "Wald"));

  if (FALSE) {
    p.hoa = res$HOA$p.values;
    p.lr = res$LR$p.values;
    plot(p.hoa, p.lr, log="xy");
    plot(p.hoa, p.lr);
    abline(0,1);
    hist(log(p.hoa/p.lr));
  }
       
  res
}

main = function() {

  rm(list=ls());
  source("de.R");

  ## 2014-02-06
  ## 2014-02-08
  res = test.multiple.contrasts();

  file.res = sprintf("../output/%s.test.4.contrasts.res.Rdata", Sys.Date());
  save(res, file=file.res);

  if (FALSE) {
    ## compare to earlier results
    res.new = res;
    file.test = "../output/2014-02-06.test.4.contrasts.res.Rdata";
    load(file.test);
    all.equal(res.new$LR$p.values, res.new$LR$p.values);
    all.equal(res.new$HOA$p.values, res.new$HOA$p.values);
    all.equal(res.new$Wald$p.values, res.new$Wald$p.values);
  }

  ## 2014-02-01
  p.values = test.pattern(c("greater", "less", "less", "less"));
  file.p = sprintf("../output/%s_uddd.p.values.Rdata", Sys.Date());
  save(p.values, file=file.p);

  p.values = test.pattern(c("less", "less", "less", "less"));
  file.p = sprintf("../output/%s_dddd.p.values.Rdata", Sys.Date());
  save(p.values, file=file.p);

  p.values = test.pattern(c("greater", "greater", "greater", "greater"));
  file.p = sprintf("../output/%s_uuuu.p.values.Rdata", Sys.Date());
  save(p.values, file=file.p);

}

##' @title summarize results from the composite test for equal log fold change over time
##' @return 
summarize.results = function()  {
  ## 2014-02-06
  ## 2014-02-08
  rm(list=ls());
  source("de.R");

  ## Load NB data
  var.name = load(file.nb.data);

  ## Load group model
  var.name = load(file.group.model);

  ## Load test results
  var.name = load(file.test);

  ## Read LR test p-values
  p.lr = res$LR$p.values;

  ## Contrasts
  contrasts = res$beta.hat[,1:4];

  ## Log fold changes at 1/6, 1, 3, 6, 12 hpi
  e = res0$beta;
  log.fc = e[,c(3, 5, 7, 9, 11)] - e[, c(2, 4, 6, 8, 10)];
  rownames(log.fc) = rownames(nb.data$counts);
  colnames(log.fc) = c("10min", "1hr", "3hr", "6hr", "12hr");
  
  or = order(p.lr);
  top = or[1:1000];
  l2fc = log.fc[top,]/log(2);

  ## clustering analysis on log fc
  library(mclust)
  args(Mclust);

  ## debug(mstepVVV);
  ## debug(estepVVV);
  ## undebug(mstepVVV);

  ## undebug(meVVV);
  ## debug(mstepVVV);
  cl = Mclust(l2fc);

  file.mclust = sprintf("../output/%s_mclust.l2fc.Rdata", Sys.Date());
  save(cl, file=file.mclust);
  file.mclust.l2fc = sprintf("../output/%s_mclust.l2fc.ps", Sys.Date());
  postscript(file.mclust.l2fc, paper="letter");
  plot(cl);
  dev.off();

  ## cl2 = Mclust(l2fc);
  ## identical(cl, cl2);

  ## cl = Mclust(l2fc, G=8);

  file.patterns.ps = sprintf("../output/%s_patterns.top1000.G%d.ps", Sys.Date(), cl$G);
  postscript(file.patterns.ps);
  ## hpi = c(1/6, 1, 3, 6, 12);
  ## matplot(hpi, t(l2fc), type="l", col=cl$classification);
  matplot(t(l2fc), type="l", col=cl$classification,
          ylab="log2 fold change");

  par(mfrow=c(2, 4));
  for (i in 1:cl$G) {
    id = cl$classification==i;
    matplot(t(l2fc[id,]), type="l", col=i, ylim=c(-7,7), main=i,
            ylab="log2 fold change");
  }
  dev.off();

  ## Save log fc and the correlation matrix (for Eugene)
  file.l2fc = sprintf("../output/%s_matrix.txt", Sys.Date());
  file.cor = sprintf("../output/%s_correlation.txt", Sys.Date());
  write.table(l2fc, file=file.l2fc, quote=FALSE);
  s = cor(t(l2fc));
  write.table(s, file=file.cor, quote=FALSE);

  ## Save results for Jeff
  obj =  cbind(l2fc, p.values=p.lr[top], classification = cl$classification);
  table(obj[,7]);
  file.patterns = sprintf("../output/%s_patterns.top1000.G%d.csv", Sys.Date(), cl$G);
  write.csv(obj, file=file.patterns, quote=FALSE);

  ## Older analyses

  top = or[1:100];

  ## Look at the patterns
  file.patterns = sprintf("../output/%s_patterns.top100.eps", Sys.Date());
  postscript(file.patterns, paper="special", width=8, height=8);
  matplot(t(log.fc[top,]), type="l", xlab="hpi", ylab="log fold change");
  dev.off();

  c1 = contrasts[top,]
  id = c1[,1]>0 & c1[,2]<0;
  mean(id);
  c1[!id,];

  file.patterns = sprintf("../output/%s_patterns.top100.others.eps", Sys.Date());
  ## postscript(file.patterns, paper="special", width=8, height=8, horizontal=FALSE);
  postscript(file.patterns, paper="special", width=8, height=8);
  matplot(t(log.fc[top,][!id,]), type="l", xlab="hpi", ylab="log fold change");
  dev.off();

  ## Gene IDs
  genes = rownames(nb.data$counts);

  ## A list of stably-expressed genes
  m = dim(nb.data$counts)[1];
  stable = as.matrix(read.table("../data/stably-expressed-genes.txt"));
  ## id = match(stable, genes);
  ids = match(stable[6:24], genes);
  ids = ids[!is.na(ids)]
  file.patterns = sprintf("../output/%s_patterns.stable.eps", Sys.Date());
  postscript(file.patterns, paper="letter", horizontal=FALSE);
  matplot(t(log.fc[ids,]), type="l", xlab="hpi", ylab="log fold change",
          ylim=c(-1,1));
  dev.off();

  library(mclust)

  ## clustering analysis of log fold changes
  cl = Mclust(lfc);
  file.mclust.lfc = sprintf("../output/%s_mclust.lfc.ps", Sys.Date());
  postscript(file.mclust.lfc);
  plot(cl);
  dev.off();

  ## clustering analysis of contrasts
  ctr = contrasts[top,];

  args(Mclust);
  cl.ctr = Mclust(ctr);
  cl.ctr.4 = Mclust(ctr, 4);

  file.mclust.ctr = sprintf("../output/%s_mclust.ctr.4.ps", Sys.Date());
  postscript(file.mclust.ctr);
  plot(cl.ctr.4);
  dev.off();

  ## PCA on log fold change
  ## boxplot(log.fc);
  pc = prcomp(log.fc[top,]);
  plot(pc);
  biplot(pc);

  ## PCA on contrasts
  pc = prcomp(ctr);
  biplot(pc);

  cl.pc = Mclust(ctr %*% pc$rotation);

  ## matplot(t(contrasts[top,]), type="l");
  plot(ctr[,1:2]);
  plot(contrasts[top, 1:2]);
  

}

##' @title compare the composite tests with the four separate tests
##' @return 
summarize.results.2 = function()  {
  ## 2014-02-06
  rm(list=ls());
  source("de.R");

  ## Load NB data
  var.name = load(file.nb.data);

  m = dim(nb.data$counts)[1];
  n = dim(nb.data$counts)[2];

  ## Load design
  var.name = load(file.design);

  ## Load dispersion estimates
  var.name = load(file.nbq);

  ## Specify the model matrix with group structure
  grp.ids = factor(rep(0:10, each=3));
  x = model.matrix(~0+grp.ids);
  p = 11;

  ## Load group model
  var.name = load(file.group.model);

  ## Load test results
  var.name = load(file.test);

  p.lr = res$LR$p.values;

  hist(p.lr);
  plot(p.lr, p.values[,2], log="xy");
  abline(0,1);

  or = order(p.lr);
  top = or[1:100];

  contrasts = res$beta.hat[,1:4];
  contrasts[top,];

  matplot(t(contrasts[top,]), type="l");

  plot(contrasts[top, 1:2]);
  
  ## Log fold changes at 1/6, 1, 3, 6, 12 hpi
  e = res0$beta;
  log.fc = e[,c(3, 5, 7, 9, 11)] - e[, c(2, 4, 6, 8, 10)];
  log.fc[top,];
  matplot(t(log.fc[top,]), type="l");


}
summarize.results.3 = function()  {
  ## 2014-02-01
  rm(list=ls());
  source("de.R");

  ## Load NB data
  file.nb.data = "../data/2014-01-29_normalized.nb.data.Rdata";
  var.name = load(file.nb.data);

  m = dim(nb.data$counts)[1];
  n = dim(nb.data$counts)[2];

  ## Load design
  file.design = sprintf("../data/2014-01-30_design.Rdata", Sys.Date());
  var.name = load(file.design);

  ## Load dispersion estimates
  file.disperison = "../output/2014-01-30_nbq.Rdata";
  var.name = load(file.disperison);

  ## Fit NB regressions to all genes
  ## res = fit.nb.regression(nb.data, disp, x);

  ## Specify the model matrix
  grp.ids = factor(rep(0:10, each=3));
  x = model.matrix(~0+grp.ids);
  p = 11;

  res = irls.nb(nb.data$counts, nb.data$eff.lib.sizes, x, nbq$estimates, beta0 = rep(NA, dim(x)[2]));

  ## Load p values
  ## file.p = "../output/2014-01-31_uddd.p.values.Rdata";
  ## var.name=load(file.p);
  ## p.values.old = p.values;

  file.p = "../output/2014-02-02_uddd.p.values.Rdata";
  var.name=load(file.p);
  p.values.0 = p.values;

  if (FALSE) {
    identical(p.values.old, p.values);
    all.equal(p.values.old, p.values);

    range(p.values.old/p.values);
    
    or = order(-apply(p.values.old-p.values, 1, max));
    cbind(p.values.old, p.values)[or[1:10],];
    (p.values.old - p.values)[or[1:10],];
  }

  file.p = "../output/2014-02-02_dddd.p.values.Rdata";
  var.name=load(file.p);
  p.values.1 = p.values;

  file.p = "../output/2014-02-02_uuuu.p.values.Rdata";
  var.name=load(file.p);
  range(p.values + p.values.1);

  hist(log(p.values[,1]));
  hist(log(p.values[,2]));
  hist(log(p.values[,3]));
  hist(log(p.values[,4]));

  plot(p.values[,1:2], log="xy");
  abline(0,1);

  pairs(log(p.values));

  ## Fisher's method for combining p-values
  fisher = pchisq(-rowSums(log(p.values)), df=2*4, lower.tail=FALSE );
  hist(fisher);
  or = order(fisher);

  n.top = 50;
  top = or[1:n.top];
  nb.data$counts[top,];

  ## Log fold changes at 1/6, 1, 3, 6, 12 hpi
  e = res$beta;
  log.fc = e[,c(3, 5, 7, 9, 11)] - e[, c(2, 4, 6, 8, 10)];
  log.fc[top,];
  matplot(t(log.fc[top,]), type="l");


}


example.0 = function() {
  ##
  data(arab);

  nb.data = prepare.nb.data(arab);

  cor(log(nb.data$rel.freq + 1e-8));

  plot.nb.data(nb.data, clip = 16, res = 50);
  
}
