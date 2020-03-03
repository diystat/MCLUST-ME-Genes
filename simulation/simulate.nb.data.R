## DEBUG=FALSE;

##' Simulate a single data set according to the specified NB regression model
##'
##' @title Simulate a single data set according to the specified NB regression model
##' @param kappa shape (size) parameter, 1/kappa is the dispersion parameter
##' @param s a n-vector, library sizes;
##' @param x a nxp design matrix
##' @param beta regression coefficients
##' @return a n-vecotr of NB counts
##' @author Yanming Di
simulate.one = function(kappa, s, x, beta) {
  mu = s * exp(x %*% beta);
  n = length(mu);
  rnbinom(n, size = kappa, m=mu);
}

##' Simulate NB data sets accroding to specified NB model
##'
##' @title Simulate NB data sets accroding to specified NB model
##' @param m number of data sets to simulate
##' @param model a list specifying the NB regression model (see the
##' parameter list of \code{\link{simulate.one}} for more details).
##' @param seed random number seed
##' @return a data matrix with m rows--each corresponding to one data set
##' @author Yanming Di
simulate.data = function(m, model, seed) {
  kappa =  model$kappa;
  s = model$s;
  x = model$x;
  beta = model$beta;

  mu = s * exp(x %*% beta);
  n = length(mu);

  set.seed(seed);

  yy = matrix(NA, m, n);
  for (i in 1:m) {
    yy[i,] = rnbinom(n, size = kappa, m=mu);
  }

  yy
}

specify.two.group.model = function(n1, n2, pi1, beta2, phi,
  s=rep(1e6, n)) { 

  n = n1+n2;

  ## Corresponding parameters in a NB regression
  kappa = 1/phi
  x = matrix(0, n, 2);
  x[,1]=1;
  x[(n1+1):n,2]=1;

  beta = numeric(2);

  ## True model
  beta[1] = log(pi1);
  beta[2] = beta2;

  ## Specify the hypothesis
  beta0 = c(NA, 0);

  if (beta[2]==0) {
    key = sprintf("n%d.%d.mu%d.phi%.2f", n1, n2, pi1*s[1], phi);
  } else {
    key = sprintf("n%d.%d.mu%d.phi%.2f.fc%.1f", n1, n2, pi1*s[1], phi, exp(beta[2]));
  }

  list(kappa=kappa, beta=beta, s=s, x=x, beta0=beta0, key=key);
}

##' Simple NB regression model with one covariate (null)
##'
##' @title 
##' @return 
##' @author Yanming Di
regression.1d.null = function(xvals, mu1, phi, s = rep(1e6, n)) {

  kappa = 1/phi

  ## design matrix
  x = model.matrix(~xvals);

  n = dim(x)[1];
  p = dim(x)[2];
  q = 1;

  ## Library sizes
  ## s =  (1:n) * 1e6;

  ## True parameter values
  beta = numeric(p);
  beta[1] = log(mu1/s[1]);
  beta[2] = 0;

  mu = s * exp(x %*% beta);
  dim(mu)=NULL;

  ## Null hypothesis
  beta0 = c(rep(NA, p-q), rep(0, q));

  ## Specify the hypothesis
  ## beta0 = c(NA, 0);

  key = sprintf("regression.mu%d.phi%.2f", mu1, phi);

  list(kappa=kappa, beta=beta, s=s, x=x, beta0=beta0, mu=mu, key=key);
}

## sim.model = regression.1d.null(c(1,2,4,8, 16, 32), phi=0.1);

##' Two treatments with two levels each. No interaction effect.
##' @title 
##' @return 
##' @author Yanming Di
##' @param n.reps 
##' @param mu0 baseline mean
##' @param fc1 effect of factor 1 (fold change)
##' @param fc2 effect of factor 2 (fold change)
##' @param phi 
two.by.two.null = function(n.reps, mu0, phi, s=rep(1e6,n), fc1=1/2, fc2=2) {
  kappa = 1/phi;
  
  p1 = 2;
  p2 = 2;
  n.reps = 3;

  ## balanced 2-way design
  dd = data.frame(a = gl(p1, p2, p1*p2*n.reps), b = gl(p2, 1, p1*p2*n.reps)) 

  ## design matrix
  x = model.matrix(~a + b + a*b, dd);

  n = dim(x)[1];
  p = dim(x)[2];
  q = (p1-1) * (p2-1);

  ## Library sizes

  ## True parameter values
  beta = numeric(p);
  beta[1] = log(mu0/s[1]);
  beta[2] = log(fc1);
  beta[3] = log(fc2);
  ## beta[2:p1] = seq(-1, 0.5, length=p1-1);
  ## beta[(p1+1):(p1+p2-1)] = seq(0.1, 1, length=p2-1);

  ## Null hypothesis
  beta0 = c(rep(NA, p-q), rep(0, q));

  ## Means
  mu = s * exp(x %*% beta);

  key = sprintf("2x2.reps%d.mu%.0f.phi%3.1f", n.reps, mu0, phi);

  list(kappa=kappa, s=s, x=x, beta=beta, beta0=beta0, mu=mu, key=key);
}

test.simulate.data = function() {
  m1 = specify.two.group.model(2, 4, 100/1e6, 0, 0.1);
  ## debug(simulate.data);
  yy = simulate.data(10000, m1, 999);
  head(yy);
  ## plot(rowMeans(yy[,1:2]), rowMeans(yy[,3:6]));
}
  
##' Simulate NB data and perform HOA, LR, and Wald tests on simulated
##' data. The dispersion is assumed unknown when performing the tests.
##'
##' @title Simulate NB data and perform HOA, LR, and Wald tests on simulated data
##' @param sim.model 
##' @param alternative 
##' @param n.sim 
##' @param sim.seed 
##' @param hoa.seed 
##' @param R 
##' @return a list of simulation parameters and p-values of the tests
##' @author Yanming Di
run.simulation.1 = function(sim.model, alternative = "two.sided", n.sim=10000, sim.seed=999, hoa.seed=999, R=10) {

  s = sim.model$s;
  x = sim.model$x;
  beta0 = sim.model$beta0;
  
  ## Simulate data
  m = n.sim;
  yy = simulate.data(m, sim.model, sim.seed);

  ## Run LR, HOA, Wald tests
  ## Dispersion unknown
  p = matrix(NA, n.sim, 3);
  colnames(p) = c("HOA", "LR", "Wald");
  msg = character(n.sim);

  set.seed(hoa.seed);

  for (i in 1:n.sim) {
    y = yy[i,];

    ## Perform LR, Wald and HOA tests
    ## args(hoa.nb.regression);
    ## debug(hoa.nb.regression);
    res = hoa.nb.regression(y, s, x, beta0, R=R,
      alternative = alternative);

    p[i,"LR"] = res$p;
    p[i,"HOA"] = res$psim;
    p[i,"Wald"] = res$p.wald;

    msg[i] = res$msg;
  }
  
  sim.pars = list(m=m, seed=sim.seed);
  hoa.pars = list(seed = hoa.seed, alternative=res$alternative, R=R);

  list (sim.model = sim.model, sim.pars = sim.pars, hoa.pars = hoa.pars, p=p, msg=msg);
}

##' Simulate NB data and perform HOA, LR, and Wald tests on simulated
##' data. The dispersion is assumed known and correct when performing
##' the tests.
##'
##' @title Simulate NB data and perform HOA, LR, and Wald tests on simulated data
##'
##' @title 
##' @param sim.model 
##' @param alternative 
##' @param n.sim 
##' @param sim.seed 
##' @param hoa.seed 
##' @return 
##' @author Yanming Di
run.simulation.2 = function(sim.model, alternative = "two.sided", n.sim=10000,
  sim.seed=999, hoa.seed=999) {

  s = sim.model$s;
  x = sim.model$x;
  beta0 = sim.model$beta0;

  ## Simulate data
  m = n.sim;
  yy = simulate.data(m, sim.model, sim.seed);

  ## Run LR, HOA, Wald tests
  ## Dispersion known
  p = matrix(NA, n.sim, 3);
  colnames(p) = c("HOA", "LR", "Wald");
  msg = character(n.sim);

  ## True dispersion value(s) will be used
  phi = 1/sim.model$kappa;

  set.seed(hoa.seed);

  for (i in 1:n.sim) {
    y = yy[i,];

    ## Perform LR, Wald and HOA tests
    res = hoa.1d(y, s, x, phi, beta0, alternative = alternative,
      print.level=0);

    p[i,"HOA"] = res$pstar;
    p[i,"LR"] = res$p;
    p[i,"Wald"] = res$p.wald;

    ## msg[i] = res$msg;
  }
  
  sim.pars = list(m=m, seed=sim.seed);
  hoa.pars = list(seed = hoa.seed, alternative=res$alternative);

  list(sim.model = sim.model, sim.pars = sim.pars, hoa.pars = hoa.pars, p=p);
}

##' Save simulation results
##'
##' @title Save simualtion results 
##' @param results a list, output from run.simulation.1
##' @return name of the saved file
##' @author Yanming Di
save.results.1 = function(results) {
  dd = Sys.Date();
  outdir = "../output";
  file.p = with(results, 
    sprintf("%s/%s.%s.seed%d.m%d.%s.pvals.phix.Rdata", outdir, dd,
            sim.model$key, sim.pars$seed, sim.pars$m,
            hoa.pars$alternative));

  with(results, save(list=objects(results),file=file.p));

  file.p
}

##' Save simulation results
##'
##' @title Save simualtion results 
##' @param results a list, output from run.simulation.2
##' @return name of the saved file
##' @author Yanming Di
save.results.2 = function(results) {
  dd = Sys.Date();
  outdir = "../output";
  file.p = with(results, 
    sprintf("%s/%s.%s.seed%d.m%d.%s.pvals.Rdata", outdir, dd,
            sim.model$key, sim.pars$seed, sim.pars$m,
            hoa.pars$alternative));

  with(results, save(list=objects(results),file=file.p));

  file.p
}

## Run simulations
sim.dr = function(sim.model, m = 10000) {

  ## alternative = "greater"
  ## dispersion unknown
  results = run.simulation.1(sim.model, "greater", m);
  file.1 = save.results.1(results);

  cat(file.1, "\n");
  cat(file.1, "\n", file="../output/file.names", append=TRUE);

  ## dispersion known
  results = run.simulation.2(sim.model, "greater", m);
  file.2 = save.results.2(results);

  cat(file.2, "\n");
  cat(file.2, "\n", file="../output/file.names", append=TRUE);

##  ## alternative = "less"
##  ## dispersion unknown
##  results = run.simulation.1(sim.model, "less", m);
##  file.3 = save.results.1(results);

##  cat(file.3, "\n");
##  cat(file.3, "\n", file="../output/file.names", append=TRUE);

##  ## dispersion known
##  results = run.simulation.2(sim.model, "less", m);
##  file.4 = save.results.2(results);
##
##  cat(file.4, "\n");
##  cat(file.4, "\n", file="../output/file.names", append=TRUE);

##  c(file.1, file.2, file.3, file.4);

  c(file.1, file.2);
}

## Fix sample size, vary fold change.
main.1 = function() {
  rm(list=ls());
  source('power.simulation.R');
  
  fc = seq(1, 2, 0.1);
  fc = 2;
  sim.model = specify.two.group.model(3, 3, 100/1e6, log(fc), 0.1);
  ## debug(simulate.data);
  ## simulate.data(10, sim.model, 999);

  ## DEBUG=TRUE;
  ## debug(hoa.nb.regression);

  for (fc in seq(1, 2, 0.1)) {
    sim.model = specify.two.group.model(3, 3, 100/1e6, log(fc), 0.1);
    ## files = sim.dr(sim.model, m=100);
    files = sim.dr(sim.model);
  }

  "Done";
}

## Fix fold change, vary sample size
main.2 = function() {
  rm(list=ls());
  source('power.simulation.R');

  ## 09/29/2013
  ## fc = 1.5;

  ## 09/29/2013
  ## fc = 1.1;

  ## 09/29/2013
  ## fc = 2.0;

  ## 09/30/2013
  fc = 1.2;

  if (FALSE) {
    ## Test runs
    
    ## n = 20;
    ## n = 15;
    ## n = 10;
    n = 5;
    sim.model = specify.two.group.model(n, n, 100/1e6, log(fc), 0.1);
    files = sim.dr(sim.model);

    ## debug(simulate.data);
    ## simulate.data(10, sim.model, 999);
  }

  for (n in 3:50) {
    sim.model = specify.two.group.model(n, n, 100/1e6, log(fc), 0.1);
    ## files = sim.dr(sim.model, m=100);
    files = sim.dr(sim.model);
  }

  if (FALSE) {
    load(files[1]);
    ## look.at.p(p);
    p2power(p);
    ## p2power(1-p);

    load(files[2]);
    ## look.at.p(p);
    p2power(p);
  }

  "Done";
}

look.at.p = function(p) {
  ## Histogram of p-values
  par(mfrow=c(2,4));
  hist(runif(sim.pars$m));
  hist(p[,"HOA"], main="HOA");
  hist(p[,"LR"], main="LR");
  hist(p[,"Wald"], main="Wald");
  ## hist(p[msg =="",2]);

  ## Compare HOA and LR p-values
  plot(p[,"HOA"], p[,"LR"],  xlab="HOA", ylab="LR", main="LR versus HOA p-values");
  abline(0, 1);

  plot(p[,"HOA"], p[,"LR"], xlab="HOA", ylab="LR", main="LR versus HOA p-values", log="xy");
  abline(0, 1);

  ## Compare HOA and Wald p-values
  plot(p[,"HOA"], p[,"Wald"],  xlab="HOA", ylab="Wald", main="Wald versus HOA p-values");
  abline(0, 1);

  plot(p[,"HOA"], p[,"Wald"], xlab="HOA", ylab="Wald", main="Wald versus HOA p-values", log="xy");
  abline(0, 1);

  ## Power calculation
  for (alpha in c(0.01, 0.05, 0.10, 0.15, 0.20)) {
    print(colMeans(p<=alpha, na.rm=TRUE));
  }

  invisible();
}

p2power = function(p, alpha = c(0.001, 0.01, 0.05, 0.10, 0.15, 0.20)) {
  n.alpha = length(alpha);
  pwr =matrix(NA, n.alpha, dim(p)[2]);
  colnames(pwr) = colnames(p);
  rownames(pwr) = alpha;
  for (i in 1:n.alpha) {
    pwr[i, ] = colMeans(p<=alpha[i], na.rm=TRUE);
  }
  pwr;
}

summarize.results = function() {

  source('power.simulation.R');
  file.names = read.table("../output/file.names", header=FALSE);

  load("../output/2013-09-19.n2.4.mu100.phi0.3.seed999.m10000.less.pvals.phix.Rdata");
  load("../output/2013-09-19.n2.4.mu100.phi0.3.seed999.m10000.less.pvals.Rdata");
  
  load(file.1);
  look.at.p(p);
  p2power(p);
  p2power(1-p);

  load(file.2);
  look.at.p(p);
  p2power(p);

  if (DEBUG) {
    ## See how much HOA p-values vary between runs
    i = 1;
    i = 4;
    i = 5;
    i = 983;
    i = 2;
    i = 9947;
    i = 7;
    y = yy[i,];

    i = 101;
    id = grep("l.hat", msg);

    id = grep("1e8", msg);

    id = grep("zsim", msg);

    id = (1:n.sim)[is.na(p[,2])];

    id = (1:n.sim)[is.na(p[,3])];

    length(id);
    p[id,];
    
    i = 4627;

    i = 2073;

    i = 2821;

    y = yy[i,];

    ## undebug(hoa.nb.regression);
    debug(hoa.nb.regression);
    res = hoa.nb.regression(y, s, x, beta0, alternative="g", R=10); 

    for (i in 1:100) {
      print(i);

      set.seed(i);
      res = hoa.nb.regression(y, s, x, beta0, alternative="g", R=10); 
    }

    res = hoa.nb.regression(y, s, x, beta0, R=10); 

    n.runs = 10;
    p.hoa =  numeric(n.runs);
    ## Large R is needed for accurate approximation.  Small R is
    ## adeqaute to get the order of magnitue right.  More runs
    ## needed to achieve accurate results for smaller p values. The
    ## approximation is more accurate when library sizes are all
    ## equal.
    R = 3000;
    R = 100;
    R = 30;

    R = 10;

    for (i in 1:n.runs) {
      res = hoa.nb.regression(y, s, x, beta0, alt="greater", R=R); 
      p.hoa[i] = res$psim;
    }
  }

  plot(p[,2], p1[,2], log="xy");
  abline(0, 1, col="magenta");

  ## Look at cases where the esitmated dispersion is  0
  id = msg != "";
  look.at.p(p[id,]);

  ## Look at cases where the esitmated dispersion is greater than 0
  look.at.p(p[!id,]);

  plot(p[id,2], p[id,1]);
  abline(0,1);

  ## DEBUG individual cases
  id.na = (1:n.sim)[is.na(p[,2])];
  i = id.na[1];

  y = yy[i,];
  debug(fit.nb.regression);
  undebug(fit.nb.regression);

  res = fit.nb.regression(y, s, x);
  res = fit.nb.regression(y, s, x, beta0);

  id = (1:n.sim)[!is.na(p[,2])];
  cor(rank(p[id,1]), rank(p[id,2]));
  
  hist(log10(p[,1]/p[,2]));
  
  par(mfrow=c(1,1));
  or = order(p[,1]);
  plot(p[or,2], p[or,1], type="b");
  plot(p[or,2], p[or,1], type="b", log="xy");
  ## lines(p[or,1], p[or,3], col="magenta");
  abline(0,1, col="cyan");

  hist(p[,4]);

  hist(rowMeans(y[,1:4]));
  hist(rowMeans(y[,2:6]));

  plot(rowMeans(y[,1:2]), rowMeans(y[,3:6]), log="xy");
}
