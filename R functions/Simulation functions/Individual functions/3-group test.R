
## Examine difference in model selection results by BIC, with 3 components
test.3group = function(seed){

  N = 300
  tau1 = 0.3
  tau2 = 0.4
  mu1 = c(-16,0)
  mu2 = c(16,0)
  mu3 = c(0,24)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  sig3 = matrix(c(36,0,0,36),nrow=2)
  k = 6
  p = 0.5

  set.seed(seed)
  data.par = sim.par3(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p)
  simdata = sim.data3(data.par)

  dat = simdata$data
  err = simdata$err

  out.mcme = nc.select(dat, err, 1:7)
  out.mclust = Mclust(dat, G=1:7, "VVV")
  
  par(mfrow=c(1,2))
  plotbic(out.mcme,"MCLUST-ME")
  plot(out.mclust,"BIC",legendArgs=list(plot=FALSE),main="MCLUST")
  par(mfrow=c(1,1))
  
  out = list(out.mcme = out.mcme$BIC, out.mclust = out.mclust$BIC)
  
}

