

sim.driver = function(N,tau,mu1,mu2,sig1,sig2,k,p,nseed){
  
  source('~/Research Project/Simulations/Sim/Functions/par_data_run.R')
  source('~/Research Project/R functions for MCME(VVV)/Other functions/fuzzyrand.R')
  source('~/Research Project/Simulations/Sim/Functions/rand.R')
  source('~/Research Project/Simulations/Sim/Functions/boundary.R')
  
  simres = list()
  randres = matrix(,nrow=nseed,ncol=8)
  
  v = seq(1,2*nseed,by=1)
  seedvec = sample(v,size=nseed)
  
  for(i in 1:nseed){
    
    print(paste("seed=",seedvec[i],sep=""))
    
    # obtain and store clustering results    
    data.par = sim.par(N,tau,mu1,mu2,sig1,sig2,k,p)
    set.seed(seedvec[i])
    simdata = sim.data(data.par)
    simrun = sim.run(simdata)
    simres = c(simres,simrun)
    
    # obtain and store rand indices
    rindex = get.rand(simrun)
    randres[i,] = rindex
    
    # create and store boundary plots
    file = paste("bdry",seedvec[i],".png",sep="")
    png(filename=file,1500,600)
    plot.boundary(simrun,seedvec[i])
    dev.off()
  }
  
  # calculate average rand index
  rand.avg = colMeans(randres)
  
  out = list(sim.result=simres, rand.raw=randres, rand.avg=rand.avg, seed=seedvec)
  save(out,file="res.RData")
}

