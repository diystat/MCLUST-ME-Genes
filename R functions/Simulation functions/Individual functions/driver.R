
# Simulation driver function for Sim 1 and Sim 2
sim.driver = function(N,tau,mu1,mu2,sig1,sig2,k,p,nseed){
 
  simres = list()
  randres = matrix(0,nrow=nseed,ncol=8)
  
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
  return(out)
}

