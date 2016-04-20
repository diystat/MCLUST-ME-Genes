
## Check if meVVV has singularity issue for any of the random seeds
check.seed = function(N,tau,mu1,mu2,sig1,sig2,k,p){
  seedvec = seq(1,200,1)
  for(i in 1:200){
    data.par = sim.par(N,tau,mu1,mu2,sig1,sig2,k,p)
    print(seedvec[i])
    set.seed(seedvec[i])
    simdata = sim.data(data.par)
  
    rand.samples = simdata$data
    z.ini = simdata$z.ini
  
    res = meVVV(rand.samples,z.ini)
  }
}

