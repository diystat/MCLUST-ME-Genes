
###### Calculate lower and upper bound for Rand index for MCLUST
###### when no measurement errors are present

## Obtain results from MCME and meVVV
rand.mevvv = function(simdata){
  rand.samples = simdata$data
  z.ini = simdata$z.ini
  errmat = simdata$err
  index = simdata$index
  k = simdata$k
  
  ## Run mevvv:
  res.mevvv = mclust::meVVV(rand.samples,z.ini)
  z.mevvv = res.mevvv$z
  
  N = nrow(z.ini)
  G = ncol(z.ini)
  
  # Discretize membership matrix:
  z1 = matrix(0,N,G)
  for(i in 1:N){
    rm.mevvv = max(z.mevvv[i,])
    for(k in 1:G){
      z1[i,k] = ifelse(z.mevvv[i,k]==rm.mevvv,1,0)
    }  
  }
  
  z.ini1 = z.ini[,1] + 1
  z1 = z1[,1] + 1

  ## Calculate Rand index:
  rand.mevvv = round(phyclust::RRand(z.ini1,z1)$Rand,4)
  return(rand.mevvv)
}


sim.driver.lu = function(N,tau,mu1,mu2,sig1,sig2,k,p,nseed){

  simres = list()
  randres = matrix(0,nrow=nseed,ncol=8)
  rr = numeric(nseed)
  v = seq(1,2*nseed,by=1)
  seedvec = sample(v,size=nseed)
  
  for(i in 1:nseed){
    
    print(paste("seed=",seedvec[i],sep=""))
    
    # obtain and store clustering results    
    data.par = sim.par(N,tau,mu1,mu2,sig1,sig2,k,p)
    set.seed(seedvec[i])
    simdata = sim.data(data.par)
    
    rr[i] = rand.mevvv(simdata)
    
  }
  
  lu = c(min(rr),max(rr))
  return(lu)
}



randLU = function(p,nseed){
  N = 200
  tau = 0.5
  mu1 = c(0,0)
  mu2 = c(8,0)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  k = 0
  #nseed = 1000

  # Run the simulation
  out = sim.driver.lu(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
  return(out)
}        

randLU(0.5,1000) # 0.5156 0.7954

    