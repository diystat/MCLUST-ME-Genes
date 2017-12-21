

# Driver function for initialization simulation
initialSim.driver = function(N,tau,mu1,mu2,sig1,sig2,k,p,nseed){

  res = matrix(0,nrow=nseed,ncol=4)
  colnames(res) = c("rt.true","rt.mclust","adjRand.true","adjRand.mclust")
  #colnames(res) = c("rt.true","rt.mclust","rt.rnd.avg",
  #  "adjRand.true","adjRand.mclust","adjRand.rnd.avg")
  
  v = seq(1,2*nseed,by=1)
  seedvec = sample(v,size=nseed)
  
  for(i in 1:nseed){
    
    print(paste("seed=",seedvec[i],sep=""))
    
    # obtain and store clustering results    
    data.par = sim.par(N,tau,mu1,mu2,sig1,sig2,k,p)
    set.seed(seedvec[i])
    simdata = sim.data(data.par)
    res[i,] = initialSim.run(simdata)
  }
  
  # calculate average rand index
  res.avg = colMeans(res)
  print(res.avg)
  
  out = list(initialSim.res = res, avg.res = res.avg)
  return(out)
}



## Obtain results from MCME and meVVV
initialSim.run = function(simdata){
  rand.samples = simdata$data
  errmat = simdata$err
  index = simdata$index
  k = simdata$k
  
  ## True membership assignment
  z.ini.true = simdata$z.ini.true

  ## Step 1: Run mevvv and obtain membership
  res.mevvv = mclust::meVVV(rand.samples,z.ini.true)
  z.ini.mclust = res.mevvv$z
  
  ## Step 2: Generate 10 random initial membership assignments
  z.ini.r = function(simdata){
    rand.samples = simdata$data
    N = nrow(rand.samples)
    z.ini.r = matrix(0,N,2)
    p = runif(1,0.1,0.9)
    ind = rbinom(N,1,p)+1
    for(i in 1:N){
      j = ind[i]
      z.ini.r[i,j] = 1
    }
    return(z.ini.r)
  }
  
  z.ini.rnd = array(0,dim=c(N,2,10))
  for(i in 1:10){
    z.ini.rnd[,,i] = z.ini.r(simdata)
  }
  
  ## Step 3: Initialize with true membership
  start = proc.time()
  res.true = mcmeVVV(rand.samples, z.ini.true, errmat)
  rt = proc.time()-start
  rt.true = round(rt[3],4)
  
  memprob.true = res.true$z
  adjRand.true = getAdjRand(simdata,memprob.true)
  print(paste("True membership: time=",rt.true,", adjRand=",adjRand.true,sep=""))

  
  ## Step 4: Initialize with MCLUST membership
  start = proc.time()
  res.mclust = mcmeVVV(rand.samples, z.ini.mclust, errmat)
  rt = proc.time()-start
  rt.mclust = round(rt[3],4)
  
  memprob.mclust = res.mclust$z
  adjRand.mclust = getAdjRand(simdata,memprob.mclust)
  print(paste("MCLUST: time=",rt.mclust,", adjRand=",adjRand.mclust,sep=""))


  ## Step 5: Initialize with 10 random memberships
  #rt.rnd = adjRand.rnd = numeric(10)
  #for(i in 1:10){
  #  z.ini = z.ini.rnd[,,i]
  #  start = proc.time()
  #  res.rnd = mcmeVVV(rand.samples, z.ini, errmat)
  #  rt = proc.time()-start
  #  rt.rnd[i] = round(rt[3],4)
    
  #  memprob.rnd = res.rnd$z
  #  adjRand.rnd[i] = getAdjRand(simdata,memprob.rnd)
  #}
  #rt.rnd.avg = mean(rt.rnd)
  #adjRand.rnd.avg = mean(adjRand.rnd)
  #print(paste("Random: time=",rt.rnd.avg,", adjRand=",adjRand.rnd.avg,sep=""))
  
  
  #out = c(rt.true=rt.true,rt.mclust=rt.mclust,rt.rnd.avg=rt.rnd.avg,
  #  adjRand.true=adjRand.true,adjRand.mclust=adjRand.mclust,
  #  adjRand.rnd.avg=adjRand.rnd.avg)
  
  out = c(rt.true=rt.true,rt.mclust=rt.mclust,
    adjRand.true=adjRand.true,adjRand.mclust=adjRand.mclust)
  return(out)
}




# Calculate adjusted Rand index
getAdjRand = function(simdata,memprob){
  # true membership
  truth = simdata$z.ini.true
  N = nrow(simdata$data)

  # Discretize membership matrix:
  clusters = matrix(0,N,2)
  for(i in 1:N){
    maxprob = max(memprob[i,])
    for(k in 1:2){
      clusters[i,k] = ifelse(memprob[i,k]==maxprob,1,0)
    }  
  }
  truth_1 = truth[,1] + 1
  clusters_1 = clusters[,1] + 1
  tmp = phyclust::RRand(truth_1,clusters_1)
  out = round(tmp$adjRand,4)
  return(out)
}



## Specifying model parameters used to generate data
sim.par = function(N,tau,mu1,mu2,sig1,sig2,k,p){
  out = list(N=N,tau=tau,mu1=mu1,mu2=mu2,sig1=sig1,sig2=sig2,k=k,p=p)
  return(out)
}



## Generate data using given parameters
sim.data = function(simpar){  
  N = simpar$N; tau = simpar$tau
  mu1 = simpar$mu1; mu2 = simpar$mu2
  sig1 = simpar$sig1; sig2 = simpar$sig2
  k = simpar$k; p = simpar$p  
  E = k*matrix(c(1,0,0,1),nrow=2)
    
  ## Randomly assign constant measurement error to 100p% of observations:
  index = rbinom(N,1,p)  
  errmat = array(0,dim=c(2,2,N))
  for(i in 1:N){
    errmat[,,i] = index[i]*E
  }
 
  ## Generate sample from mixture distribution:
  U = runif(N)
  z.ini.true = matrix(0,N,2)
  rand.samples = matrix(0,N,2)
  for(i in 1:N){
    if(U[i]<tau){
      rand.samples[i,] = MASS::mvrnorm(1,mu1,(sig1+errmat[,,i]))
      z.ini.true[i,] = c(1,0)
    } else{
      rand.samples[i,] = MASS::mvrnorm(1,mu2,(sig2+errmat[,,i]))
      z.ini.true[i,] = c(0,1)
    }
  }

  out = list(data=rand.samples, z.ini.true=z.ini.true, err=errmat, index=index, k=k)
  return(out)
}



## Check if mcmeVVV has singularity issue for any of the random seeds
check.seed = function(N,tau,mu1,mu2,sig1,sig2,k,p){
  seedvec = seq(1,200,1)
  for(i in 1:200){
    data.par = sim.par(N,tau,mu1,mu2,sig1,sig2,k,p)
    print(seedvec[i])
    set.seed(seedvec[i])
    simdata = sim.data(data.par)
  
    rand.samples = simdata$data
    z.ini.true = simdata$z.ini.true
  
    res = meVVV(rand.samples,z.ini.true)
  }
}



############################################################################

############################################################################



  N = 200
  tau = 0.5
  mu1 = c(0,0)
  mu2 = c(8,0)
  sig1 = matrix(c(64,0,0,64),nrow=2)
  sig2 = matrix(c(16,0,0,16),nrow=2)
  k = 36
  nseed = 100

  library(mclust)
  
  # Check if singularity occurs for any seed
  check.seed(N,tau,mu1,mu2,sig1,sig2,k,p)

  # Run the simulation
  out = initialSim.driver(N,tau,mu1,mu2,sig1,sig2,k,p,nseed)
  saveRDS(out,"iniSim_res_noRND.rds")


rt = out$initialSim.res[,1:2]
plot(rt[,2]-rt[,1],ylab="mclust-true",main="Pairwise runtime difference")
abline(0,0,col="red",lty="dashed")
# Test if mclust initialization has shorter runtime
t.test(rt[,1],rt[,2],alternative="greater",paired=TRUE) # p-value = 0.05247

adjR = out$initialSim.res[,3:4]
plot(adjR[,2]-adjR[,1],ylab="mclust-true",main="Pairwise adj Rand difference")
abline(0,0,col="red",lty="dashed")
# Test if mclust initialization leads to less accurate classification
t.test(adjR[,1],adjR[,2],alternative="greater",paired=TRUE) # p-value = 0.0002185


par(mfrow=c(1,2))
plot(rt[,2]-rt[,1],ylab="mclust-true",main="Pairwise runtime difference")
abline(0,0,col="red",lty="dashed")

plot(adjR[,2]-adjR[,1],ylab="mclust-true",main="Pairwise adj Rand difference")
abline(0,0,col="red",lty="dashed")
par(mfrow=c(1,1))



