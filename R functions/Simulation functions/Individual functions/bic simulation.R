
## Simulate 3-component data
## Specifying model parameters used to generate data
sim.par3 = function(N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p){
  out = list(N=N,tau1=tau1,tau2=tau2,mu1=mu1,mu2=mu2,mu3=mu3,
    sig1=sig1,sig2=sig2,sig3=sig3,k=k,p=p)
  return(out)
}


## Generate data using given parameters
sim.data3 = function(simpar){  
  N = simpar$N; tau1 = simpar$tau1; tau2 = simpar$tau2
  mu1 = simpar$mu1; mu2 = simpar$mu2; mu3 = simpar$mu3
  sig1 = simpar$sig1; sig2 = simpar$sig2; sig3 = simpar$sig3
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
  z.ini = matrix(0,N,2)
  rand.samples = matrix(0,N,2)
  for(i in 1:N){
    if(U[i]<tau1){
      rand.samples[i,] = MASS::mvrnorm(1,mu1,(sig1+errmat[,,i]))
      z.ini[i,] = c(1,0)
    } else if(U[i]<(tau1+tau2)){
      rand.samples[i,] = MASS::mvrnorm(1,mu2,(sig2+errmat[,,i]))
      z.ini[i,] = c(0,1)
    } else{
      rand.samples[i,] = MASS::mvrnorm(1,mu3,(sig3+errmat[,,i]))
    }
  }

  out = list(data=rand.samples, z.ini=z.ini, err=errmat, index=index, k=k)
  return(out)
  
}




## Implement selection of number of clusters with BIC ##

nc.select = function(data, err, nc=1:5, d=1, itmax=Inf, lb=1e-3){
  ## "maxClnum" lets user specify the maximum number of clusters
  ## in consideration.
  
  G = Gmin = min(nc)
  Gmax = max(nc)
  bic.vec = numeric(length(nc))
  N = nrow(data)
  hcTree = mclust::hc(data)
  k = 1
  res = list()
  
  while(G <= Gmax){
    print(paste("Components = ", G, sep=""))
    
    # Generate initial membership matrix with hierarchical agglomeration
    cl = mclust::hclass(hcTree, G)
    z.ini = matrix(0,N,G)
    for(i in 1:N){
      for(j in 1:G){
        z.ini[i,j] = ifelse(cl[i]==j,1,0)
      }
    }

    # Run MCLUST-ME with chosen initial membership
    out = mcmeVVV(data, z.ini, err, d, itmax, lb)
    res = c(res,out)
    bic.vec[k] = out$BIC
    
    G = G + 1
    k = k + 1
  }
  
  ## Allow user to view results of the optimal model chosen
  n.opt = which(bic.vec==max(bic.vec))
  l = 14*n.opt - 13
  u = 14*n.opt
  res.opt = res[l:u]
  
  ## Output BIC values for each choice of number of clusters
  tbl = cbind(nc,bic.vec)

  out = list(G=nc[n.opt], bic=max(bic.vec), res.opt=res.opt, BIC=tbl,
    nc=nc, bicvec=bic.vec, res=res)
  return(out)
}  
 

## Plot BIC for each choice of number of clusters
plotbic = function(nc.result,title){
  nc = nc.result$nc
  bic.vec = nc.result$bicvec
  
  # get the range for the x and y axis
  xrange = range(nc)
  yrange = range(bic.vec)
  # set up the plot
  plot(xrange, yrange, type="n", xlab="number of components",
     ylab="BIC", xaxt="n", main=title)
  # add lines
  lines(nc, bic.vec, type="b", lwd=1.5, col="red")
  axis(1, at=nc)
}  
  


## Examine difference in model selection results by BIC, with 3 components
test.3group = function(seed,N,tau1,tau2,mu1,mu2,mu3,sig1,sig2,sig3,k,p){
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

